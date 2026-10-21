#!/usr/bin/env python3
"""Sequential CP2K loading/relaxation/single-point jobs; Python 3, no packages.

Default is a dry run. Use --run to launch CP2K in the active cluster allocation.
Each completed stage is fingerprinted. Failed stages are never advanced.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_xyz(path):
    """Return the last complete frame; reject truncated or malformed output."""
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    offset, last = 0, None
    while offset < len(lines):
        if not lines[offset].strip():
            offset += 1
            continue
        n = int(lines[offset])
        if n < 1 or offset + n + 2 > len(lines):
            raise ValueError(f"Incomplete XYZ frame: {path}")
        last = []
        for line in lines[offset + 2:offset + n + 2]:
            parts = line.split()
            if len(parts) < 4 or parts[0] not in ("C", "Si"):
                raise ValueError(f"Bad XYZ atom: {line}")
            xyz = [float(x) for x in parts[1:4]]
            if not all(math.isfinite(x) for x in xyz):
                raise ValueError("Non-finite coordinate")
            last.append([parts[0], *xyz])
        offset += n + 2
    if last is None:
        raise ValueError(f"Empty XYZ file: {path}")
    return last


def write_xyz(path, atoms, comment):
    content = f"{len(atoms)}\n{comment}\n" + "".join(
        f"{a[0]:2s} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in atoms)
    if path.exists() and path.read_text(encoding="utf-8") != content:
        raise RuntimeError(f"Refusing to replace different coordinates: {path.name}")
    path.write_text(content, encoding="utf-8", newline="\n")


def validate_atoms(atoms, cfg, reference, step):
    if len(atoms) != len(reference) or [a[0] for a in atoms] != [a[0] for a in reference]:
        raise ValueError("Atom count/order differs from initial structure")
    for atom in atoms:
        for j in range(3):
            if not 5.0 < atom[j + 1] < cfg["cell_angstrom"][j] - 5.0:
                raise ValueError("Atom approaches box boundary (<5 A vacuum). Review cell size.")
    moving = set(cfg["grip_ids"])
    for atom_id in cfg["bottom_ids"] + cfg["grip_ids"]:
        target = reference[atom_id - 1][1:].copy()
        if atom_id in moving:
            target[2] += step * cfg["increment_angstrom"]
        if max(abs(atoms[atom_id - 1][j + 1] - target[j]) for j in range(3)) > 2e-5:
            raise ValueError(f"Fixed boundary moved incorrectly: atom {atom_id}, step {step}")


def load_next(previous, cfg):
    atoms = [a.copy() for a in previous]
    for atom_id in cfg["grip_ids"]:
        atoms[atom_id - 1][3] += cfg["increment_angstrom"]
    return atoms


def check_log(path, geometry=False):
    text = path.read_text(encoding="utf-8", errors="replace")
    if re.search(r"SCF run\s+NOT converged|SCF run did not converge|ABORT|GEOMETRY OPTIMIZATION NOT CONVERGED|MAXIMUM NUMBER OF OPTIMIZATION STEPS REACHED", text, re.I):
        raise RuntimeError(f"Failed/nonconverged calculation: {path.name}")
    if "PROGRAM ENDED AT" not in text or not re.search(r"SCF run converged", text, re.I):
        raise RuntimeError(f"Missing normal completion/SCF convergence: {path.name}")
    if geometry and "GEOMETRY OPTIMIZATION COMPLETED" not in text:
        raise RuntimeError(f"Missing geometry convergence: {path.name}")
    values = re.findall(r"ENERGY\|\s+Total FORCE_EVAL.*?energy.*?([-+]?\d+\.\d+(?:[EeDd][-+]?\d+)?)\s*$", text, re.M)
    if not values:
        raise RuntimeError(f"No ENERGY| result: {path.name}")
    return float(values[-1].replace("D", "E").replace("d", "e"))


def stage(name, dependencies, outputs, command, geometry):
    """Resume only a stage with a valid receipt and unchanged input/output files."""
    inp, out = ROOT / f"{name}.inp", ROOT / f"{name}.out"
    receipt = ROOT / f"{name}.done.json"
    deps = [inp, ROOT / "run_quasistatic.py", ROOT / "loading_protocol.json", *dependencies]
    fingerprints = {p.name: digest(p) for p in deps}
    if receipt.exists():
        old = json.loads(receipt.read_text(encoding="utf-8"))
        if old["inputs"] != fingerprints or old["outputs"] != {p.name: digest(p) for p in [out, *outputs]}:
            raise RuntimeError(f"Files changed since completed stage {name}; use a separate run folder.")
        check_log(out, geometry)
        print(f"Verified completed stage: {name}", flush=True)
        return
    # Interrupted/failed output is preserved, not silently reused or overwritten.
    occupied = [p.name for p in [out, ROOT / f"{name}.launcher.log", *outputs] if p.exists()]
    if occupied:
        raise RuntimeError(f"Unverified outputs for {name}: {occupied}. Archive this failed stage's outputs before retrying.")
    print(f"Running {name} ...", flush=True)
    with (ROOT / f"{name}.launcher.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(command + ["-i", inp.name, "-o", out.name], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        raise RuntimeError(f"CP2K exit code {result.returncode}: see {name}.launcher.log / {name}.out")
    check_log(out, geometry)
    for p in outputs:
        if not p.is_file() or p.stat().st_size == 0:
            raise RuntimeError(f"Missing CP2K output: {p.name}")
    receipt.write_text(json.dumps({"inputs": fingerprints, "outputs": {p.name: digest(p) for p in [out, *outputs]}}, indent=2), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="store_true", help="Launch calculations (otherwise check files and show plan only)")
    parser.add_argument("--end", type=int, default=20, help="Last loading step (0..20); 0 is baseline only")
    parser.add_argument("--cp2k-command", default="cp2k.psmp", help='E.g. "mpirun -np 30 cp2k.psmp"; use allocation-specific launcher')
    args = parser.parse_args()
    cfg = json.loads((ROOT / "loading_protocol.json").read_text(encoding="utf-8"))
    if not 0 <= args.end <= cfg["steps"]:
        parser.error("--end must be between 0 and 20")
    ref = read_xyz(ROOT / "step_0.xyz")
    validate_atoms(ref, cfg, ref, 0)
    if digest(ROOT / "initial_vdm.xyz") != cfg["original_sha256"] or digest(ROOT / "step_0.xyz") != cfg["prepared_sha256"]:
        raise RuntimeError("Initial structure changed; regenerate inputs and review atom groups")
    data = [ROOT / s for s in ("BASIS_MOLOPT", "GTH_POTENTIALS", "dftd3.dat")]
    manifest = json.loads((ROOT / "data_sources.json").read_text(encoding="utf-8"))
    for p in data:
        if digest(p) != manifest["files"][p.name]["sha256"]:
            raise RuntimeError(f"Data-file hash mismatch: {p.name}")
    for n in range(args.end + 1):
        for name in (f"step_{n}.inp", f"sp_{n}.inp"):
            if not (ROOT / name).is_file():
                raise FileNotFoundError(name)
    print(f"{len(ref)} atoms; grip={len(cfg['grip_ids'])}; bottom={len(cfg['bottom_ids'])}; free={len(ref)-len(cfg['grip_ids'])-len(cfg['bottom_ids'])}")
    print(f"Baseline relaxation + SP, then {args.end} x 0.2 A load/relax/SP; endpoint {args.end * .2:.1f} A")
    if not args.run:
        print("DRY RUN ONLY. Later step_N.xyz files are created from the preceding converged geometry, not from rigid snapshots.")
        return
    command = shlex.split(args.cp2k_command)
    if not command:
        parser.error("Empty CP2K command")
    rows = []
    for n in range(args.end + 1):
        source = ROOT / f"step_{n}.xyz"
        if n:
            previous = read_xyz(ROOT / f"relaxed_{n-1:02d}.xyz")
            validate_atoms(previous, cfg, ref, n - 1)
            write_xyz(source, load_next(previous, cfg), f"Step {n}; grip displacement {n * .2:.1f} A; inherited relaxed interior")
        validate_atoms(read_xyz(source), cfg, ref, n)
        project = f"Smooth_Relax_{n:02d}"
        wfn = ROOT / f"{project}-RESTART.wfn"
        trajectory = ROOT / f"trajectory_{n:02d}.xyz"
        inherited = [ROOT / f"Smooth_Relax_{n-1:02d}-RESTART.wfn"] if n else []
        stage(f"step_{n}", [*data, source, *inherited], [trajectory, wfn], command, True)
        atoms = read_xyz(trajectory)
        validate_atoms(atoms, cfg, ref, n)
        relaxed = ROOT / f"relaxed_{n:02d}.xyz"
        write_xyz(relaxed, atoms, f"Converged geometry after prescribed grip displacement {n * .2:.1f} A")
        stage(f"sp_{n}", [*data, relaxed, wfn], [], command, False)
        rows.append([n, f"{n * .2:.1f}", check_log(ROOT / f"sp_{n}.out")])
        with (ROOT / "energy_vs_displacement.csv").open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["step", "grip_displacement_angstrom", "cp2k_energy_hartree"])
            writer.writerows(rows)
    print("Requested load/relax/SP sequence completed.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError, KeyError) as exc:
        print(f"STOP: {exc}", file=sys.stderr)
        sys.exit(1)
