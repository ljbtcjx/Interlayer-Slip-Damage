"""Reproducible preparation record. Refuses to overwrite generated inputs."""
from pathlib import Path
import csv
import hashlib
import json
import math
import shutil

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SMOOTH = ROOT / "平滑面拉伸 - 提交"
EDGE = ROOT / "边缘面拉伸 - 提交"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, text):
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite: {path}")
    path.write_text(text, encoding="utf-8", newline="\n")


lines = (SMOOTH / "initial_vdm.xyz").read_text(encoding="utf-8-sig").splitlines()
atoms = [[p[0], *map(float, p[1:4])] for p in (line.split() for line in lines[2:] if line.strip())]
assert len(atoms) == int(lines[0]) == 242
sheet = list(range(97, 147))
assert all(atoms[i-1][0] == "C" and 24 < atoms[i-1][3] < 24.2 for i in sheet)
ymax = max(atoms[i-1][2] for i in sheet)
grip = [i for i in sheet if atoms[i-1][2] >= ymax - 1.0]
bottom = [i+1 for i, a in enumerate(atoms) if a[3] < 8.5]
assert len(grip) == 10 and len(bottom) == 36 and not set(grip) & set(bottom)
mins = [min(a[j+1] for a in atoms) for j in range(3)]
maxs = [max(a[j+1] for a in atoms) for j in range(3)]
shift = [12 - x for x in mins]
cell = [math.ceil(maxs[j] - mins[j] + 24 + (4 if j == 2 else 0)) for j in range(3)]
prepared = [[a[0], *[a[j+1] + shift[j] for j in range(3)]] for a in atoms]
max_force = 0.02 / (27.211386245988 / 0.529177210903)

# Back up the old smooth input before replacing the copied, invalid edge template.
backup = SMOOTH / "原文件备份"
backup.mkdir(exist_ok=True)
old = SMOOTH / "step_1.inp"
if (SMOOTH / "loading_protocol.json").exists():
    raise RuntimeError("Inputs already generated. Do not rerun in a calculation directory.")
if old.exists():
    target = backup / "step_1.inp"
    if target.exists():
        raise FileExistsError(target)
    shutil.copy2(old, target)

manifest = {"release": "CP2K v2025.2", "download_date": "2026-10-20", "files": {}}
for name in ("BASIS_MOLOPT", "GTH_POTENTIALS", "dftd3.dat"):
    path = HERE / name
    assert path.stat().st_size > 10000
    manifest["files"][name] = {"url": f"https://raw.githubusercontent.com/cp2k/cp2k/v2025.2/data/{name}", "bytes": path.stat().st_size, "sha256": sha(path)}
manifest["files"]["POTENTIAL"] = {**manifest["files"]["GTH_POTENTIALS"], "note": "Byte-identical alias for legacy edge inputs"}
for folder in (SMOOTH, EDGE):
    for name in ("BASIS_MOLOPT", "GTH_POTENTIALS", "dftd3.dat", "POTENTIAL"):
        source = HERE / ("GTH_POTENTIALS" if name == "POTENTIAL" else name)
        dest = folder / name
        if dest.exists() and sha(dest) != sha(source):
            raise RuntimeError(f"Existing different data file: {dest}")
        shutil.copy2(source, dest)
    save(folder / "data_sources.json", json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
save(HERE / "data_sources.json", json.dumps(manifest, indent=2, ensure_ascii=False) + "\n")
save(SMOOTH / "step_0.xyz", f"{len(prepared)}\nOriginal geometry translated into nonperiodic box; no geometry relaxation yet\n" + "".join(f"{a[0]:2s} {a[1]:.10f} {a[2]:.10f} {a[3]:.10f}\n" for a in prepared))
cfg = {
    "original_file": "initial_vdm.xyz", "original_sha256": sha(SMOOTH / "initial_vdm.xyz"),
    "prepared_file": "step_0.xyz", "prepared_sha256": sha(SMOOTH / "step_0.xyz"),
    "atom_numbering": "1-based XYZ row order; no atom reordering",
    "sheet_ids": sheet, "grip_ids": grip, "bottom_ids": bottom,
    "grip_selection": "Top-Y strip of carbon sheet (width 1.0 angstrom); fixed XYZ within relaxation",
    "bottom_selection": "Original z < 8.5 angstrom; fixed XYZ throughout",
    "load_axis": "+Z", "increment_angstrom": 0.2, "steps": 20,
    "total_grip_displacement_angstrom": 4.0,
    "coordinate_translation_angstrom": shift, "cell_angstrom": cell,
    "periodic": "NONE", "poisson_solver": "WAVELET",
    "max_force_hartree_per_bohr": max_force,
    "note": "step_0 is a separate zero-displacement relaxation, not one of the 20 loading increments"
}
save(SMOOTH / "loading_protocol.json", json.dumps(cfg, indent=2, ensure_ascii=False) + "\n")
with (SMOOTH / "atom_groups.csv").open("w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["xyz_row_id", "element", "group", "original_x", "original_y", "original_z"])
    for i, a in enumerate(atoms, 1):
        group = "moving_grip" if i in grip else "fixed_bottom" if i in bottom else "free_carbon_sheet" if i in sheet else "free_substrate"
        w.writerow([i, a[0], group, *a[1:]])


def inp(n, single=False):
    project = f"Smooth_{'SP' if single else 'Relax'}_{n:02d}"
    coordinates = f"relaxed_{n:02d}.xyz" if single else f"step_{n}.xyz"
    prev = n if single else n - 1
    restart = f"    WFN_RESTART_FILE_NAME Smooth_Relax_{prev:02d}-RESTART.wfn\n" if prev >= 0 else ""
    output = ""
    if single:
        output = """    &PRINT
      &HIRSHFELD ON
      &END HIRSHFELD
      &MULLIKEN ON
      &END MULLIKEN
      &MOMENTS ON
        MAX_MOMENT 1
      &END MOMENTS
      &E_DENSITY_CUBE ON
        STRIDE 2 2 2
      &END E_DENSITY_CUBE
      &ELF_CUBE ON
        STRIDE 2 2 2
      &END ELF_CUBE
    &END PRINT
"""
    motion = ""
    if not single:
        lists = "\n".join("      LIST " + " ".join(map(str, group[j:j+15])) for group in (bottom, grip) for j in range(0, len(group), 15))
        motion = f"""
&MOTION
  &GEO_OPT
    TYPE MINIMIZATION
    OPTIMIZER BFGS
    MAX_ITER 300
    ! CP2K force components in Hartree/bohr; equivalent to 0.02 eV/angstrom.
    MAX_FORCE {max_force:.12E}
    RMS_FORCE 2.0E-4
    MAX_DR 3.0E-3
    RMS_DR 1.5E-3
  &END GEO_OPT
  &CONSTRAINT
    &FIXED_ATOMS
      COMPONENTS_TO_FIX XYZ
{lists}
    &END FIXED_ATOMS
  &END CONSTRAINT
  &PRINT
    &TRAJECTORY ON
      FORMAT XYZ
      FILENAME =trajectory_{n:02d}.xyz
      ADD_LAST NUMERIC
      &EACH
        GEO_OPT 1
      &END EACH
    &END TRAJECTORY
  &END PRINT
&END MOTION
"""
    return f"""! Finite smooth interface; grip +Z displacement = {n * .2:.1f} angstrom.
! Coordinates for loading steps depend on preceding converged relaxation.
&GLOBAL
  PROJECT {project}
  RUN_TYPE {'ENERGY_FORCE' if single else 'GEO_OPT'}
  PRINT_LEVEL MEDIUM
  PREFERRED_FFT_LIBRARY FFTSG
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
{restart}    CHARGE 0
    UKS TRUE
    MULTIPLICITY 1
    &QS
      EPS_DEFAULT 1.0E-12
    &END QS
    &MGRID
      CUTOFF [Ry] 400
      REL_CUTOFF [Ry] 60
    &END MGRID
    &POISSON
      PERIODIC NONE
      POISSON_SOLVER WAVELET
    &END POISSON
    &SCF
      MAX_SCF 200
      EPS_SCF 1.0E-6
      SCF_GUESS {'RESTART' if prev >= 0 else 'ATOMIC'}
      IGNORE_CONVERGENCE_FAILURE FALSE
      ADDED_MOS 100
      &DIAGONALIZATION
        ALGORITHM STANDARD
      &END DIAGONALIZATION
      &SMEAR
        METHOD FERMI_DIRAC
        ! Numerical electronic smearing; this is not a 300 K ionic MD run.
        ELECTRONIC_TEMPERATURE [K] 300
      &END SMEAR
      &MIXING
        METHOD BROYDEN_MIXING
        ALPHA 0.1
        NBROYDEN 8
      &END MIXING
      &OUTER_SCF
        MAX_SCF 10
        EPS_SCF 1.0E-6
      &END OUTER_SCF
      &PRINT
        &RESTART ON
        &END RESTART
      &END PRINT
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
      &VDW_POTENTIAL
        POTENTIAL_TYPE PAIR_POTENTIAL
        &PAIR_POTENTIAL
          TYPE DFTD3(BJ)
          PARAMETER_FILE_NAME dftd3.dat
          REFERENCE_FUNCTIONAL PBE
        &END PAIR_POTENTIAL
      &END VDW_POTENTIAL
    &END XC
{output}  &END DFT
  &SUBSYS
    &CELL
      ABC [angstrom] {cell[0]} {cell[1]} {cell[2]}
      PERIODIC NONE
    &END CELL
    &TOPOLOGY
      COORD_FILE_FORMAT XYZ
      COORD_FILE_NAME {coordinates}
    &END TOPOLOGY
    &KIND C
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
    &KIND Si
      BASIS_SET DZVP-MOLOPT-SR-GTH
      POTENTIAL GTH-PBE-q4
    &END KIND
  &END SUBSYS
&END FORCE_EVAL
{motion}"""


for n in range(21):
    for single in (False, True):
        path = SMOOTH / f"{'sp' if single else 'step'}_{n}.inp"
        if n == 1 and not single:
            path.write_text(inp(n, single), encoding="utf-8", newline="\n")
        else:
            save(path, inp(n, single))
print(json.dumps({"cell_angstrom": cell, "grip_ids": grip, "bottom_count": len(bottom), "input_files": 42}, ensure_ascii=False))
