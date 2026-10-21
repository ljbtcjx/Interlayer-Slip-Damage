"""Static input checks and mock-process workflow tests; NOT a DFT calculation."""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import re
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "平滑面拉伸 - 提交"
spec = importlib.util.spec_from_file_location("runner", SOURCE / "run_quasistatic.py")
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
cfg = json.loads((SOURCE / "loading_protocol.json").read_text(encoding="utf-8"))
original = runner.read_xyz(SOURCE / "initial_vdm.xyz")
ref = runner.read_xyz(SOURCE / "step_0.xyz")
assert len(ref) == 242 and len(cfg["grip_ids"]) == 10 and len(cfg["bottom_ids"]) == 36
assert not set(cfg["grip_ids"]) & set(cfg["bottom_ids"])
assert all(ref[i-1][0] == "C" for i in cfg["grip_ids"])
for a, b in zip(original, ref):
    assert a[0] == b[0]
    assert all(abs(b[j+1] - a[j+1] - cfg["coordinate_translation_angstrom"][j]) < 1e-9 for j in range(3))
for n in range(21):
    for prefix in ("step", "sp"):
        text = (SOURCE / f"{prefix}_{n}.inp").read_text(encoding="utf-8")
        stack = []
        for raw in text.splitlines():
            line = raw.split("!")[0].strip()
            if line.startswith("&END"):
                assert line.split()[1] == stack.pop(), line
            elif line.startswith("&"):
                stack.append(line.split()[0][1:])
        assert not stack
        assert text.count("PERIODIC NONE") == 2
        assert "TYPE DFTD3(BJ)" in text
        assert "ABC [angstrom] 40 35 47" in text
        if prefix == "step":
            ids = [int(i) for line in text.splitlines() if line.strip().startswith("LIST ") for i in line.split()[1:]]
            assert ids == cfg["bottom_ids"] + cfg["grip_ids"]
            assert f"COORD_FILE_NAME step_{n}.xyz" in text
            assert "RUN_TYPE GEO_OPT" in text
        else:
            assert f"COORD_FILE_NAME relaxed_{n:02d}.xyz" in text
            assert "RUN_TYPE ENERGY_FORCE" in text
for element in ("C", "Si"):
    assert re.search(rf"^\s*{element}\s+[^\n]*\bDZVP-MOLOPT-SR-GTH\b", (SOURCE / "BASIS_MOLOPT").read_text(), re.M)
    assert re.search(rf"^\s*{element}\s+[^\n]*\bGTH-PBE-q4\b", (SOURCE / "GTH_POTENTIALS").read_text(), re.M)
for folder in (SOURCE, HERE.parent / "边缘面拉伸 - 提交"):
    manifest = json.loads((folder / "data_sources.json").read_text(encoding="utf-8"))
    for name, info in manifest["files"].items():
        assert hashlib.sha256((folder / name).read_bytes()).hexdigest() == info["sha256"]
print("PASS: 42 input structures, atom groups, coordinate translation, CP2K data hashes/basis names")

fake_program = '''from pathlib import Path
import re, sys
args = sys.argv
inp = Path(args[args.index('-i')+1])
out = Path(args[args.index('-o')+1])
t = inp.read_text()
geo = 'RUN_TYPE GEO_OPT' in t
name = re.search(r'^  PROJECT (.+)$', t, re.M).group(1)
source = Path(re.search(r'COORD_FILE_NAME (\\S+)', t).group(1))
if Path('fail_here').exists() and inp.name == 'step_1.inp':
    out.write_text('SCF run NOT converged\\nPROGRAM ENDED AT mock\\n')
    sys.exit(0)
if geo:
    rows = source.read_text().splitlines()
    p = rows[2].split()
    p[1] = f'{float(p[1]) + .001:.10f}'
    rows[2] = ' '.join(p)
    target = re.search(r'FILENAME =(\\S+)', t).group(1)
    Path(target).write_text('\\n'.join(rows)+'\\n')
Path(name+'-RESTART.wfn').write_text('MOCK WFN; NOT PHYSICAL DATA')
out.write_text('SCF run converged in 5 steps\\n' + ('GEOMETRY OPTIMIZATION COMPLETED\\n' if geo else '') + 'ENERGY| Total FORCE_EVAL ( QS ) energy [hartree]: -100.000000\\nPROGRAM ENDED AT mock\\n')
'''

temp = Path(tempfile.mkdtemp(prefix="workflow_test_", dir=HERE)).resolve()
try:
    def fixture(name):
        d = temp / name
        d.mkdir()
        for p in SOURCE.iterdir():
            if p.is_file() and (p.suffix == ".inp" or p.name in ("initial_vdm.xyz", "step_0.xyz", "run_quasistatic.py", "loading_protocol.json", "data_sources.json", "BASIS_MOLOPT", "GTH_POTENTIALS", "dftd3.dat")):
                shutil.copy2(p, d / p.name)
        (d / "mock_cp2k.py").write_text(fake_program)
        return d

    def run(d):
        command = f'"{Path(sys.executable).as_posix()}" "{(d / "mock_cp2k.py").as_posix()}"'
        return subprocess.run([sys.executable, str(d / "run_quasistatic.py"), "--run", "--cp2k-command", command], capture_output=True, text=True)

    good = fixture("success")
    result = run(good)
    assert result.returncode == 0, result.stdout + result.stderr
    final = runner.read_xyz(good / "relaxed_20.xyz")
    runner.validate_atoms(final, cfg, ref, 20)
    assert abs(final[0][1] - ref[0][1] - .021) < 1e-8, "Free atom updates were not inherited"
    assert len(list(good.glob("step_*.done.json"))) == 21
    assert len(list(good.glob("sp_*.done.json"))) == 21
    assert len((good / "energy_vs_displacement.csv").read_text().splitlines()) == 22
    resumed = run(good)
    assert resumed.returncode == 0 and resumed.stdout.count("Verified completed stage") == 42, resumed.stderr
    (good / "step_1.inp").write_text((good / "step_1.inp").read_text() + "\n! changed\n")
    assert run(good).returncode != 0, "Changed inputs must not reuse completed results"
    bad = fixture("failure")
    (bad / "fail_here").touch()
    failed = run(bad)
    assert failed.returncode != 0 and not (bad / "step_2.xyz").exists()
    partial = temp / "truncated.xyz"
    partial.write_text("2\nincomplete\nC 0 0 0\n")
    try:
        runner.read_xyz(partial)
    except ValueError:
        pass
    else:
        raise AssertionError("Truncated XYZ accepted")
    print("PASS: mock 20-step load/relax/SP, 4 A grip displacement, inherited free atoms, resume, changed-input refusal, SCF failure stop, truncated XYZ rejection")
finally:
    # Verify the absolute target is a dedicated test directory inside this workspace.
    if temp.parent == HERE.resolve() and temp.name.startswith("workflow_test_"):
        shutil.rmtree(temp)
    else:
        raise RuntimeError(f"Unsafe temporary cleanup path: {temp}")
print("NO CP2K/DFT calculations were performed.")
