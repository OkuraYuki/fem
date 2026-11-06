"""Reproduce the nine cases in a new, non-existing output directory.
Usage: python3 tools/run_poweroff_suite.py /absolute/new/output/path
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess,sys,json
root=Path(__file__).resolve().parents[1]
if len(sys.argv)!=2:raise SystemExit(__doc__)
dest=Path(sys.argv[1]).resolve();dest.mkdir(parents=True,exist_ok=False)
spec=json.loads((root/'output/poweroff_20261026/cases.json').read_text())
(dest/'cases.json').write_text(json.dumps(spec,indent=2))
runner=dest/'run_case'
subprocess.check_call(['g++','-O2','-std=c++17',str(root/'tools/poweroff_case.cpp'),str(root/'src/electrostatic_analyzer.cpp'),'-o',str(runner)])
def run(c):
 d=dest/c[0];d.mkdir()
 cmd=[str(runner),str(root/'models/model_entyu_1'),str(d)]+[str(v) for v in c[1:]]
 with (d/'run.log').open('w') as log:
  code=subprocess.call(cmd,stdout=log,stderr=subprocess.STDOUT)
 print(c[0],code,flush=True)
 return code
with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(run,spec['cases']))
if any(results):raise SystemExit(1)
