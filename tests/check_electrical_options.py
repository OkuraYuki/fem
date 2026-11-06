"""CLI checks; argv[1] is the candidate executable, argv[2] is the results directory."""
from pathlib import Path
import subprocess,tempfile,csv,json,sys,math
exe=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]);results=[]
config='''electrode_material 1
series_resistance_ohm ideal
disconnect_time_s 3e-9
off_resistance_ohm 900
external_capacitance_F 0
coordinate_scale_to_m 1
initial_state steady
'''
for name in ['valid','no_resistor','misaligned_switch','missing_off_resistance','invalid_scale','unknown_option','invalid_symmetry','material_conflict']:
 d=Path(tempfile.mkdtemp(prefix='fem-electrical-'+name+'-'))
 (d/'in.dat').write_text('9 3 0 0 1\n1 2 3 4\n2 5 6 7\n3 4 8 9\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n2 0 0\n1 1 0\n1 0 1\n-1 1 0\n0 1 1\n')
 (d/'sin.dat').write_text('0 3 3 0 0\n1 1 2\n2 2 1\n3 3 3\n5.8e7 5.8e7 5.8e7 8.9e-12 8.9e-12 8.9e-12\n.01 .01 .01 1e-10 1e-10 1e-10\n5.8e7 5.8e7 5.8e7 8.9e-12 8.9e-12 8.9e-12\n0\n')
 (d/'tmate.dat').write_text('1\n3 0\n');(d/'pulse.dat').write_text('1 0 1 0 3e-9 0 100e-9\n');(d/'sina.dat.t').write_text('0 0\n10 0\n1e-9 0 0\n0 0\n')
 c=config
 if name=='no_resistor':c=c.replace('off_resistance_ohm 900','off_resistance_ohm open')
 if name=='misaligned_switch':c=c.replace('disconnect_time_s 3e-9','disconnect_time_s 3.5e-9')
 if name=='missing_off_resistance':c=c.replace('off_resistance_ohm 900\n','')
 if name=='invalid_scale':c=c.replace('coordinate_scale_to_m 1','coordinate_scale_to_m 0')
 if name=='unknown_option':c+='invented_option 3\n'
 if name=='invalid_symmetry':(d/'symmetry.dat').write_text('axis z\n')
 if name=='material_conflict':(d/'tmate.dat').write_text('1\n2 0\n')
 (d/'circuit.dat').write_text(c)
 p=subprocess.run([str(exe)],cwd=str(d),stdout=subprocess.PIPE,stderr=subprocess.PIPE,universal_newlines=True)
 (d/'run.log').write_text(p.stdout+p.stderr)
 valid=name in ['valid','no_resistor'];passed=p.returncode==(0 if valid else 1)
 if valid and passed:
  with (d/'circuit_history.csv').open() as f:rows=list(csv.DictReader(f))
  G=.01/9;C=1e-10/9;expected=1
  passed=len(rows)==11
  for i,row in enumerate(rows):
   if i>3:expected/=1+1e-9*(G+(0 if name=='no_resistor' else 1./900))/C
   passed=passed and abs(float(row['electrode_V'])-expected)<2e-9 and abs(float(row['electrode_charge_C'])-C*expected)<2e-18
 else:
  passed=passed and not (d/'tempa.dat001').exists()
 results.append({'case':name,'passed':bool(passed),'exit_code':p.returncode,'directory':str(d)})
(out/'electrical_cli_checks.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
if not all(r['passed'] for r in results):raise SystemExit(1)
