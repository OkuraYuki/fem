from pathlib import Path
import subprocess,tempfile,json
root=Path(__file__).resolve().parents[1]
template={
'in.dat':'4 1 0 0 1\n1 2 3 4\n0 0 0\n1 0 0\n0 1 0\n0 0 1\n',
'sin.dat':'2 1 1 0 0\n2 2 1 1\n3 4 4 0\n1 1 1\n.01 .01 .01 1e-10 1e-10 1e-10\n0\n',
 'tmate.dat':'0\n', 'sina.dat.t':'0 0 0\n5 0 0\n1e-9 0 0\n0 0\n',
 'voltage_difference.dat':'3 1 free_potential\n'}
results=[]
for name in ['valid','time0','bad_config','nonzero_flux','unsupported_bc']:
 d=Path(tempfile.mkdtemp(prefix='fem-cli-'+name+'-'))
 for f,contents in template.items():(d/f).write_text(contents)
 if name=='time0':(d/'sina.dat.t').write_text('# comment\n0 0 0\n5 0 0 # nstep\n1e-9 2e-9 0 # dt,time0\n0 0\n')
 if name=='bad_config':(d/'sina.dat.t').write_text('0 0\n')
 if name=='nonzero_flux':(d/'sin.dat').write_text((d/'sin.dat').read_text().replace('2 2 1 1','2 2 5 1'))
 if name=='unsupported_bc':(d/'sin.dat').write_text((d/'sin.dat').read_text().replace('2 2 1 1','2 2 99 1'))
 p=subprocess.run([str(root/'bin/fem')],cwd=d,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 (d/'run.log').write_text(p.stdout+p.stderr)
 expected=0 if name in ['valid','time0'] else 1
 ok=(p.returncode==expected)
 if name in ['valid','time0']:
  vals=[float(x) for x in (d/'tempa.dat001').read_text().split()]
  ok=ok and abs(vals[0]-1/3)<5e-6
  if name=='time0':ok=ok and 'final_time: 3.000000e-09' in (d/'analysis_summary.txt').read_text()
 else:ok=ok and not (d/'tempa.dat001').exists()
 results.append({'case':name,'passed':ok,'returncode':p.returncode,'directory':str(d)})
(root/'output/poweroff_20261026/cli_checks.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
if not all(x['passed'] for x in results):raise SystemExit(1)
