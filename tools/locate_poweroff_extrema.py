from pathlib import Path
import csv
import numpy as np
root=Path(__file__).resolve().parents[1];out=root/'output/poweroff_20261026'
with (root/'models/model_entyu_1/in.dat').open() as f:
 h=f.readline().split();n,e=int(h[0]),int(h[1]);raw=np.fromstring(f.read(),sep=' ')
xyz=raw[4*e:4*e+3*n].reshape(-1,3)
rows=[]
for case,folder in [('old',root/'output/model_entyu_1'),('fixed',out/'base')]:
 for step in [1,24,25,50,75,100,150]:
  p=folder/f'tempa.dat{step:03}'
  if not p.exists():continue
  data=np.frombuffer(p.read_bytes().replace(b'\n',b'').replace(b'\r',b''),dtype='S12').astype(float)
  if len(data)!=n:continue
  for kind,idx in [('minimum',int(np.argmin(data))),('maximum',int(np.argmax(data)))]:
   rows.append({'case':case,'time_ns':step,'extremum':kind,'node':idx+1,'potential_V':float(data[idx]),'x_mesh':float(xyz[idx,0]),'y_mesh':float(xyz[idx,1]),'z_mesh':float(xyz[idx,2])})
with (out/'extrema_locations.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
for r in rows:
 if r['time_ns']==25:print(r)
