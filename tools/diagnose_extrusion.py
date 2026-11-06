from pathlib import Path
import numpy as np,json
root=Path(__file__).resolve().parents[1];p=root/'models/model_entyu_1';out=root/'output/poweroff_20261029'
with (p/'in.dat').open() as f:
 h=f.readline().split();n,e=int(h[0]),int(h[1]);a=np.fromstring(f.read(),sep=' ')
c=a[:4*e].astype(np.int64).reshape(-1,4)-1;x=a[4*e:4*e+3*n].reshape(-1,3)
s=[z for z in (p/'sin.dat').read_text().splitlines() if z.strip()];nb,nr=map(int,s[0].split()[:2]);m=np.zeros(e,int)
for z in s[1+nb:1+nb+nr]:
 first,last,mat=map(int,z.split());m[first-1:last]=mat
q=np.sort(c%2409,axis=1);distinct=1+(q[:,1:]!=q[:,:-1]).sum(axis=1)
info={'projection_distinct_xy_vertices':{str(k):int(np.count_nonzero(distinct==k)) for k in np.unique(distinct)},'interfaces':{}}
for mat in [1,3]:
 faces={}
 for tet in c[m==mat]:
  for skip in range(4):
   face=tuple(sorted(tet[np.arange(4)!=skip]));faces[face]=faces.get(face,0)+1
 boundary=np.array([k for k,v in faces.items() if v==1]);points=x[boundary];normal=np.cross(points[:,1]-points[:,0],points[:,2]-points[:,0]);side=np.ptp(points[:,:,2],axis=1)>0
 info['interfaces'][str(mat)]={'boundary_faces':len(boundary),'side_faces':int(side.sum()),'side_faces_with_nonzero_z_normal':int(np.count_nonzero(np.abs(normal[side,2])>1e-7)),'max_abs_side_z_normal':float(np.max(np.abs(normal[side,2])))}
phi=np.frombuffer((root/'output/poweroff_20261026/base/tempa.dat001').read_bytes().replace(b'\n',b''),dtype='S12').astype(float).reshape(151,2409)
span=np.ptp(phi,axis=0);i=int(np.argmax(span));info['step1_z_variation']={'max_range_V':float(span[i]),'cross_section_node_1based':i+1,'min_V':float(phi[:,i].min()),'max_V':float(phi[:,i].max()),'range_at_peak_node884_V':float(span[883])}
(out/'extrusion_geometry.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))
