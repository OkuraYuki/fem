from pathlib import Path
import numpy as np, json, hashlib
root=Path(__file__).resolve().parents[1]
model=root/'models/model_entyu_1'; out=root/'output/poweroff_20261026'
with (model/'in.dat').open() as f:
    header=f.readline().split(); n,e=int(header[0]),int(header[1]); unit=float(header[4])
    data=np.fromstring(f.read(),sep=' ')
conn=data[:4*e].astype(np.int64).reshape(-1,4)-1
xyz=data[4*e:4*e+3*n].reshape(-1,3)
lines=[x for x in (model/'sin.dat').read_text().splitlines() if x.strip()]
b,nb,nm=map(int,lines[0].split()[:3]);bc=np.array([list(map(float,x.split())) for x in lines[1:1+b]])
ranges=np.array([list(map(int,x.split())) for x in lines[1+b:1+b+nb]])
mat=np.array([list(map(float,x.split())) for x in lines[1+b+nb:1+b+nb+nm]])
assignment=np.zeros(e,dtype=int);coverage=np.zeros(e,dtype=int)
for start,end,m in ranges:assignment[start-1:end]=m;coverage[start-1:end]+=1
members={m:np.unique(conn[assignment==m]) for m in range(1,nm+1)}
vol=[]
for start in range(0,e,50000):
    points=xyz[conn[start:start+50000]]; jac=np.transpose(points[:,1:]-points[:,:1],(0,2,1));vol.append(np.linalg.det(jac)/6)
vol=np.concatenate(vol)
probe=[181500,181254,180954,180733,180684]
info={
 'nodes':n,'elements':e,'header_unit_unapplied':unit,'coordinates_min':xyz.min(axis=0).tolist(),'coordinates_max':xyz.max(axis=0).tolist(),
 'connectivity_in_range':bool(np.all((conn>=0)&(conn<n))), 'coverage_min':int(coverage.min()),'coverage_max':int(coverage.max()),
 'material_element_counts':{str(m):int(np.count_nonzero(assignment==m)) for m in members},
 'positive_finite_materials':bool(np.all(np.isfinite(mat)) and np.all(mat>0)),
 'conductivity_permittivity':mat.tolist(), 'min_abs_volume':float(np.abs(vol).min()),'max_abs_volume':float(np.abs(vol).max()),
 'degenerate_elements':int(np.count_nonzero(np.abs(vol)<1e-20)), 'negative_orientation_elements':int(np.count_nonzero(vol<0)),
 'boundary_type_counts':{str(int(k)):int(np.count_nonzero(bc[:,2]==k)) for k in np.unique(bc[:,2])},
 'pulse_ground_shared_nodes':int(len(np.intersect1d(members[1],members[3]))),
 'pulse_explicit_ground_conflicts':int(len(np.intersect1d(members[1],bc[:,:2].astype(int).ravel()-1))),
 'probe_nodes':[{ 'node':v,'xyz':xyz[v-1].tolist(),'materials':[m for m in members if v-1 in members[m]]} for v in probe],
 'sha256':{f:hashlib.sha256((model/f).read_bytes()).hexdigest() for f in ['in.dat','sin.dat','tmate.dat','pulse.dat','sina.dat.t']}
}
(out/'mesh_audit.json').write_text(json.dumps(info,indent=2));print(json.dumps(info,indent=2))
