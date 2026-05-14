#!/usr/bin/env python3
"""Verify if fixed potential values appear in model_7 steps 1-8"""

from pathlib import Path

outdir = Path('output/model_7')

# Check for fixed values in different steps
for step in [1, 2, 3, 8, 9, 10, 200]:
    fname = outdir / f'tempa.dat{step:03d}'
    if not fname.exists():
        print(f'{fname.name}: NOT FOUND')
        continue
    
    text = fname.read_text()
    tokens = text.split()
    
    # Count occurrences of fixed values
    count_pos_1 = sum(1 for t in tokens if t.strip().upper().replace('D','E') == '1.0000E+00')
    count_neg_1 = sum(1 for t in tokens if t.strip().upper().replace('D','E') == '-1.0000E+00')
    
    print(f'{fname.name}: +1.0E:{count_pos_1:3d}, -1.0E:{count_neg_1:3d}  (total tokens: {len(tokens)})')

print('\n--- Detailed check for step 1 ---')
fname = outdir / 'tempa.dat001'
text = fname.read_text()
tokens = text.split()

# Find indices where fixed values appear
pos_indices = []
neg_indices = []
for i, t in enumerate(tokens):
    tu = t.strip().upper().replace('D','E')
    if tu == '1.0000E+00':
        pos_indices.append(i)
    elif tu == '-1.0000E+00':
        neg_indices.append(i)

print(f'Indices of +1.0000E+00: {pos_indices[:10]}')
print(f'Indices of -1.0000E+00: {neg_indices[:10]}')

# Check if nodes match expected material assignments
print('\n--- Material block info from sin.dat ---')
print('Material 1 (fixed at -1.0V): elements ~291109-316818')
print('Material 2 (fixed at +1.0V): elements ~291133-316842')
print('Each 4-node tetrahedron contributes ~4 nodes to fixed_node set')
