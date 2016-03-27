#!/usr/bin/python
"""
svd.py

    A multiconformer SVD fitting program.

    Author: Justin Lorieau, (c) 2012
"""

import sys
from collections import OrderedDict
from scipy import linalg, array, dot, inf
from math import sqrt
from mollib import *

truncate = lambda s,length,suffix='...': s if len(s) < length else \
    s[:length-len(suffix)].rsplit(' ', 1)[0]+ suffix

scaling_factors= {('N', 'HN'): 1.0000,
                  ('CA', 'HA'): 0.48463594}


params = {}
current_arg = None
for arg in sys.argv[1:]:
    if arg[0] == '-':
        current_arg = arg
        continue
    if current_arg is None:
        continue

    param_list = params.setdefault(current_arg, [])
    param_list.append(arg)

for opt in ('-pdb', '-dc'):
    if opt not in params:
        print "\n".join(("svd.py",
                         "\t-pdb\tA list of PDB files.",
                         "\t-dc\tA dipolar rdc table."))
        exit()

# Read in the rdcs
class RDC(object):

    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
        self.rdc = float(self.rdc)
        self.resid_i = int(self.resid_i)
        self.resid_j = int(self.resid_j)
        
        # Determine the normalization scaling_factor for this dipolar coupling
        name_i, name_j = self.name_i, self.name_j
        if (name_i.strip('#'),name_j.strip('#')) in scaling_factors:
            self.scaling = scaling_factors[name_i.strip('#'),name_j.strip('#')]
        elif (name_j.strip('#'),name_i.strip('#')) in scaling_factors:
            self.scaling = scaling_factors[name_j.strip('#'),name_i.strip('#')]
        else:
            self.scaling = 1.0
        
        # A list of conformers and their dipoles that contribute to this
        # rdc
        self.conformers = OrderedDict()

coupling_key = lambda rdc: (
    "%s%s-%s%s" % (rdc.resid_i, rdc.atom_i,
                   rdc.resid_j, rdc.atom_j))
rdcs = []

with open(params['-dc'][0]) as dc_file:
    for line in dc_file.readlines():
        match = re.search(r'^\s*(?P<resid_i>[0-9]+)'
                           '\s*[A-Z]{3}'
                           '\s*(?P<name_i>[A-Z\#]{1,4})'
                           '\s*(?P<resid_j>[0-9]+)'
                           '\s*[A-Z]{3}'
                           '\s*(?P<name_j>[A-Z\#]{1,4})'
                           '\s*(?P<rdc>-?[0-9\.]+)', line)
        if match is None:
            continue
        rdc = RDC(**match.groupdict())
        rdcs.append(rdc)
        
# Read in the dipoles from the PDB files
class Dipole(object):

    
    def __init__(self, structure, atom_i, atom_j, scaling=1.0):
        self.atom_i = atom_i
        self.atom_j = atom_j
        self.structure = structure
        self.scaling = scaling

    @property
    def x(self):
        return self.atom_j.x - self.atom_i.x
    
    @property
    def y(self):
        return self.atom_j.y - self.atom_i.y
    
    @property
    def z(self):
        return self.atom_j.z - self.atom_i.z

    @property
    def r(self):
        return sqrt(self.x*self.x + self.y*self.y + self.z*self.z)

    @property
    def cos_x(self):
        return self.x/self.r

    @property
    def cos_y(self):
        return self.y/self.r

    @property
    def cos_z(self):
        return self.z/self.r

#dipole_key = lambda dipole
structures = []

for pdb_file in params['-pdb']:
    structure = Molecule(pdb_file)
    structures.append(structure)

def get_dipoles(structure, resid_i, name_i, resid_j, name_j):
    """Returns a list of dipole objects that are in the given a structure, and 
    selected by the atom name_i and name_j name strings (i.e. 'HA') and residue
    resid (i.e. 6). 
    
    Works with '#' selections--'HAX' would select 'HA1' and HA2'.
    """
    # Use the first change only
    chain_name = structure.keys()[0]
    mol = structure[chain_name]
    
    if name_i.endswith('#'):
        name_i = name_i.strip('#')
        atoms_i = [a for a in mol[resid_i].atoms if a.name.startswith(name_i)]
    else:
        atoms_i = [a for a in mol[resid_i].atoms if a.name == name_i]
    
    if name_j.endswith('#'):
        name_j = name_j.strip('#')
        atoms_j = [a for a in mol[resid_i].atoms if a.name.startswith(name_j)]
    else:
        atoms_j = [a for a in mol[resid_i].atoms if a.name == name_j]
    
    return [Dipole(structure=structure, atom_i=atom_i, atom_j=atom_j)
            for atom_i in atoms_i for atom_j in atoms_j]
    
# Go through the dipolar rdcs and find the corresponding dipoles in the
# structures
for rdc in rdcs:
    for structure in structures:
        # Use the first chain only
        chain = structure.keys()[0]
        resid_i, name_i = rdc.resid_i, rdc.name_i
        resid_j, name_j = rdc.resid_j, rdc.name_j

        rdc.conformers[structure.name] = get_dipoles(structure, 
                                                     resid_i, name_i, 
                                                     resid_j, name_j)

# Construct the dipole vector and structure matrix A.
D = []
A = []

for rdc in rdcs:
    D.append(rdc.rdc*rdc.scaling)
    A_line = []
    for name, dipoles in rdc.conformers.items():
        
        A_line.extend([sum([0.5*(3.*d.cos_z**2 - 1.) for d in dipoles]),
                       sum([0.5*(d.cos_x**2 - d.cos_y**2) for d in dipoles]),
                       sum([2.*d.cos_x*d.cos_y for d in dipoles]),
                       sum([2.*d.cos_x*d.cos_z for d in dipoles]),
                       sum([2.*d.cos_y*d.cos_z for d in dipoles])])
    A.append(A_line)

# Perform the linear algebra on this vector and matrix
A=array(A)
D=array(D)
U, w, V = linalg.svd(A, full_matrices=False)

# Remove inf items from the inverted 'w' vector
w_inv = 1./w
for i in range(w_inv.size):
    w_inv[i] = 0. if w_inv[i] == inf else w_inv[i]
w_inv = linalg.diagsvd(w_inv, w.shape[0], w.shape[0])

# Calculate the Saupe matrix
A_inv = dot(V.transpose(), dot(w_inv, U.transpose()))
Saupe = dot(A_inv, D) # Saupe matrix

# Reconstruct the 3x3 Saupe matrix for each conformer and diagonalize
S_xyz = []
Da, Dr, Rh = [], [], []
for x in xrange(0, len(Saupe), 5):
    s = Saupe[x:x+5]
    s_xyz = array([[ -0.5*(s[0]-s[1]), s[2],              s[3],],
                   [s[2]             , -0.5*(s[0]+s[1] ), s[4]],
                   [s[3],              s[4],              s[0]]])
    s_xyz = linalg.eigvals(s_xyz).real
    S_xyz.append(s_xyz)
    
    xx,yy,zz = [i for i in sorted(abs(s_xyz))]
    da = max(s_xyz)/2. if max(s_xyz) == zz else min(s_xyz)/2.
    dr = (yy-xx)/3.
    Da.append(da)
    Dr.append(dr)
    Rh.append(dr/abs(da))

# Calculate the predicted RDCs, and scale them to their original values with
# the scaling factor.
D_pred = dot(A, Saupe)
delta_D = []
for i, rdc in enumerate(rdcs):
    delta_D.append( D_pred[i] - rdc.rdc*rdc.scaling)
    D_pred[i] = D_pred[i]/rdc.scaling
    
rms = sqrt(sum([d**2 for d in delta_D])/(len(delta_D)-1.))

# Calculate the multi-conformer Q-factor
N_RDC = len(delta_D)
Q = sqrt(sum([delta**2 for delta in delta_D])/ \
    (N_RDC * sum([da**2 *(4.+3.*rh**2)/5. for da,rh in zip(Da, Rh)])))

### Print out the results ###
print '\n'.join(['Saupe({}): {}'.format(count ,Saupe[x:x+5])
                 for count, x in enumerate(xrange(0, len(Saupe), 5), 1)])

print '\n'.join(['S_xyz({}): {}'.format(count, s_xyz)
                 for count, s_xyz in enumerate(S_xyz, 1)])
print '\n'.join(['Da({c}): {Da:.3f}Hz, Dr({c}): {Dr:.3f}, Rh({c}): {Rh:.3f}'.format(
    c=count,Da=da,Dr=dr,Rh=rh) 
    for count, (da,dr,rh) in enumerate(zip(Da, Dr, Rh),1)]) 
print 'Q-factor: {:.1f}%'.format(Q*100.)
print 'RMS (Hz):', rms

header = "{:<20} {:>10} / {:<10} {:>12} {:>12} {:>12} {:>12}".format(
    'conformer:', 'atom(i)', 'atom(j)', 'D_obs (Hz)', 'D_pred (Hz)', 
    'Delta (Hz)', 'Scale')
header = '\n'.join(('\n', '='*len(header), header, '='*len(header)))
print header
for rdc, pred in zip(rdcs, D_pred):
    dipole_labels = []
    for name, dipoles in rdc.conformers.items():
        for dipole in dipoles:
            dipole_labels.append("{name:<20} {i:>10} / {j:<10}".format(
                name=truncate(name,20),
                i=dipole.atom_i, j=dipole.atom_j))

    
    print "{} {:>12.3f} {:>12.3f} {:>12.3f} {:>12.3f}".format(dipole_labels[0], 
        rdc.rdc, pred, rdc.rdc-pred, rdc.scaling)
    for dipole_label in dipole_labels[1:]:
        print dipole_label
        



