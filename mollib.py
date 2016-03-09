"""
MolLib for Python
=================

Author: J Lorieau

Copyright 2011

MolLib are a simply Python toolset for manipulating molecular
structures in Python. A Molecule is a dict of Chains, which is a dict of
Residues, which is a dict of Atoms. The Molecule is constructed with helper
read and write functions in the Molecule class, and more sophisticated behavior
can be added to the Molecule, Chain, Residue and Atom classes by deriving them
and changing the chain_class, residue_class and atom_class of Molecule.

>>> mol=Molecule(pdb_code='2OED')
>>> print mol
Molecule:    1 chains, 56 residues, 862 atoms.
>>> print mol['A'][1], mol['A'][1].atom_size,
M1 19
>>> print mol['A'][2]['CA'], mol['A'][2]['CA'].mass
Q2-CA 12.01
>>> print mol.mass, 'Da'
6206.75 Da
>>> print "({:.3f}, {:.3f}, {:.3f})".format(*mol.center_of_mass)
(0.133, -0.347, -0.002)
>>> mol.rotate_zyz(0,90,0)
"""

import re
from itertools import imap, ifilter
from itertools import chain as ichain
from math import cos, sin, sqrt, pi
import numpy as np

import urllib
import os.path

### Utility Functions

def convert(s):
    """Converts a string 's' into either an integer, float or string"""
    if isinstance(s, str):
        s=s.strip()
    else:
        return None
    for t in int, float, str:
        try:
            return t(s)
        except ValueError:
            continue
    return None


### MolLib Implementation

class PrimitiveMetaClass(type):
    """A Primitive Metaclass for properly assigning and collecting attributes,
    such as 'optional'"""
    def __new__(meta, classname, bases, classDict):

        # This code adds 'optional' tuples from parent classes to the
        # 'optional' attribute of this class.
        parent_optional = tuple(*[getattr(base, 'optional')
                                 for base in bases
                                 if hasattr(base, 'optional')])
        classDict['optional'] = classDict['optional'] + parent_optional \
            if 'optional' in classDict else parent_optional

        return type.__new__(meta, classname, bases, classDict)

class Primitive(object):
    "The base for all objects."

    __metaclass__ = PrimitiveMetaClass

    __slots__ = ()
    optional = ()

    def __init__(self, **kwargs):

        # Check that the required arguments have been specified
        req_kwargs = [kw for kw in self.__slots__ if kw not in self.optional]
        assert all(kw in kwargs for kw in req_kwargs), \
            "All of the following parameters are needed: {}"\
            .format(req_kwargs)

        # Assign the values
        [setattr(self, kw, value) for kw,value in kwargs.items()
         if kw in self.__slots__]

        super(Primitive, self).__init__()


class Atom(Primitive):
    "An atom in a residue."

    __slots__ = ('number', 'name', 'x', 'y', 'z', 'charge', 'element',
                 'residue', 'chain', 'molecule')
    optional  = ('charge', 'residue', 'chain', 'molecule')

    atom_Mw = {'H': 1.01, 'C': 12.01, 'N': 14.01, 'O': 16.00, 'Na': 22.99,
               'P': 30.97, 'S': 32.07, 'Cl': 35.45}

    def __repr__(self):
        return u"{}-{}".format(self.residue,self.name) if self.residue else \
            u"{}".format(self.name)

    @property
    def mass(self):
        return self.atom_Mw[self.element]


class Residue(dict):
    "A residue in a chain."

    one_letter_codes = {'ALA': 'A', 'GLY': 'G', 'SER': 'S', 'THR': 'T',
                        'MET': 'M', 'CYS': 'C', 'ILE': 'I', 'LEU': 'L',
                        'VAL': 'V', 'PHE': 'F', 'TYR': 'Y', 'TRP': 'W',
                        'ASN': 'N', 'GLN': 'Q', 'ASP': 'D', 'GLU': 'E',
                        'HIS': 'H', 'PRO': 'P', 'ARG': 'R', 'LYS': 'K'}

    def __init__(self, name, number, *args, **kwargs):
        # The amino-acid must be one of the 20 known amino-acids.
        name = name.upper()
        assert name in self.one_letter_codes, \
            "Amino-acid {} not recognized".format(name)

        self.name = name                                # full name, i.e. MET
        self.letter =  self.one_letter_codes[self.name] # letter code, ex: M
        self.number = number
        super(Residue, self).__init__(*args, **kwargs)

    def __repr__(self):
        return u"{}{}".format(self.letter, self.number)

    @property
    def atoms(self):
        """Returns an iterator over all atoms in this residue,
        sorted by atom number"""
        return (atom for atom in sorted(self.values(), key=lambda a:a.number))

    @property
    def atom_size(self):
        return len(list(self.atoms))


class Chain(dict):
    "A chain in a molecule."

    def __init__(self, id, *args, **kwargs):
        self.id = id
        super(Chain, self).__init__(*args, **kwargs)

    def __repr__(self):
        return u"{}".format(self.id)

    @property
    def residues(self):
        """Returns an iterator over all residues in this chain,
        sorted by residue number"""
        return (residue for residue in
                sorted(self.values(), key=lambda a: a.number))

    @property
    def residue_size(self):
        return len(list(self.residues))

    @property
    def atoms(self):
        """Returns an iterator over all atoms in this chain,
        sorted by atom number"""
        return (a for a in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda a: a.number))

    @property
    def atom_size(self):
        return len(list(self.atoms))

class Molecule(dict):
    "A class for molecular structures."

    # The following class-level attributes are used to customize the base or
    # derived Chain, Residue and Atom classes used in the molecule
    chain_class = Chain
    residue_class = Residue
    atom_class = Atom

    def __init__(self, filename=None, pdb_code=None, *args, **kwargs):
        if filename is not None:        
            self.read_pdb(filename)
        elif pdb_code is not None:
            self.fetch_pdb(pdb_code)
        else:
            raise
        super(Molecule, self).__init__(*args, **kwargs)

    def __repr__(self):
        return (u"Molecule:"
                 "    {} chains, {} residues, {} atoms." \
                 .format(self.chain_size, self.residue_size, self.atom_size))


    ### Basic Accessors and Mutators ###

    @property
    def chain_size(self):
        """Returns the number of chains.

        >>> mol=Molecule(pdb_code='1HTM') # Influenza hemagglutinin, 6 subunits
        >>> print mol.chain_size
        6
        """
        return len(self)

    @property
    def residues(self):
        """Returns an iterator over all residues in this molecule,
        sorted by residue number.
        
        >>> mol=Molecule(pdb_code='2KXA')
        >>> print [r.number for r in mol.residues]
        [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]
        """
        return (r for r in sorted(ichain(*[r.values() for r in self.values()]),
                                  key=lambda r: r.number))

    @property
    def residue_size(self):
        return len(list(self.residues))

    @property
    def atoms(self):
        """Returns an iterator over all atoms in this molecule,
        sorted by atom number"""
        return (a for a in sorted(ichain(*[r.values() for r in
                ichain(*[c.values() for c in self.values()])]),
                key=lambda a:a.number))

    @property
    def atom_size(self):
        return len(list(self.atoms))


    ### Molecular Properties ###
    #TODO: replace atoms with different isotopes

    @property
    def mass(self):
        """ Returns the mass of the molecule.

        >>> mol = Molecule(pdb_code='2KXA')
        >>> print mol.mass
        2445.07
        """
        return sum(a.mass for a in self.atoms)

    @property
    def center_of_mass(self):
        """ Returns the center-of-mass x,y,z vector of the molecule.

        >>> mol = Molecule(pdb_code='2KXA')
        >>> print "{:.3f} {:.3f} {:.3f}".format(*mol.center_of_mass)
        16.970 0.070 0.122
        """        
        
        x,y,z=(0,0,0)
        m_total = 0.
        for atom in self.atoms:
            mass = atom.mass
            m_total += mass
            x += atom.x*mass
            y += atom.y*mass
            z += atom.z*mass
        return (x/m_total, y/m_total, z/m_total)


    ### Mutator Functions ###

    def center(self):
        """Centers a molecule about its center_of_mass.
        
        >>> mol = Molecule(pdb_code='2KXA')
        >>> print "{:.3f} {:.3f} {:.3f}".format(*mol.center_of_mass)
        16.970 0.070 0.122
        >>> mol.center()
        >>> print "{:.3f} {:.3f} {:.3f}".format(*mol.center_of_mass)
        0.000 0.000 0.000
        """
        com = self.center_of_mass
        for atom in self.atoms:
            atom.x -= com[0]
            atom.y -= com[1]
            atom.z -= com[2]

    def rotate_zyz(self, alpha, beta, gamma):
        "Rotates a molecule by the Euler z-y-z angles in degrees."

        # Trig.
        sin_a = sin(alpha*pi/180.)
        sin_b = sin(beta*pi/180.)
        sin_g = sin(gamma*pi/180.)

        cos_a = cos(alpha*pi/180.)
        cos_b = cos(beta*pi/180.)
        cos_g = cos(gamma*pi/180.)

        m = np.matrix([[-sin_a * sin_g + cos_a * cos_b * cos_g,
                        cos_a * sin_g + sin_a * cos_b * cos_g,
                        -sin_b * cos_g],
                       [-sin_a * cos_g - cos_a * cos_b * sin_g,
                        cos_a * cos_g - sin_a * cos_b * sin_g,
                        sin_b * sin_g],
                       [cos_a * sin_b,
                        sin_a * sin_b,
                        cos_b]])
        for atom in self.atoms:
            v = np.matrix([atom.x, atom.y, atom.z]).T
            v_new = np.dot(m,v)
            v_new = v_new.tolist()
            atom.x, atom.y, atom.z = v_new[0][0], v_new[1][0], v_new[2][0]

        return None

    ### Read and Write Methods ###
    # TODO: PDBs loaded from online are stored and loaded from tmp.

    def write_pdb(self, filename):
        "Write data to a PDB file."

        with open(filename, 'w') as f:

            # Populate 'ATOM' lines

            atom_line = ("{line_type:6}"
                         "{atom_num:>5} "
                         "{atom_name:<4}"
                         "{alt_loc:1}"
                         "{res_name:3} "
                         "{chain:1}"
                         "{res_number:4}"
                         "{icode:1}   "
                         "{x:8.3f}"
                         "{y:8.3f}"
                         "{z:8.3f}"
                         "{occupancy:6.2f}"
                         "{B_factor:6.2f}"
                         "          "
                         "{element:>2}"
                         "{charge:2}"
                         "\n")

            for count, atom in enumerate(self.atoms, 1):
                atom_parms = {'line_type':'ATOM',
                              'atom_num':count,
                              'atom_name':atom.name,
                              'alt_loc':'',
                              'res_name':atom.residue.name,
                              'chain':atom.chain,
                              'res_number':atom.residue.number,
                              'icode':'',
                              'x':atom.x,
                              'y':atom.y,
                              'z':atom.z,
                              'occupancy':1,
                              'B_factor':0,
                              'element':atom.element,
                              'charge':atom.charge}
                f.write(atom_line.format(**atom_parms))

    def read_pdb(self, filename):
        "Reads in data from a PDB file."
        
        with open(filename) as f:
            self.read_stream(f)

    def fetch_pdb(self, pdb_code):
        """Downloads/fetches a pdb file online."""
        url = 'http://ftp.rcsb.org/download/{}.pdb'.format(pdb_code)
        path = os.path.join('/tmp',pdb_code) + '.pdb'

        if not os.path.isfile(path):
            urllib.urlretrieve(url, path)
        self.read_pdb(path)
        
    def read_stream(self, stream):
        "Reads in data from a string stream."

        self.clear()

        pdb_line = re.compile((r"ATOM  (?P<number>[\s\d]{5}) "
                               "(?P<name>[\s\w]{4})"
                               "(?P<alt_loc>[\w\s])"
                               "(?P<residue_name>[\w\s]{3}) "
                               "(?P<chain>[\s\w]{1})"
                               "(?P<residue_number>[\s\w]{4})"
                               "(?P<icode>[\w\s])   "
                               "(?P<x>[\d\s\.\-]{8})"
                               "(?P<y>[\d\s\.\-]{8})"
                               "(?P<z>[\d\s\.\-]{8})"
                               "(?P<occupancy>[\d\s\.\-]{6})"
                               "(?P<B_factor>[\d\s\.\-]{6})          "
                               "(?P<element>[\s\w]{2})"
                               "(?P<charge>[\d\s\.\-]{2})?"))

        # Find the ATOM lines and pull out the necessary data
        atom_generator = ifilter(None, imap(pdb_line.match, 
                                            stream.readlines()))

        # Retrieve a set from the match objects
        for match in atom_generator:
            groupdict = {field_name:convert(field_value)
                         for field_name, field_value
                         in match.groupdict().items()}

            # create Chain, if it doesn't already exist
            id = groupdict['chain']
            if id not in self:
                chain = self.chain_class(id=id)
                chain.molecule = self
                self[id] = chain

            chain = self[id]

            # create Residue, if it doesn't already exist
            number, name = (groupdict[i] for i in ('residue_number',
                            'residue_name'))
            if number not in chain:
                try:
                    residue = self.residue_class(number=number, name=name)
                    residue.chain = chain
                    residue.molecule = self
                    chain[number] = residue
                except:
                    continue
            residue = chain[number]

            # create the Atom. The following code overwrites atoms duplicate
            # in atom name
            name = groupdict['name']
            atom_dict = {k:v for k,v in groupdict.items()
                         if k in Atom.__slots__}
            atom = self.atom_class(**atom_dict)
            atom.residue = residue
            atom.chain = chain
            atom.molecule = self
            residue[name] = atom


#### TESTS ####
import unittest

class TestMolLib(unittest.TestCase):

    pass


if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
#    unittest.main()
