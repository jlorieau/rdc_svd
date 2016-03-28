---
title: MolLib for Python
author: Justin L Lorieau
tags:
    - software
---

# Overview

MolLib is a simple Python module for manipulation biomolecular structures.

- Github: <https://github.com/jlorieau/mollib>

- Download: [mollib.py](https://raw.githubusercontent.com/
jlorieau/mollib/master/mollib.py)

## Features

- Common manipulations of biomolecules for biophysics, like translation, Euler
  rotations
- *Easily extendable* data types (Atom, Residue, Chain, Molecule) through object
  inheritance to add functionality and data to base objects.
- Loading, fetching and writing PDB files

## Releases

- v2.0: 20160311
    - Added Python 3 compatibility
- v1.0: 20160309

# Documentation

## Example - Loading and Manipulating Molecules

~~~ python
>>> mol=Molecule('2OED')
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
>>> mol.write_pdb('2OED_mollib.pdb')
~~~

## Methods (iterators)

Molecule.residues
: Returns an iterator for all residues in the molecule

Molecule.atoms
: Returns an iterator for all atoms in the molecule


## Methods (accessor)

Molecule.chain_size
: Returns the number of chains in the biomolecule

Molecule.residue_size
: Returns the number of residues in the biomolecule

Molecule.atom_size
: Return the number of atoms in the biomolecule

Molecule.mass
: Return the mass (in Da) of the biomolecule

Molecule.center_of_mass
: Returns the center of mass vector of the biomolecule

## Methods (mutator)

Molecule.center()
: Centers the molecule at the center of mass.

Molecule.rotate_zyz(alpha, beta, gamma)
: Rotates the molecule using Euler Z-Y-Z angles (degrees)

## Methods (file)

Molecule.write_pdb(filename)
: Write the molecule in PDB format to the specified filename/path

Molecule.read_pdb(filename)
: Reads the PDB filename/path into the current molecule

Molecule.fetch\_pdb(pdb\_code)
: Fetches and caches (to /tmp) the PDB file for the specified pdb_code

~~~ python
>>> print "({:.3f}, {:.3f}, {:.3f})".format(*mol.center_of_mass)
(0.133, -0.347, -0.002)
~~~
