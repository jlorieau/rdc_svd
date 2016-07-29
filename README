---
title: Residual Dipolar Coupling SVD Fit Software for Python
author: Justin L Lorieau
tags:
    - software
---

# Overview

The rdc_svd software fits residual dipolar couplings (RDCs) to a structure from
the protein databank (PDB) using singular value decomposition (SVD).

- Github: <https://github.com/jlorieau/rdc_svd>

## Features

- Conducts an SVD fit using the same input format as DC in NMRPipe.
- Capable of dealing with sum couplings, like in methylenes, using
  the '#' atom name. (ex: 'HA#' to fit both HA1 and HA2 of glycines.)
- Capable of fitting multiple conformers from multiple PDB files
  together.
- Custom weightings of RDC values in SVD fit.
- Detailed output of the Saupe tensor and fitting statistics.

## Releases

- v1.1: 20160728.
  [v1.1.tar.gz](https://github.com/jlorieau/rdc_svd/archive/v1.1.tar.gz)
  [v1.1.zip](https://github.com/jlorieau/rdc_svd/archive/v1.1.zip)
- v1.0: 20160727

# Documentation

## Sample usage for a 3-conformer RDC dataset

~~~
$ ./rdc_svd.py -pdb test/*.pdb -dc test/hafp23_g8a.inp

Saupe(1): [ 1.24046551  0.30456653 -0.12318243  0.81267299  0.40941609]
Saupe(2): [ -7.93815366   6.16633859   6.11673696  -4.13554146 -10.0421545 ]
Saupe(3): [  5.84082969  10.23609654  -9.57914833  13.20301772   9.35618703]
S_xyz(1): [ 1.6127692  -0.54233787 -1.07043132]
S_xyz(2): [ 14.63281043  -0.13731538 -14.49549505]
S_xyz(3): [-22.06437559  17.37816129   4.6862143 ]
Da(1): 0.806Hz, Dr(1): 0.176, Rh(1): 0.218
Da(2): 7.316Hz, Dr(2): 4.786, Rh(2): 0.654
Da(3): -11.032Hz, Dr(3): 4.231, Rh(3): 0.383
Q-factor: 3.5%
RMS (Hz): 0.450861536825


================================================================================================
conformer:              atom(i) / atom(j)      D_obs (Hz)  D_pred (Hz)   Delta (Hz)        Scale
================================================================================================
test/2lwa_struc-a...       F3-N / F3-H              4.700        5.059       -0.359        1.000
test/2lwa_struc-b...       F3-N / F3-H
test/2lwa_struc-c...       F3-N / F3-H

...

test/2lwa_struc-a...     W21-CA / W21-HA            3.200        2.441        0.759        0.485
test/2lwa_struc-b...     W21-CA / W21-HA
test/2lwa_struc-c...     W21-CA / W21-HA
test/2lwa_struc-a...     Y22-CA / Y22-HA           -1.100       -2.686        1.586        0.485
test/2lwa_struc-b...     Y22-CA / Y22-HA
test/2lwa_struc-c...     Y22-CA / Y22-HA
~~~