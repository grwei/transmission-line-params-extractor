# Transmission Line Parameters Extractor

Matlab implementation of patent [Transmission-line simulators and methods](https://patents.google.com/patent/US8892414B1/en).

![CI](https://github.com/grwei/transmission-line-params-extractor/workflows/CI/badge.svg?branch=matlab)

## Overview

- [ ] add content

## Functions

- [`s2rlgc_t.m`](s2rlgc_t.m): Extract TL RLGC params from S params
- [`rlgc2s_t.m`](rlgc2s_t.m): Calculate S params from TL's RLGC params
- [`check_consistence.m`](check_consistence.m): Accept the RLGC of the transmission line as input. Call [`rlgc2s_t.m`](rlgc2s_t.m) to obtain the corresponding S parameters, and then call [`s2rlgc_t.m`](s2rlgc_t.m) to re-extract RLGC. Comparing two sets of RLGC, they should be exactly the same.

## Tests

1. [`test0.m`](test0.m): 25.4mm Single-line Transmission Line
2. [`test1.m`](test1.m): 25.4mm Coupled Transmission Line
   1. [`test1_1.m`](test1_1.m): Using larger frequency intervals
3. [`test2.m`](test2.m): 20mm Four-line Transmission Line
   1. [`test2_1.m`](test2_1.m): removing singular frequencies
   2. [`test2_2.m`](test2_2.m): 10mm, with freq range 10M ~ 70G, HFSS `Automatically use casual materials` enabled.

## Contact me

E-mail: 313017602@qq.com
