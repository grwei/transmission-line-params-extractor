# Transmission Line Parameters Extractor

Matlab implementation of patent [Transmission-line simulators and methods](https://patents.google.com/patent/US8892414B1/en).

![CI](https://github.com/grwei/transmission-line-params-extractor/workflows/CI/badge.svg?branch=matlab)

## Overview

- [ ] add content

## In progress

> Need to add:
>
> - [ ]  Resonance treatment
> - [ ]  Low-frequency conversion
> - [ ]  Passivity check
> - [ ]  Causality enforcement

## Functions

- [`s2rlgc_t.m`](s2rlgc_t.m): Extract TL RLGC params from S params
- [`rlgc2s_t.m`](rlgc2s_t.m): Calculate S params from TL's RLGC params
- [`check_consistence.m`](check_consistence.m): Accept the RLGC of the transmission line as input. Call [`rlgc2s_t.m`](rlgc2s_t.m) to obtain the corresponding S parameters, and then call [`s2rlgc_t.m`](s2rlgc_t.m) to re-extract RLGC. Comparing two sets of RLGC, they should be exactly the same.

## Tests

1. [`test0.m`](test0.m): 500mil Single-line Transmission Line (Polar Si9000e)
2. [`test1.m`](test1.m): 500mil Coupled Transmission Line (Polar Si9000e)
   1. [`test1_1.m`](test1_1.m): (good) Using larger frequency intervals
   2. [`test1_2.m`](test1_2.m): (bad) Add randn noise to S, observe the impact on extracted-RLGC
   3. [`test1_3.mlx`](test1_3.mlx): (Very bad) Try to clean outliers or smooth of RLGC.
      - [使用实时编辑器任务清理杂乱数据并找到极值](https://ww2.mathworks.cn/help/matlab/data_analysis/cleandatawithliveeditortasks.html)
      - [数据平滑和离群值检测](https://ww2.mathworks.cn/help/matlab/data_analysis/data-smoothing-and-outlier-detection.html)
3. [`test2.m`](test2.m): 200mil Four-line Transmission Line (Ansys 2020R1)
   1. [~~`test2_1.m`~~](test2_1.m): (**_deprecated_**) removing singular frequencies
   2. [~~`test2_2.m`~~](test2_2.m): (**_deprecated_**) 10mm, with freq range 10M ~ 70G, HFSS `Automatically use casual materials` enabled.
   3. [`test2_3.m`](test2_3.m): HFSS_Welement -> ADS -> S
   4. [`test2_4.mlx`](test2_4.mlx): (not good) Try to clean outliers or smooth of RLGC. (Ansys 2020R1)
4. [`test3.m`](test3.m): 200mil 16-line Transmission Line (Ansys 2020R1)
   1. [`test3_1.m`](test3_1.m): HFSS_Welement -> ADS -> S
   2. [`test3_2.mlx`](test3_2.mlx): (bad) Try to clean outliers or smooth of RLGC. (Ansys 2020R1)

## Contact us

E-mail: 313017602@qq.com

## Reference

- using [Semantic Versioning](https://semver.org/)
