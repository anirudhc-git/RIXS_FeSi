# RIXS_FeSi
C++ code for calculating the resonant inelastic Xray scattering (RIXS) spectrum of FeSi. Accompanies arXiv:2301.02677

Please follow these steps in order:
1. Extract **wannier90_hr.zip**

2. Compile refine.cpp as
**clang++ refine.cpp -o refine** or **g++ refine.cpp -o refine** and execute using **./refine wannier90_hr.dat refine_new.dat 52 32 729** (The usage of the program is **./refine 'In_file_name' 'out_file_name' offset no_of_bands no_of_hops**. This program can be used for other Wannier90 data files as well). This will ensure that we store the Wannier90 data in a modified format in **refine_new.dat** suitable for easy handling.

3. TODO
