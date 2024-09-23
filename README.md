# RIXS_FeSi
C++ code for calculating the resonant inelastic Xray scattering (RIXS) spectrum of FeSi. Accompanies arXiv:2301.02677

Please compile refine.cpp first as:
g++ refine.cpp -o refine

and execute using:
./refine wannier90_hr.dat refine_new.dat 52 32 729
(this is simply ./refine 'In_file_name' 'out_file_name' offset no_of_bands no_of_hops)

This will ensure that we store the Wannier90 data in a modified format suitable for easy handling.
