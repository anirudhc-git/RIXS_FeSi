# RIXS_FeSi
C++ code for calculating the resonant inelastic Xray scattering (RIXS) spectrum of FeSi. Accompanies arXiv:2301.02677

Please follow these steps in order:
1. Extract **wannier90_hr.zip**.

2. Compile refine.cpp as
**clang++ refine.cpp -o refine** or **g++ refine.cpp -o refine** and execute using **./refine wannier90_hr.dat refine_new.dat 52 32 729** (The usage of the program is **./refine 'In_file_name' 'out_file_name' offset no_of_bands no_of_hops**. This program can be used for other Wannier90 data files as well). This will ensure that we store the Wannier90 data in a modified format in **refine_new.dat** suitable for easy handling.

3. The bandstructure and DOS have been plotted in **bands_DOS_final.pdf**, highlighting some possible scattering pathways from occupied to unoccupied channels.

4. Compile MatEls.cpp as **clang++ MatEls.cpp -o MatELs -llapack -lblas -Xclang -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include** (on Mac. Modify appropriately if the location of OpenMP folders are different on your local machine. You may need to install/configure it with Homebrew and Xcode). Alternatively, run **g++ MatEls.cpp -o MatELs -llapack -lblas -fopenmp**. Please ensure that **zheev.h** and **refine_new.dat** are included in the same folder while compiling/running **MatEls**.
