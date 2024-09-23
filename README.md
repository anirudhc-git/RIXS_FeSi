# RIXS_FeSi
C++ code for calculating the resonant inelastic Xray scattering (RIXS) spectrum of FeSi. Accompanies arXiv:2301.02677

Please follow these steps in order:
1. Extract **wannier90_hr.zip**.

2. Compile refine.cpp as
**clang++ refine.cpp -o refine** or **g++ refine.cpp -o refine** and execute using **./refine wannier90_hr.dat refine_new.dat 52 32 729** (The usage of the program is **./refine 'In_file_name' 'out_file_name' offset no_of_bands no_of_hops**. This program can be used for other Wannier90 data files as well). This will ensure that we store the Wannier90 data in a modified format in **refine_new.dat** suitable for easy handling.

3. The bandstructure and DOS have been plotted in **bands_DOS_final.pdf**, highlighting some possible scattering pathways from occupied to unoccupied channels.

4. Compile MatEls.cpp as **clang++ MatEls.cpp -o MatELs -llapack -lblas -Xclang -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include** (on Mac. Modify appropriately if the location of OpenMP folders are different on your local machine. You may need to install/configure it with Homebrew and Xcode). Alternatively, run **g++ MatEls.cpp -o MatELs -llapack -lblas -fopenmp**. Please ensure that **zheev.h** and **refine_new.dat** are included in the same folder while compiling/running **MatEls**.

5. To run the code for a single input case **./MatEls refine_new.dat matsTTH150TH120.dat 32 729 24 150 120 708.8**. Include an & at the end to run in the background. It took about 24 minutes on my Macbook Pro (M1). If your computer supports more threads, change THREAD_NUM from 6 to higher on line 40 of MatEls.cpp. For running a batch together **(./MatEls refine_new.dat matsTTH150TH10.dat 32 729 24 150 10 708.8;  ./MatEls refine_new.dat matsTTH150TH30.dat 32 729 24 150 30 708.8; ./MatEls refine_new.dat matsTTH150TH45.dat 32 729 24 150 45 708.8; ./MatEls refine_new.dat matsTTH150TH68.dat 32 729 24 150 68 708.8; ./MatEls refine_new.dat matsTTH150TH120.dat 32 729 24 150 120 708.8; ./MatEls refine_new.dat matsTTH70TH10.dat 32 729 24 70 10 708.8; ./MatEls refine_new.dat matsTTH70TH30.dat 32 729 24 70 30 708.7; ./MatEls refine_new.dat matsTTH70TH60.dat 32 729 24 70 60 708.6)&**

6. The data generated contains the appropriate occupied/unoccupied states, the energy levels and appropriate transition matrix elements needed for the RIXS calculation. Now compile the program to compute the RIXS spectrum from the generated data on Mac as **clang++ rixs_matels.cpp -o rixs_matels** (or replace clang++ with g++ for GCC). Try executing ./rixs_matels to see the usage.

7. To generate the RIXS data for 2θ = 150°, execute **(./rixs_matels matsTTH150TH10.dat 1.13 0.8 0.075 1.; ./rixs_matels matsTTH150TH30.dat 1.13 0.8 0.075 1.; ./rixs_matels matsTTH150TH45.dat 1.13 0.8 0.075 1.; ./rixs_matels matsTTH150TH68.dat 1.13 0.8 0.075 1.; ./rixs_matels matsTTH150TH120.dat 1.13 0.8 0.075 1.)&**

8. For 2θ = 70°, we need to modify the the inputs since omega_in is slightly different
**(./rixs_matels matsTTH70TH10.dat 1.13 0.8 0.075 1.; ./rixs_matels matsTTH70TH30.dat 1.03 0.8 0.075 1.; ./rixs_matels matsTTH70TH60.dat 0.93 0.8 0.075 1.)&**
