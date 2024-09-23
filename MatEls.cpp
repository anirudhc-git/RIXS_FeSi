/*
Program for computing the matrix elements in a specified Q-direction for FeSi.
Author: Anirudh Chandrasekaran.
Date 17 August 2023.

*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>
#include "zheev.h"
#include <omp.h>

using namespace std;

// We define some global constants here
const complex<double> I(0.0,1.0);
const int vec_dim = 3;
const double Pi(acos(-1.));
const double r2(sqrt(2.));
const double kb(8.617333e-5);
const double expt(2.71828);

// The Fermi Dirac distribution
double nF(double ee, double Tt)
{
    double tt = (1./(pow(1. * expt, 1. * ee/(1. * kb * Tt)) + 1.));
    return tt;
}

complex<double> Complex(double a, double b)
{
    return (1. * a + 1. * b * I);
}

// Number of parallel threads to be used
#define THREAD_NUM 6
// #define THREAD_NUM 24

// A function that returns the dipole matrix elements <valence,spin|e.r|core> in the form of an array
// These elements are computed in the associated mathematica file and copied here.
inline double chi(complex<double> ***Tm, complex<double> ex, complex<double> ey, complex<double> ez) // Returns polarisation matrix elements
{
    // first index - spin, second - valence, third - core.
    Tm[0][0][0] = -0.10976273017251136*(ex - Complex(0.,1.)*ey);
    Tm[0][0][1] = 1.0976273017251135*ez;
    Tm[0][1][0] = -0.3802292508725261*ez;
    Tm[0][1][1] = -0.3802292508725261*ex;
    Tm[0][2][0] = Complex(0.,0.3802292508725261)*ez;
    Tm[0][2][1] = -0.3802292508725261*ey;
    Tm[0][3][0] = 0.2851719381543946*(ex + Complex(0.,1.)*ey);
    Tm[0][3][1] = 0.;
    Tm[0][4][0] = 0.2851719381543946*(Complex(0.,-1.)*ex + ey);
    Tm[0][4][1] = 0.;
    Tm[1][0][0] = -1.0976273017251135*ez;
    Tm[1][0][1] = -0.10976273017251136*(ex + Complex(0.,1.)*ey);
    Tm[1][1][0] = 0.3802292508725261*ex;
    Tm[1][1][1] = -0.3802292508725261*ez;
    Tm[1][2][0] = 0.3802292508725261*ey;
    Tm[1][2][1] = Complex(0.,-0.3802292508725261)*ez;
    Tm[1][3][0] = 0.;
    Tm[1][3][1] = 0.2851719381543946*(ex - Complex(0.,1.)*ey);
    Tm[1][4][0] = 0.;
    Tm[1][4][1] = 0.2851719381543946*(Complex(0.,1.)*ex + ey);

    return 0;
}


// In bash do: export OMP_NUM_THREADS=5 to set number of threads.

int main(int argc, char* argv[])
{
    if(argc != 9)
    {
        cout<<"\n Usage: MatEls  'in_file'  'out_file'  no_of_bands  no_of_hops  half_grid_size  TTH  TH omega_in_eV";

        return 0;
    }


    // Storing input and output filenames on to the respective strings:
    string out_file_name(argv[2]);
    string in_file_name(argv[1]);

    // Defining some essential variables, some of which are passed as command line arguments.
    // hmlt_dim is the number of bands (not including spin degeneracy). This is same as the dimension of the k-space Hamiltonian.
    // no_of_hops is the number of hopping terms stored in the Wannier90 output file. 
    // size of grid of k points used in the discretized Brillouin zone is (2 * half_grid_size + 1)^3.
    // Gamma_g is the electron inverse life time that we want to use for computing the absorption spectrum.
    // data_length and no_of_colums will be explained below...
    int hmlt_dim, no_of_hops, no_of_plt_pts, data_length, no_of_columns, half_grid_size;
    double epsilon;

    // double Gamma_g;

    hmlt_dim = atoi(argv[3]);
    no_of_hops = atoi(argv[4]);
    // no_of_plt_pts = atoi(argv[5]);
    half_grid_size = atoi(argv[5]);
    // Gamma_g = atof(argv[7]);

    // th_i and two_th are theta and two_theta which describe the momentum transfer direction.
    double th_i, two_th, omega_in;

    two_th = atof(argv[6]);
    th_i = atof(argv[7]);
    omega_in = atof(argv[8]);

    // double omega_g;
    // omega_g = atof(argv[10]); // = 0.5 for example.

    no_of_columns = 7; // 7 columns in modified data are i, j, r_x, r_y, r_z, Re(H_ij), Im(H_ij)

    data_length = ((hmlt_dim + 1) * hmlt_dim * no_of_hops)/2;

    // Mu is the chemical potential offset.
    double Mu;

    // a data buffer to store the hopping terms that are obtianed from in_file, which is itself obtained from the output of Wannier90.
    double **data_buff;
    data_buff = new double*[data_length];

    for(int i = 0; i < data_length; i++)
    {
        data_buff[i] = new double[no_of_columns];
    }

    ifstream f1;
    f1.open(in_file_name.c_str());

    ofstream f2;
    f2.open(out_file_name.c_str(),ios::trunc);

    for(int i = 0; i < data_length; i++)
    {
        for(int j = 0; j < no_of_columns; j++)
        {
            f1>>data_buff[i][j];
        }
    }

    // To store the Rixs vs omega_in data:
    // double *RIXS;
    // RIXS = new double[no_of_plt_pts];

    // This sets Mu to be exactly in the gap in the middle of the bands as it is expected for FeSi. See the band plots for details.
    Mu = 9.55;

    // eta is the broadening factor
    // double eta;

    // eta = 0.05; // 50 meV

    /*
    for(int i = 0; i < no_of_plt_pts; i++)
    {
        RIXS[i] = 0.;
    }
    */

    int TOT_THRDS = THREAD_NUM;

    // double w_max = 5.; // in eV.
    double q_vec[3];

    // Ignore the following, commented out code block.
    // epsilon = 0.16235;
/*    
    q_vec[0] = 2. * Pi * epsilon / sqrt(2.);
    q_vec[1] = 0.;
    q_vec[2] = 2. * Pi * epsilon / sqrt(2.);
*/    
    // q_vec[0] = 2. * Pi / sqrt(3.);
    // q_vec[1] = 2. * Pi / sqrt(3.);
    // q_vec[2] = 2. * Pi / sqrt(3.);
/*
    q_vec[0] = 2. * Pi * 0.06;
    q_vec[1] = 0.;
    q_vec[2] = 2. * Pi * 0.49;

    complex<double> ex_in_g, ey_in_g, ez_in_g, ex_out_g, ey_out_g, ez_out_g;

    ex_in_g = -cos(th_i); // TODO
    ey_in_g = sin(th_i); // TODO
    ez_in_g = 0.; // TODO

    ex_out_g = -cos(two_th - th_i);
    ey_out_g = -sin(two_th - th_i);
    ez_out_g = 0.;

    ex_in_g = -sin(th_i); // TODO
    ez_in_g = cos(th_i); // TODO
    ey_in_g = 0.; // TODO
    

    ex_out_g = sin(two_th - th_i);
    ey_out_g = 0.;
    ez_out_g = cos(two_th - th_i);
*/

    double pimag; // This is the magnitude of the ingoing momentum. This is calculated as follows:
    // p = E/c for light.
    // pimag = ((ingoing photon energy in eV / charge of e-) / speed of light in m/s) / (2 * Pi * hbar / lattice spacing)
    // Now pimag is in the units of length of the Brillouoin zone (which equals 2 * Pi * hbar / lattice spacing)
    pimag = ((omega_in) * 1.602e-19 / 3e8) / 1.4962e-24;

    // for defining q vectors for calculation we use extra factors of 2*Pi since we want the functions defining the band hamiltonian (containing sines and cosines) to be periodic over one BZ
    q_vec[0] = 2. * Pi * pimag * (cos(th_i * Pi / 180.) - cos((two_th - th_i) * Pi / 180.));
    q_vec[1] = 0.;
    q_vec[2] = 2. * Pi * pimag * (sin(th_i * Pi / 180.) + sin((two_th - th_i) * Pi / 180.));

    complex<double> ex_in_g, ey_in_g, ez_in_g, ex_out_g, ey_out_g, ez_out_g;
    
    // we now define the ingoing and outgoing polarization vectors
    ex_in_g = -sin(th_i * Pi / 180.); // TODO
    ez_in_g = cos(th_i * Pi / 180.); // TODO
    ey_in_g = 0.; // TODO
    

    ex_out_g = sin((two_th - th_i) * Pi / 180.);
    ey_out_g = 0.;
    ez_out_g = cos((two_th - th_i) * Pi / 180.);

    // Text header of the output file:
    f2<<"# Numerically calculated RIXS matrix elements on the full BZ with grid size "<<(2*half_grid_size+1)<<"^3"<<", Q_dir = 2*Pi*("<<(q_vec[0]/(2. * Pi))<<", "<<(q_vec[1]/(2. * Pi))<<", "<<(q_vec[2]/(2. * Pi))<<")\r\n";
    f2<<"# Ingoing polarization (LH) = ("<<ex_in_g<<", "<<ey_in_g<<", "<<ez_in_g<<")\r\n";
    f2<<"# TTH = "<<two_th<<", TH = "<<th_i<<"\r\n";
    /*
    q_vec[0] = epsilon * q_vec[0];
    q_vec[1] = epsilon * q_vec[1];
    q_vec[2] = epsilon * q_vec[2];
    */

    omp_set_dynamic(0);
    #pragma omp parallel num_threads(THREAD_NUM)              
    {
   
        int threads = THREAD_NUM;

        // we now define local variables for each thread and equate them to the global ones in main.
        // This ensures that the threads can run independently
        // double Gamma; // core hole life time

        // double Tt; // Temperature

        // Tt = 75.;

        // We split the Brillouin zone into equal segments along the k_x direction. Each thread goes through points on one of the segments.
        int interval = (int)(2 * half_grid_size/threads);

        // cout<<"\n Number of points in each thread = "<<interval;

        // In any particular thread, k_x takes values between k1_min and k1_max, that define the particular segment of the BZ.
        // other variable will be defined below
        int k1_min, k1_max, k1, k2, k3, j, l, m, n, i, old_ind, u, no_val_orbs, no_cor_orbs;
        int la, la_p;

        // variables to be used as appropriate indices for summing over.
        int val1, cor, val2;

        double temp, del_omega, omega, Mu_loc;

        // local copies in each thread
        // omega = omega_g; 
        // Gamma = Gamma_g; 
        Mu_loc = Mu;

        complex<double> ex_in, ey_in, ez_in, ex_out, ey_out, ez_out;
        ex_in = ez_in_g;
        ey_in = ex_in_g;
        ez_in = ey_in_g;

        ex_out = ez_out_g;
        ey_out = ex_out_g;
        ez_out = ey_out_g;

        // We define two arrays to store the band eigenvalues at k and k + q.
        double *evals_kq;
        evals_kq = new double[hmlt_dim];

        double *evals_k;
        evals_k = new double[hmlt_dim];

        // Defining arrays for storing Hamiltonian and eigenvectors
        complex<double> **hmlt, **hmlt2, **evec_k, **evec_kq;
        hmlt = new complex<double>*[hmlt_dim];
        hmlt2 = new complex<double>*[hmlt_dim];
        evec_k = new complex<double>*[hmlt_dim];
        evec_kq = new complex<double>*[hmlt_dim];

        complex<double> **scat_mat;
        scat_mat = new complex<double>*[hmlt_dim];

        for(i = 0; i < hmlt_dim; i++)
        {
            hmlt[i] = new complex<double>[hmlt_dim];
            hmlt2[i] = new complex<double>[hmlt_dim];
            scat_mat[i] = new complex<double>[hmlt_dim];
            evec_k[i] = new complex<double>[hmlt_dim];
            evec_kq[i] = new complex<double>[hmlt_dim];
        }

        double k_vec[3], r_vec[3], kq_vec[3], q_vec_local[3], w_max_local;

        // local copy of q vector for each thread.
        q_vec_local[0] = q_vec[2];
        q_vec_local[1] = q_vec[0];
        q_vec_local[2] = q_vec[1];

        // w_max_local = w_max;

        // The k1 sum is to be parallelized

        int cur_thread = omp_get_thread_num();

        k1_min = -half_grid_size + cur_thread * interval;
        k1_max = -half_grid_size + (cur_thread + 1) * interval;

        int half_grid_size_local = half_grid_size;
        int hmlt_dim_local = hmlt_dim;
        // int no_of_plt_pts_local = no_of_plt_pts;

        // double *RIXS_local;
        // RIXS_local = new double[no_of_plt_pts_local];
        complex<double> ctemp1, ctemp2;

/*
        for(i = 0; i < no_of_plt_pts_local; i++)
        {
    
            RIXS_local[i] = 0.;
        }
*/

        bool scat;
        int no_metal_atms;

        no_metal_atms = 4;

        no_val_orbs = 5; // spin is treated separately
        no_cor_orbs = 2;

        complex<double> ***Tm_in = new complex<double>**[2]; // two spin species
        for(i = 0; i < 2; i++)
        {
            Tm_in[i] = new complex<double>*[no_val_orbs];

            for(j = 0; j < no_val_orbs; j++)
            {
                Tm_in[i][j] = new complex<double>[no_cor_orbs];
            }
        }

        complex<double> ***Tm_out_LH = new complex<double>**[2]; // two spin species
        for(i = 0; i < 2; i++)
        {
            Tm_out_LH[i] = new complex<double>*[no_val_orbs];

            for(j = 0; j < no_val_orbs; j++)
            {
                Tm_out_LH[i][j] = new complex<double>[no_cor_orbs];
            }
        }

        complex<double> ***Tm_out_LV = new complex<double>**[2]; // two spin species
        for(i = 0; i < 2; i++)
        {
            Tm_out_LV[i] = new complex<double>*[no_val_orbs];

            for(j = 0; j < no_val_orbs; j++)
            {
                Tm_out_LV[i][j] = new complex<double>[no_cor_orbs];
            }
        }

        chi(Tm_out_LH, ex_out, ey_out, ez_out); // storing the matrix elements of dipole operator
        chi(Tm_out_LV, 0., 0., 1.); // For outgoing photon, we sum over horizontal and vertical polarizations.
        chi(Tm_in, ex_in, ey_in, ez_in);

// **************************************************************************

        for(k1 = k1_min; k1 < k1_max; k1++)
        {
            k_vec[0] = (Pi * 1. * k1)/(1. * half_grid_size_local);

            for(k2 = -half_grid_size_local; k2 < half_grid_size_local; k2++)
            {
                k_vec[1] = (Pi * 1. * k2)/(1. * half_grid_size_local);
                
                for(k3 = -half_grid_size_local; k3 < half_grid_size_local; k3++)
                {
                    k_vec[2] = (Pi * 1. * k3)/(1. * half_grid_size_local);

                    
                    kq_vec[0] = k_vec[0] + q_vec_local[0];
                    kq_vec[1] = k_vec[1] + q_vec_local[1];
                    kq_vec[2] = k_vec[2] + q_vec_local[2];
                    

                    // Now we define the Hamiltonian at k and compute its eigenvalues.
                    // We use the hopping terms that we loaded from the in_file
                    // This scheme for defining the Hamiltonian seems to be correct, as reflected by the band structure generated using this.
                    old_ind = 0;

                    for(j = 0; j < hmlt_dim_local; j++)
                    {
                        for(l = j; l < hmlt_dim_local; l++)
                        {
                            hmlt[j][l] = 0. + 0. * I;
                            hmlt2[j][l] = 0. + 0. * I;

                            for(m = 0; m < no_of_hops; m++)
                            {
                                // old_ind is the 'line' in the file. 2 reflects the fact that position of lattice site starts at (2+1)th column.
                                r_vec[0] = data_buff[old_ind][2];
                                r_vec[1] = data_buff[old_ind][3];
                                r_vec[2] = data_buff[old_ind][4];

                                hmlt[j][l] = hmlt[j][l] + (cos(k_vec[0] * r_vec[0] + k_vec[1] * r_vec[1] + k_vec[2] * r_vec[2]) - I * sin(k_vec[0] * r_vec[0] + k_vec[1] * r_vec[1] + k_vec[2] * r_vec[2])) * (data_buff[old_ind][5] + I * data_buff[old_ind][6]);
                                hmlt2[j][l] = hmlt2[j][l] + (cos(kq_vec[0] * r_vec[0] + kq_vec[1] * r_vec[1] + kq_vec[2] * r_vec[2]) - I * sin(kq_vec[0] * r_vec[0] + kq_vec[1] * r_vec[1] + kq_vec[2] * r_vec[2])) * (data_buff[old_ind][5] + I * data_buff[old_ind][6]);
                                
                                old_ind++;

                            }

                            hmlt[l][j] = conj(hmlt[j][l]);
                            hmlt2[l][j] = conj(hmlt2[j][l]);
                        }
                    }

                    zheev(hmlt, hmlt_dim_local, evals_k, evec_k);

                    zheev(hmlt2, hmlt_dim_local, evals_kq, evec_kq);

                    /*
                    if(k1 == 1 && k2 == 1 && k3 == 1)
                    {
                        cout<<"\n Done! "<<evals_k[2]<<" and "<<evals_k[28]<<"\n";
                    }
                    */

                    scat = false;

                    // la and la_p are band indices lambda and lambda'                    
                    
                        
                    for(la = 0; la < hmlt_dim_local; la++)
                    {
                        // if we had only one metal atom per unit cell, we wouldn't need the l sum,...
                        // or if there were no non-metal orbitals participating in the TBM as valence sum would...
                        // take care of all the metal orbitals.
                    
                        for(la_p = 0; la_p < hmlt_dim_local; la_p++)
                        {
                            temp = 0.;    
                            ctemp1 = 0.;
                            ctemp2 = 0.;

                            if(evals_k[la_p] < Mu_loc && evals_kq[la] >= Mu_loc) // If scattering is possible
                            {
                                scat = true;

                                for(l = 0; l < no_metal_atms; l++)
                                {
                                    for(val1 = 0; val1 < no_val_orbs; val1++)
                                    {
                                        for(val2 = 0; val2 < no_val_orbs; val2++)
                                        {
                                            for(cor = 0; cor < no_cor_orbs; cor++)
                                            {
                                                // only one core index to sum since same core orbital. But two valence summation,
                                                // one for ingoing and one for outgoing.

                                                // also ctemp1 is for Linear horizontal outgoing polarization, and ctemp2 is linear vertical
                                                ctemp1 = ctemp1 + (Tm_in[0][val1][cor] + Tm_in[1][val1][cor]) * conj(Tm_out_LH[0][val2][cor] + Tm_out_LH[1][val2][cor]) * evec_k[val2 + 5 * l][la_p] * conj(evec_kq[val1 + 5 * l][la]);
                                                ctemp2 = ctemp2 + (Tm_in[0][val1][cor] + Tm_in[1][val1][cor]) * conj(Tm_out_LV[0][val2][cor] + Tm_out_LV[1][val2][cor]) * evec_k[val2 + 5 * l][la_p] * conj(evec_kq[val1 + 5 * l][la]);
                                            }
                                        }
                                    }                            
                                }
                            
                                ctemp1 = 1. * ctemp1;
                                ctemp2 = 1. * ctemp2;

                                temp = temp + norm(ctemp1) + norm(ctemp2);

                                #pragma omp critical
                                {
                                    f2<<(evals_k[la_p]-Mu_loc)<<" "<<(evals_kq[la]-Mu_loc)<<" "<<temp<<"\r\n";
                                }
                                

                                //  * (1. - nF(evals_kq[la], Tt)) * nF(evals_k[la_p], Tt) within the norm for T dependence.
                                // scat_mat[la][la_p] = 1./(omega - (evals_kq[la] - Mu) + I * Gamma);  

                                        // Here we add the result of each of the threads on to the RIXS absorption spectrum.
                                /*                                #pragma omp critical
                                {
                                    
                                    for(i = 0; i < no_of_plt_pts_local; i++)
                                    {
                                        RIXS[i] = RIXS[i] + RIXS_local[i];
                                    }
                                    

                                    
                                }
                                */
                      
                            }
                        }

                    }


                    /*
                    if(scat == true)
                    {
                        for(n = 0; n < no_of_plt_pts_local; n++)
                        {
                            temp = 0.;
                            del_omega = 1. * n * w_max/ no_of_plt_pts_local;
                            
                            for(la = 0; la < hmlt_dim_local; la++)
                            {
                                for(la_p = 0; la_p < hmlt_dim_local; la_p++)
                                {
                                        temp = temp + 1. * norm(scat_mat[la][la_p]) / ((del_omega - (evals_kq[la] - evals_k[la_p]))*(del_omega - (evals_kq[la] - evals_k[la_p])) + eta * eta);
                                        // there is a gamma/pi overall factor but we skip it since the units are arbitrary.
                                }
                            }
                            
                                                                              
                            RIXS_local[n] = RIXS_local[n] + temp;
                        }

                    }

                    */

                }
            }
        }


        for(i = 0; i < 2; i++)
        {
            for(j = 0; j < no_val_orbs; j++)
            {
                delete[] Tm_in[i][j];
                delete[] Tm_out_LH[i][j];
                // delete[] Tm_lin[i][j];
                delete[] Tm_out_LV[i][j];
            }
        }

        for(i = 0; i < 2; i++)
        {
            delete[] Tm_in[i];
            delete[] Tm_out_LH[i];
            // delete[] Tm_lin[i];
            delete[] Tm_out_LV[i];
        }

        for(i = 0; i < hmlt_dim; i++)
        {
            delete[] hmlt[i];
            delete[] hmlt2[i];
            delete[] scat_mat[i];
            delete[] evec_k[i];
            delete[] evec_kq[i];
        }

        delete[] hmlt;  delete[] hmlt2;  delete[] evals_k;  delete[] scat_mat; //   delete[] RIXS_local;
        delete[] evec_k;    delete[] evec_kq;
        delete[] Tm_in; delete[] Tm_out_LH; // delete[] Tm_lin; 
        delete[] Tm_out_LV; 

    }

    // pref ensures that the results don't depend on number of grid points
    // double pref;

    // pref = 1./(Pi * pow(half_grid_size,3));


    // double omg;

    /*
    for(int i = 0; i < no_of_plt_pts; i++)
    {
        omg = -w_max + 2. * i * w_max/ no_of_plt_pts;
        RIXS[i] = RIXS[i] * pref;

        f2<<omg<<" "<<RIXS[i]<<"\r\n";
    }
    */


    // Time for cleanup!

    f1.close();
    f2.close();

    for(int i = 0; i < data_length; i++)
    {
        delete[] data_buff[i];
    }


    delete[] data_buff;  // delete[] RIXS;

    return 0;
}


/*
    Tm[0][0][0] = -0.109763 * (ex - (0. + 1. * I) * ey);
    Tm[0][0][1] = 1.09763 * ez;
    Tm[0][1][0] = -0.380229 * ez;
    Tm[0][1][1] = -0.380229 * ex;
    Tm[0][2][0] = (0. + 0.380229 * I) * ez;
    Tm[0][2][1] = -0.380229 * ey;
    Tm[0][3][0] = 0.285172 * (ex + (0. + 1. * I) * ey);
    Tm[0][3][1] = 0.;
    Tm[0][4][0] = 0.285172 * ((0. - 1. * I) * ex + ey);
    Tm[0][4][1] = 0.;

    Tm[1][0][0] = -1.09763 * ez;
    Tm[1][0][1] = -0.109763 * (ex + (0. + 1. * I) * ey);
    Tm[1][1][0] = 0.380229 * ex;
    Tm[1][1][1] = -0.380229 * ez;
    Tm[1][2][0] = 0.380229 * ey;
    Tm[1][2][1] = (0. - 0.380229 * I) * ez;
    Tm[1][3][0] = 0.;
    Tm[1][3][1] = 0.285172 * (ex - (0. + 1. * I) * ey);
    Tm[1][4][0] = 0.;
    Tm[1][4][1] = 0.285172 * ((0. + 1. * I) * ex + ey);

    */
