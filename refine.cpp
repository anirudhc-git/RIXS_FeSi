/* 
Program for refining the output of Wannier90. 
Author: Anirudh Chandrasekaran.
Date: 22 January 2021.

The output of this program is stored in the following format:
<row_i>   <column_j>    <R_x>   <R_y>   <R_z>   <Re[H_{ij}]>    <Im[H_{ij}]>
with iteration over R_z, followed by R_y, R_x, column_j and finally column_i. This is done keeping in mind the row major
nature of data storage in C++. While performing calculations, we would like to store this data in a buffer to speed up access.

Furthermore, since H_{ij}(k) = \sum_{R} e^{-i k.R} H_{ij}(R) and we need H_{ij}(k) = H_{ji}*(k), we need to only compute the 
upper triangular part of the matrix H_{ij}. So we store only data such that j >= i. 
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>

using namespace std;
const complex<double> I(0.0,1.0);
const int vec_dim = 3;
const double Pi(acos(-1.));

int main(int argc, char* argv[])
{
    if(argc != 6)
    {
        cout<<"\n"<<"Usage: refine 'In_file_name' 'out_file_name' offset no_of_bands no_of_hops \n";
        return 0;
    }

    string out_file_name(argv[2]);
    string in_file_name(argv[1]);

    cout<<"\n Output file name: "<<out_file_name<<"\n";
    
    int offset, hmlt_dim, old_length, new_length, no_of_hop, no_of_colums;

    // offset is the number of lines in input file where macro info is stored. Data begins at line # (offset + 1).
    // hmlt_dim is the dimension of the Hamiltonian, and also the number of bands.
    // old_length is the effective length of data (no. of lines) in input.
    // new_length is the length of the 'refined' data file.
    // no_of_hop is the number of hopping terms (or equivalently, the number of lattice sites listed in the old data).

    offset = atoi(argv[3]);
    hmlt_dim = atoi(argv[4]);
    no_of_hop = atoi(argv[5]);

    old_length = hmlt_dim * hmlt_dim * no_of_hop;

    cout<<"\n Offset: "<<offset<<",\t number of bands: "<<hmlt_dim<<",\t old length: "<<old_length<<"\n";

    no_of_colums = 7;

    new_length = ((hmlt_dim + 1) * hmlt_dim * no_of_hop)/2;
    // This is obtained as follows: there are (N^2 - N)/2 + N upper triangular entries to calculate, including diagonal
    // elements. Each requires a sum over each of the no_of_hop hopping terms with appropriate phase factor.

    int i, j, k, ind_old, ind_new;

    // Define buffer and input contents of file on to it.

    double **buff_old, **buff_new;
    buff_old = new double*[old_length];
    buff_new = new double*[new_length];

    for(i = 0; i < old_length; i++)
    {
        buff_old[i] = new double[no_of_colums];
    }

    for(i = 0; i < new_length; i++)
    {
        buff_new[i] = new double[no_of_colums];
    }


    double *evals_k;
    evals_k = new double[hmlt_dim];

    complex<double> **hmlt;
    hmlt = new complex<double>*[hmlt_dim];

    for(i = 0; i < hmlt_dim; i++)
    {
        hmlt[i] = new complex<double>[hmlt_dim];
    }

    ifstream f1;
    f1.open(in_file_name.c_str());

    ofstream f2;
    f2.open(out_file_name.c_str());

    // A temporary string and code to skip to (offset + 1)'th line.
    string temp;

    for(i = 0; i < offset; i++)
    {
        getline(f1,temp);
    }

/*  Sample code to print first line.  
    getline(f1,temp);
    cout<<"\n Contents of first (effective) line of '"<<in_file_name<<"': "<<temp<<"\n"; 
*/

    // Copy contents of the file in to the old buffer.
    for(i = 0; i < old_length; i++)
    {
        for(j = 0; j <no_of_colums; j++)
        {
            f1>>buff_old[i][j];
        }
    }

    ind_new = 0;

    // A sample code to print the (effective) first line of the data. 
    cout<<"\n ("<<buff_old[0][0]<<", "<<buff_old[0][1]<<", "<<buff_old[0][2]<<") \t ["<<buff_old[0][3]<<"]["<<buff_old[0][4]<<"] \t"<<buff_old[0][5]<<" + i"<<buff_old[0][6]<<"\n";

    for(i  = 0; i < hmlt_dim; i++)
    {
        for(j = i; j < hmlt_dim; j++)
        {
            for(k = 0; k < no_of_hop; k++)
            {
                // within the block for a particular hopping term at R, the location [i][j] is dim*j + i since i changes before j.
                // the (k + 1)'th hopping term has k terms before each corresponding to a series of i, j values of total number hmlt_dim^2
                ind_old = k * hmlt_dim * hmlt_dim + hmlt_dim * j + i;

                buff_new[ind_new][0] = buff_old[ind_old][3];
                buff_new[ind_new][1] = buff_old[ind_old][4];
                buff_new[ind_new][2] = buff_old[ind_old][0];
                buff_new[ind_new][3] = buff_old[ind_old][1];
                buff_new[ind_new][4] = buff_old[ind_old][2];
                buff_new[ind_new][5] = buff_old[ind_old][5];
                buff_new[ind_new][6] = buff_old[ind_old][6];

                ind_new++;

            }
        }
    }

    ind_new = 0;

    for(i  = 0; i < hmlt_dim; i++)
    {
        for(j = i; j < hmlt_dim; j++)
        {
            for(k = 0; k < no_of_hop; k++)
            {
                f2<<buff_new[ind_new][0]<<"  ";
                f2<<buff_new[ind_new][1]<<"    ";
                f2<<buff_new[ind_new][2]<<" ";
                f2<<buff_new[ind_new][3]<<" ";
                f2<<buff_new[ind_new][4]<<"    ";
                f2<<buff_new[ind_new][5]<<" ";
                f2<<buff_new[ind_new][6]<<"\r\n";

                ind_new++;

            }
        }
    }

    // Time for cleanup!

    f1.close();
    f2.close();

    for(i = 0; i < old_length; i++)
    {
        delete[] buff_old[i];
    }

    for(i = 0; i < new_length; i++)
    {
        delete[] buff_new[i];
    }

    for(i = 0; i < hmlt_dim; i++)
    {
        delete[] hmlt[i];
    }

    delete[] buff_new;  delete[] buff_old;  delete[] hmlt;  delete[] evals_k;


}