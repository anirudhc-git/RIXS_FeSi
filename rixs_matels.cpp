#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <fstream>
#include <stdlib.h>

using namespace std;

const complex<double> I(0.0,1.0);
const double Pi(acos(-1.));

int main(int argc, char *argv[])
{
    if(argc != 6)
    {
        cout<<"\n Usage: ./rixs_matels <filename> omega_in gamma eta z_renorm";

        return 0;
    }

    double omega_in = atof(argv[2]);
    double gamma = atof(argv[3]);
    double eta = atof(argv[4]);
    double z_renorm = atof(argv[5]);

    cout<<"\n omega_in = "<<omega_in<<", gamma = "<<gamma<<", eta = "<<eta<<", z_renorm = "<<z_renorm<<"\n";

    const int nos_dw_pts = 100;
    // const int nos_Gamma_pts = 11;

    // double Gamma_min = 0.1;
    // double Gamma_max = 1.1;
    double dw_max = 6.;
    // double Gamma_curr;
    double dw_curr;

    string in_file_name(argv[1]);
    string line, word;
    string out_file_name = "rixs_"+in_file_name;

    // cout<<"\n Outfile name: "<<out_file_name<<"\n";

    ifstream f1;
    f1.open(in_file_name.c_str());

    ofstream f2;
    f2.open(out_file_name.c_str(), ios::trunc);

    double E1, E2, Mats;

    double rixs_arr[nos_dw_pts];

    for(int i = 0; i < nos_dw_pts; i++)
    {
        rixs_arr[i] = 0.;
    }


    for(int i = 0; i < 3; i++)
    {
        getline(f1, line, '\r');
        
        cout<<"line "<<i+1<<" "<<line<<"\n";
        
        f2<<line<<"\n";
    }


    int half_grid_size = 24;
    double pref = 1./(Pi * pow(half_grid_size,3));

    int i, j;

    while (getline(f1, line))
    {
        getline(f1, line);

        istringstream iss(line);

        iss>>word;
        E1 = stod(word);

        iss>>word;
        E2 = stod(word);

        iss>>word;
        Mats = stod(word);

        for(i = 0; i < nos_dw_pts; i++)
        {
            dw_curr = 1. * i * dw_max / (nos_dw_pts - 1);

            rixs_arr[i] += pref * Mats * eta / (norm(omega_in - z_renorm * E2 + I * gamma) * ((z_renorm * (E2 - E1) - dw_curr)*(z_renorm * (E2 - E1) - dw_curr) + eta * eta));
        }

    }

    f2<<"# omega_in = "<<omega_in<<", gamma = "<<gamma<<", eta = "<<eta<<", z_renorm = "<<z_renorm<<" \r\n";
    f2<<"# ";

    f2<<"\r\n";
    
    for(i = 0; i < nos_dw_pts; i++)
    {
        dw_curr = 1. * i * dw_max / (nos_dw_pts - 1);

        f2<<dw_curr;

        f2<<" "<<rixs_arr[i];
        f2<<"\r\n";
    }


    f1.close();
    f2.close();

    return 0;

}