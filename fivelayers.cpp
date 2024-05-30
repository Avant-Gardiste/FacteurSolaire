#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include<chrono>
#include "constantes.h"

using namespace Eigen;
using namespace std;

enum Gaz{AIR, ARGON, KRYPTON, XENON, NUM_GAZ};

vector<double> calculHg(double T, double dT, double s, const std::vector<double>& fractions, int inclinaison)
{   
    double Gr = 0;
    double Pr = 0;
    double Nu = 0;
    double hg = 0;
    int ipos = 1;

    (inclinaison < 60) ? ipos = 2 : ipos;
    (inclinaison < 30) ? ipos = 3 : ipos;

    // Propriétés de base pour chaque gaz
    std::vector<double> rho(NUM_GAZ), muy(NUM_GAZ), lambda(NUM_GAZ), c(NUM_GAZ);
    double final_rho = 0, final_muy = 0, final_lambda = 0, final_c = 0;
    
    // Voir Annexe F (page 43) - NF EN ISO 52022-3
    rho[AIR] = 1.189 - 0.0044 * (T - 293.15); // 293.15 correspond à 20°C en Kelvin 
    muy[AIR] = 0.00001811 + 0.00000005 * (T - 293.15);
    lambda[AIR] = 0.02576 + 0.00008 * (T - 293.15);
    c[AIR] = 1008.0;

    rho[ARGON] = 1.640 - 0.006 * (T - 293.15);
    muy[ARGON] = 0.00002228 + 0.000000064 * (T - 293.15);
    lambda[ARGON] = 0.01734 + 0.00005 * (T - 293.15);
    c[ARGON] = 519.0;

    rho[KRYPTON] = 3.43 - 0.013 * (T - 293.15);
    muy[KRYPTON] = 0.0000247 + 0.00000007 * (T - 293.15);
    lambda[KRYPTON] = 0.00926 + 0.000026 * (T - 293.15);
    c[KRYPTON] = 245.0;

    rho[XENON] = 5.495 - 0.0209 * (T - 293.15);
    muy[XENON] = 0.00002299 + 0.000000074 * (T - 293.15);
    lambda[XENON] = 0.00546 + 0.000017 * (T - 293.15);
    c[XENON] = 161.0;

    // Calcul des propriétés finales du mélange
    for (int i = 0; i < NUM_GAZ; i++) {
        final_rho += rho[i] * fractions[i];
        final_muy += muy[i] * fractions[i];
        final_lambda += lambda[i] * fractions[i];
        final_c += c[i] * fractions[i];
    }
 
    // Calcul de Grashof, Prandtl, et Nusselt selon la norme NF EN 673 (Page 7-8 dans la norme NF EN 673)
    Gr = 9.81 * pow(s, 3) * dT * pow(final_rho, 2) / (T * pow(final_muy, 2)); // Equation (7) page 8 (NF EN 673)
    Pr = final_muy * final_c / final_lambda; // Equation (8) page 8 (NF EN 673)
    
    (ipos == 1) ? Nu = fmax(0.035 * pow(Gr * Pr, 0.38), 1) : Nu;  // Equation (6) page 8 (NF EN 673) (A = 0.035, n=0.38)
    (ipos == 2) ? Nu = fmax(0.16 * pow(Gr * Pr, 0.28), 1) : Nu; 
    (ipos == 3) ? Nu = fmax(0.1 * pow(Gr * Pr, 0.31), 1) : Nu; 

    // Calcul de hg
    hg = Nu * final_lambda / s;  // Equation (5) page 7 (NF EN 673)

    vector<double> results;
    results.push_back(hg);   
    results.push_back(final_rho); 
    results.push_back(final_muy);  
    results.push_back(final_c); 
    return results;
}

double ems_corrige(double ems)
{   
    double ems_cor = 1.1887 * ems - 0.4967 * ems * ems + 0.2452 * ems * ems * ems;
    return ems_cor;
}

MatrixXd CreateTemperatureMatrix(double Te, double Ti, double hc_e, double hc_i, vector<double>ems, vector<double>ems_prime, vector<double>s){
    int nombre_couches = ems.size();
    MatrixXd temperatureMatrix = MatrixXd::Zero(nombre_couches, nombre_couches);

    double deltaT;
    double he, hi;
    double hr_1_2, hr_2_1, hr_2_3, hr_3_2, hr_3_4, hr_4_3, hr_4_5, hr_5_4;
    double hg_1, hg_2, hg_3, hg_4;
    double s_1, s_2, s_3, s_4;

    he = hc_e + 4.0*_StefanBoltzmann*pow(Te,3);
    hi = hc_i + 4.0*_StefanBoltzmann*pow(Ti, 3);

    if(nombre_couches==2)
    {   
        deltaT = 5.0;
        s_1 = s[0];

        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hr_1_2 = 4.0 * (1 / (1/ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_1 = 4.0 * (1 / (1/ems_prime[1] + 1/ems[0] - 1)) * _StefanBoltzmann * pow(_T0,3);

        temperatureMatrix(0,0) = -(he+hr_1_2+hg_1);
        temperatureMatrix(1,0) = hr_2_1;
        temperatureMatrix(0,1) = hr_1_2;
        temperatureMatrix(1,1) = -(hr_2_1+hi+hg_1);
    }

    else if (nombre_couches==3)
    {   
        deltaT = 2.5;
        s_1 = s[0];
        s_2 = s[1];
        
        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hg_2 = calculHg(_T0, deltaT, s_2, {1, 0, 0, 0}, 90)[0];
        hr_1_2 = 4.0 * (1 / (1/ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_3 = 4.0 * (1 / (1/ems_prime[1] + 1/ems[2] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_3_2 = 4.0 * (1 / (1/ems_prime[2] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);

        temperatureMatrix(0,0) = -(he+hr_1_2+hg_1);
        temperatureMatrix(1,0) = hr_1_2;
        temperatureMatrix(2,0) = 0.0;

        temperatureMatrix(0,1) = hr_1_2;
        temperatureMatrix(1,1) = -(hr_1_2+hr_2_3+hg_1+hg_2);
        temperatureMatrix(2,1) = hr_3_2;

        temperatureMatrix(0,2) = 0.0;
        temperatureMatrix(1,2) = hr_2_3;
        temperatureMatrix(2,2) = -(hi+hr_3_2+hg_2);
    }

    else if(nombre_couches==4){
        deltaT = 2.5;
        s_1 = s[0];
        s_2 = s[1];
        s_3 = s[2];

        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hg_2 = calculHg(_T0, deltaT, s_2, {1, 0, 0, 0}, 90)[0];
        hg_3 = calculHg(_T0, deltaT, s_3, {1, 0, 0, 0}, 90)[0];

        hr_1_2 = 4.0 * (1 / (1 / ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_1 = 4.0 * (1 / (1 / ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_3 = 4.0 * (1 / (1 / ems_prime[1] + 1/ems[2] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_3_2 = 4.0 * (1 / (1 / ems_prime[1] + 1/ems[2] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_3_4 = 4.0 * (1 / (1 / ems_prime[2] + 1/ems[3] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_4_3 = 4.0 * (1 / (1 / ems_prime[2] + 1/ems[3] - 1)) * _StefanBoltzmann * pow(_T0,3);
        
        temperatureMatrix(0,0) = -(he+hr_1_2+hg_1);
        temperatureMatrix(1,0) = hr_2_1;
        temperatureMatrix(2,0) = 0;
        temperatureMatrix(3,0) = 0;

        temperatureMatrix(0,1) = hr_1_2;
        temperatureMatrix(1,1) = -(hr_1_2+hr_2_3+hg_1+hg_2);
        temperatureMatrix(2,1) = hr_3_2;
        temperatureMatrix(3,1) = 0;

        temperatureMatrix(0,2) = 0;
        temperatureMatrix(1,2) = hr_2_3;
        temperatureMatrix(2,2) = -(hr_3_2+hr_3_4+hg_2+hg_3);
        temperatureMatrix(3,2) = hr_4_3;

        temperatureMatrix(0,3) = 0;
        temperatureMatrix(1,3) = 0;
        temperatureMatrix(2,3) = hr_3_4;
        temperatureMatrix(3,3) = -(hi+hr_4_3+hg_3);
    }

    else if(nombre_couches == 5){
        deltaT = 2.5;
        s_1 = s[0];
        s_2 = s[1];
        s_3 = s[2];
        s_4 = s[3];

        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hg_2 = calculHg(_T0, deltaT, s_2, {1, 0, 0, 0}, 90)[0];
        hg_3 = calculHg(_T0, deltaT, s_3, {1, 0, 0, 0}, 90)[0];
        hg_4 = calculHg(_T0, deltaT, s_4, {1, 0, 0, 0}, 90)[0];

        hr_1_2 = 4.0 * (1 / (1 / ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_1 = 4.0 * (1 / (1 / ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_3 = 4.0 * (1 / (1 / ems_prime[1] + 1/ems[2] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_3_2 = 4.0 * (1 / (1 / ems_prime[1] + 1/ems[2] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_3_4 = 4.0 * (1 / (1 / ems_prime[2] + 1/ems[3] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_4_3 = 4.0 * (1 / (1 / ems_prime[2] + 1/ems[3] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_4_5 = 4.0 * (1 / (1 / ems_prime[3] + 1/ems[4] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_5_4 = 4.0 * (1 / (1 / ems_prime[3] + 1/ems[4] - 1)) * _StefanBoltzmann * pow(_T0,3);
        
        temperatureMatrix(0,0) = -(he+hr_1_2+hg_1);
        temperatureMatrix(1,0) = hr_2_1;
        temperatureMatrix(2,0) = 0;
        temperatureMatrix(3,0) = 0;
        temperatureMatrix(4,0) = 0;

        temperatureMatrix(0,1) = hr_1_2;
        temperatureMatrix(1,1) = -(hr_1_2+hr_2_3+hg_1+hg_2);
        temperatureMatrix(2,1) = hr_3_2;
        temperatureMatrix(3,1) = 0;
        temperatureMatrix(4,0) = 0;

        temperatureMatrix(0,2) = 0;
        temperatureMatrix(1,2) = hr_2_3;
        temperatureMatrix(2,2) = -(hr_3_2+hr_3_4+hg_2+hg_3);
        temperatureMatrix(3,2) = hr_4_3;
        temperatureMatrix(4,0) = 0;

        temperatureMatrix(0,3) = 0;
        temperatureMatrix(1,3) = 0;
        temperatureMatrix(2,3) = hr_3_4;
        temperatureMatrix(3,3) = -(hr_4_3+hr_5_4+hg_3+hg_4);
        temperatureMatrix(4,0) = hr_5_4;

        temperatureMatrix(0,4) = 0;
        temperatureMatrix(1,4) = 0;
        temperatureMatrix(2,4) = 0;
        temperatureMatrix(3,4) = hr_4_5;
        temperatureMatrix(4,4) = -(hi+hr_5_4+hg_4);

    }
    
    return temperatureMatrix;
}

VectorXd createTemperatureVector(double Es, double hc_e, double hc_i, double Te, double Ti, vector<double> alphaE, vector<double> s){
    int nombre_couches = alphaE.size();
    VectorXd temperatureVect(nombre_couches);

    double deltaT;
    double he, hi;
    double hg_1, hg_2, hg_3, hg_4;
    double s_1 = s[0];
    double s_2 = s[1];
    double s_3 = s[2];
    double s_4 = s[3];

    he = hc_e + 4.0*_StefanBoltzmann*pow(Te,3);
    hi = hc_i + 4.0*_StefanBoltzmann*pow(Ti, 3);

    if(nombre_couches==2){
        deltaT = 5.0; // deltaT = 5/N (N : NombreEspaceEntreLames = 1)
        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];        
        temperatureVect << -alphaE[0]*Es - he*Te - hg_1*_T0,
                    -alphaE[1]*Es - hi*Ti - hg_1*_T0;
    }
    else if(nombre_couches==3){
        deltaT = 2.5; // (N : NombreEspaceEntreLames = 2)
        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hg_2 = calculHg(_T0, deltaT, s_2, {1, 0, 0, 0}, 90)[0];
        temperatureVect << -alphaE[0]*Es - he*Te - hg_1*_T0,
                    -alphaE[1]*Es - hg_1*_T0 - hg_2*_T0,
                    -alphaE[2]*Es - hi*Ti - hg_2*_T0;
    }
    else if (nombre_couches==4)
    {
        deltaT = 2.5;
        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hg_2 = calculHg(_T0, deltaT, s_2, {1, 0, 0, 0}, 90)[0];
        hg_3 = calculHg(_T0, deltaT, s_3, {1, 0, 0, 0}, 90)[0];

        temperatureVect << -alphaE[0]*Es - he*Te - hg_1*_T0,
                    -alphaE[1]*Es - hg_1*_T0 - hg_2*_T0,
                    -alphaE[2]*Es - hg_2*_T0 - hg_3*_T0,
                    -alphaE[3]*Es - hg_3*_T0 - hi*Ti;
    }
    else if (nombre_couches==5)
    {
        deltaT = 2.5;
        hg_1 = calculHg(_T0, deltaT, s_1, {1, 0, 0, 0}, 90)[0];
        hg_2 = calculHg(_T0, deltaT, s_2, {1, 0, 0, 0}, 90)[0];
        hg_3 = calculHg(_T0, deltaT, s_3, {1, 0, 0, 0}, 90)[0];
        hg_4 = calculHg(_T0, deltaT, s_4, {1, 0, 0, 0}, 90)[0];

        temperatureVect << -alphaE[0]*Es - he*Te - hg_1*_T0,
                    -alphaE[1]*Es - hg_1*_T0 - hg_2*_T0,
                    -alphaE[2]*Es - hg_2*_T0 - hg_3*_T0,
                    -alphaE[3]*Es - hg_3*_T0 - hg_4*_T0,
                    -alphaE[4]*Es - hg_4*_T0 - hi*Ti;
    }
        
    return temperatureVect;
}

MatrixXd createQthMatrix(vector<double> ems, vector<double> ems_prime, vector<double> tau_th)
{   
    if(ems.size() != ems_prime.size() || ems.size() != tau_th.size()){
        throw std::invalid_argument("Size Error");
    }

    int nombre_couches = ems.size();
    int size = ems.size() * 2;
    MatrixXd thermalMatrix = MatrixXd::Zero(size, size);
    thermalMatrix.diagonal().setOnes();

    for(int i=0; i<nombre_couches; ++i){
        if(nombre_couches==2)
        {
            thermalMatrix(0, 3) = -(1 - ems_prime[0] - tau_th[0]);
            thermalMatrix(1, 3) = -tau_th[0];
            thermalMatrix(2, 0) = -tau_th[1];
            thermalMatrix(3, 0) = -(1 - ems[1] - tau_th[1]);
        }
        else if(nombre_couches==3)
        {
            thermalMatrix(0, 3) = -(1 - ems_prime[0] - tau_th[0]);
            thermalMatrix(1, 3) = -tau_th[0];
            thermalMatrix(2, 0) = -tau_th[1];
            thermalMatrix(3, 0) = -(1 - ems[1] - tau_th[1]);
            
            thermalMatrix(2, 5) = -(1 - ems_prime[1] - tau_th[1]);
            thermalMatrix(3, 5) = -tau_th[1];
            thermalMatrix(4, 2) = -tau_th[2];
            thermalMatrix(5, 2) = -(1 - ems[2] - tau_th[2]);
        }
        else if (nombre_couches==4)
        {
            thermalMatrix(0, 3) = -(1 - ems_prime[0] - tau_th[0]);
            thermalMatrix(1, 3) = -tau_th[0];
            thermalMatrix(2, 0) = -tau_th[1];
            thermalMatrix(3, 0) = -(1 - ems[1] - tau_th[1]);
            
            thermalMatrix(2, 5) = -(1 - ems_prime[1] - tau_th[1]);
            thermalMatrix(3, 5) = -tau_th[1];
            thermalMatrix(4, 2) = -tau_th[2];
            thermalMatrix(5, 2) = -(1 - ems[2] - tau_th[2]);

            thermalMatrix(4, 7) = -(1 - ems_prime[2] - tau_th[2]);
            thermalMatrix(5, 7) = -tau_th[2];
            thermalMatrix(6, 4) = -tau_th[3];
            thermalMatrix(7, 4) = -(1 - ems[3] - tau_th[3]);
        }
        else if (nombre_couches==5)
        {
            thermalMatrix(0, 3) = -(1 - ems_prime[0] - tau_th[0]);
            thermalMatrix(1, 3) = -tau_th[0];
            thermalMatrix(2, 0) = -tau_th[1];
            thermalMatrix(3, 0) = -(1 - ems[1] - tau_th[1]);
            
            thermalMatrix(2, 5) = -(1 - ems_prime[1] - tau_th[1]);
            thermalMatrix(3, 5) = -tau_th[1];
            thermalMatrix(4, 2) = -tau_th[2];
            thermalMatrix(5, 2) = -(1 - ems[2] - tau_th[2]);

            thermalMatrix(4, 7) = -(1 - ems_prime[2] - tau_th[2]);
            thermalMatrix(5, 7) = -tau_th[2];
            thermalMatrix(6, 4) = -tau_th[3];
            thermalMatrix(7, 4) = -(1 - ems[3] - tau_th[3]);

            thermalMatrix(6, 9) = -(1 - ems_prime[3] - tau_th[3]);
            thermalMatrix(7, 9) = -tau_th[3];
            thermalMatrix(8, 6) = -tau_th[4];
            thermalMatrix(9, 6) = -(1 - ems[4] - tau_th[4]);
        }     
    }        
    return thermalMatrix;
}

VectorXd createQthVector(vector<double> ems, vector<double> ems_prime, vector<double> tau_th, vector<double> q_th, vector<double> q_th_prime, vector<double> temperature)
{    
    int nombre_couches = temperature.size();
    VectorXd constantVector(temperature.size() * 2);

    for(int i=0; i<nombre_couches; ++i){
        if(nombre_couches==2)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * _StefanBoltzmann * pow(temperature[0], 4),
                            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * _StefanBoltzmann * pow(temperature[0], 4),
                            (1 - ems_prime[1] - tau_th[1]) * q_th_prime[2] + ems_prime[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            tau_th[1] * q_th_prime[2] + ems[1] * _StefanBoltzmann * pow(temperature[1], 4);
        }
    
        else if (nombre_couches==3)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * _StefanBoltzmann * pow(temperature[0], 4),
                            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * _StefanBoltzmann * pow(temperature[0], 4),
                            ems_prime[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            ems[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            (1 - ems_prime[2] - tau_th[2]) * q_th_prime[3] + ems_prime[2] * _StefanBoltzmann * pow(temperature[2], 4),
                            tau_th[2] * q_th_prime[3] + ems[2] * _StefanBoltzmann * pow(temperature[2], 4);
        }

        else if (nombre_couches==4)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * _StefanBoltzmann * pow(temperature[0], 4),
                            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * _StefanBoltzmann * pow(temperature[0], 4),
                            ems_prime[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            ems[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            ems_prime[2] * _StefanBoltzmann * pow(temperature[2], 4),
                            ems[2] * _StefanBoltzmann * pow(temperature[2], 4),
                            (1 - ems_prime[3] - tau_th[3]) * q_th_prime[4] + ems_prime[3] * _StefanBoltzmann * pow(temperature[3], 4),
                            tau_th[3] * q_th_prime[4] + ems[3] * _StefanBoltzmann * pow(temperature[3], 4);
        }

        else if (nombre_couches==5)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * _StefanBoltzmann * pow(temperature[0], 4),
                            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * _StefanBoltzmann * pow(temperature[0], 4),
                            ems_prime[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            ems[1] * _StefanBoltzmann * pow(temperature[1], 4),
                            ems_prime[2] * _StefanBoltzmann * pow(temperature[2], 4),
                            ems[2] * _StefanBoltzmann * pow(temperature[2], 4),
                            ems_prime[3] * _StefanBoltzmann * pow(temperature[3], 4),
                            ems[3] * _StefanBoltzmann * pow(temperature[3], 4),
                            (1 - ems_prime[4] - tau_th[4]) * q_th_prime[5] + ems_prime[4] * _StefanBoltzmann * pow(temperature[4], 4),
                            tau_th[4] * q_th_prime[5] + ems[4] * _StefanBoltzmann * pow(temperature[4], 4);
        }    
    }
    
    return constantVector;
}

VectorXd calculVitesse(vector<double> T, double Ti, double hauteur, double largeur, double s) //vector<double>& tau_th)
{   
    /* Cette fonction est utilisé pour calculer hc et Tg en fonction de la vitesse */
    int nombre_couches = T.size();

    double Z1 = 1.041;
    double Z2 = 8.536;
        
    double T1 = T[0];
    double T2 = T[1];
    double T3 = (nombre_couches >= 3) ? T[2] : 0;
    double T4 = (nombre_couches >= 4) ? T[3] : 0;
    double T5 = (nombre_couches == 5) ? T[4] : 0;

    double Tgi;
    double rho, muy, lambda, c, hg;

    double Tmp_1_2, Tmp_2_3, Tmp_3_4, Tmp_4_5;
    double deltaT_1_2, deltaT_2_3, deltaT_3_4, deltaT_4_5;

    //double Te = 298.15;

    double vitesse = 1.0;  // Initialisation de la vitesse
    double tol = 0.0000001;     // Tolérance de convergence pour la vitesse
    double vitesse_precedente;
    bool continueLoop = true;
    int max_iter = 100;     // Pour éviter la boucle infinie
    int iter = 0;

    double DPB, DPHP, DPZ, DPT;
    double A, B, C;
    double Htp;
    double Tg;
    double Tsor;
    double hc;
    double Tgjv;

    if(nombre_couches == 2){
        Tmp_1_2 = (T1 + T2) / 2.0;
        Tgi = (Tmp_1_2 + Ti) / 2.0;
        deltaT_1_2 = fabs(T2 - T1);
        rho = calculHg(_T0, deltaT_1_2, s, {1, 0, 0, 0}, 90)[1];
        muy = calculHg(_T0, deltaT_1_2, s, {1, 0, 0, 0}, 90)[2];
        c = calculHg(_T0, deltaT_1_2, s, {1, 0, 0, 0}, 90)[3];
        hg = calculHg(_T0, deltaT_1_2, s, {1, 0, 0, 0}, 90)[0];
    }

    else if (nombre_couches == 3)
    {
        Tmp_2_3 = (T2 + T3) / 2.0;
        Tgi = (Tmp_2_3 + Ti) / 2.0;
        deltaT_2_3 = fabs(T3 - T2);
        rho = calculHg(_T0, deltaT_2_3, s, {1, 0, 0, 0}, 90)[1];
        muy = calculHg(_T0, deltaT_2_3, s, {1, 0, 0, 0}, 90)[2];
        c = calculHg(_T0, deltaT_2_3, s, {1, 0, 0, 0}, 90)[3];
        hg = calculHg(_T0, deltaT_2_3, s, {1, 0, 0, 0}, 90)[0];
    }

    else if (nombre_couches == 4)
    {
        Tmp_3_4 = (T3 + T4) / 2.0;
        Tgi = (Tmp_3_4 + Ti) / 2.0;
        deltaT_3_4 = fabs(T4 - T3);
        rho = calculHg(_T0, deltaT_3_4, s, {1, 0, 0, 0}, 90)[1];
        muy = calculHg(_T0, deltaT_3_4, s, {1, 0, 0, 0}, 90)[2];
        c = calculHg(_T0, deltaT_3_4, s, {1, 0, 0, 0}, 90)[3];
        hg = calculHg(_T0, deltaT_3_4, s, {1, 0, 0, 0}, 90)[0];
    }

    else if (nombre_couches == 5)
    {
        Tmp_4_5 = (T4 + T5) / 2.0;
        Tgi = (Tmp_4_5 + Ti) / 2.0;
        deltaT_4_5 = fabs(T5 - T4);
        rho = calculHg(_T0, deltaT_4_5, s, {1, 0, 0, 0}, 90)[1];
        muy = calculHg(_T0, deltaT_4_5, s, {1, 0, 0, 0}, 90)[2];
        c = calculHg(_T0, deltaT_4_5, s, {1, 0, 0, 0}, 90)[3];
        hg = calculHg(_T0, deltaT_4_5, s, {1, 0, 0, 0}, 90)[0];
    }
    
    DPB = (rho * pow(vitesse, 2)) / 2;
    DPHP = (12 * muy * hauteur * vitesse) / pow(s, 2);
    DPZ = rho * vitesse * ((Z1 + Z2) / 2);
    
    while (continueLoop && iter < max_iter) {
        vitesse_precedente = vitesse;         
        DPT = 1.189 * 293 * 9.81 * hauteur * fabs(Tgi - Ti) / (Tgi * Ti);
        
        /* A = DPB + DPZ;
        B = DPHP;
        C = -DPT; */
        A = 0.5*rho*(1+Z1+Z2);
        B = 12.0 * muy * hauteur / pow(s,2);
        C = -DPT;

        vitesse = (-B+pow(B*B - 4*A*C,0.5))/2/A; 

        if(nombre_couches == 2){
            hc = 2.0 * hg + 4.0 * vitesse;
            Htp = rho * c * s * vitesse / (2.0 * hc);
            Tsor = Tmp_1_2 - (Tmp_1_2 - Ti) * exp(-hauteur / Htp);
            Tg = Tmp_1_2 - (Htp / hauteur) * (Tsor - Ti);
            hg = calculHg(Tg, deltaT_1_2, s, {1, 0, 0, 0}, 90)[0]; 
        }
        
        else if(nombre_couches == 3)
        {
            hc = 2 * hg + 4 * vitesse;
            Htp = rho * c * s * vitesse / (2 * hc);
            Tsor = Tmp_2_3 - (Tmp_2_3 - Ti) * exp(-hauteur / Htp);
            Tg = Tmp_2_3 - (Htp / hauteur) * (Tsor - Ti);
            hg = calculHg(Tg, deltaT_2_3, s, {1, 0, 0, 0}, 90)[0];
        }

        else if (nombre_couches == 4)
        {
            hc = 2 * hg + 4 * vitesse;
            Htp = rho * c * s * vitesse / (2 * hc);
            Tsor = Tmp_3_4 - (Tmp_3_4 - Ti) * exp(-hauteur / Htp);
            Tg = Tmp_3_4 - (Htp / hauteur) * (Tsor - Ti);
            hg = calculHg(Tg, deltaT_3_4, s, {1, 0, 0, 0}, 90)[0];
        }

        else if (nombre_couches == 5)
        {
            hc = 2 * hg + 4 * vitesse;
            Htp = rho * c * s * vitesse / (2 * hc);
            Tsor = Tmp_4_5 - (Tmp_4_5 - Ti) * exp(-hauteur / Htp);
            Tg = Tmp_4_5 - (Htp / hauteur) * (Tsor - Ti);
            hg = calculHg(Tg, deltaT_4_5, s, {1, 0, 0, 0}, 90)[0];
        }
        
        // Mise à jour de Tgi pour le prochain calcul de DPT
        Tgi = Tg;
        // Vérifier la condition de sortie de la boucle
        if (fabs(vitesse - vitesse_precedente) < tol) {
            continueLoop = false;
        }

        // Mise à jour de Tgi pour le prochain calcul de DPT
        Tgi = Tg;
        // Vérifier la condition de sortie de la boucle
        if (fabs(vitesse - vitesse_precedente) < tol) {
            continueLoop = false;
        }
        iter++;
    }

    VectorXd result(2);
    result[0] = hc;
    result[1] = Tg;

    return result;
}

VectorXd bilan(double Es, vector<double> temperatureVector,vector<double> ems, vector<double> ems_prime, vector<double> alphaE, vector<double> tau_th, 
            double hc_e, double hc_i, double Te, double Ti, double hauteur, double largeur, vector<double> s, int espace)
{   
    int nombre_couches = temperatureVector.size();

    double s1 = s[0];
    double s2 = s[1];
    double s3 = (nombre_couches >= 4) ? s[2] : 0;
    double s4 = (nombre_couches == 5) ? s[3] : 0;

    double T0 = 293.15;
    double T1 = temperatureVector[0];
    double T2 = temperatureVector[1];
    double T3 = (nombre_couches >= 3) ? temperatureVector[2] : 0;
    double T4 = (nombre_couches >= 4) ? temperatureVector[3] : 0;
    double T5 = (nombre_couches == 5) ? temperatureVector[4] : 0;

    double qth_0, qth_1, qth_2, qth_3, qth_4, qth_5;
    double qth_0_prime, qth_1_prime, qth_2_prime, qth_3_prime, qth_4_prime, qth_5_prime;

    double bilan_1, bilan_2, bilan_3, bilan_4, bilan_5;
    double df1_dT1, df1_dT2, df1_dT3, df1_dT4, df1_dT5;
    double df2_dT1, df2_dT2, df2_dT3, df2_dT4, df2_dT5;
    double df3_dT1, df3_dT2, df3_dT3, df3_dT4, df3_dT5;
    double df4_dT1, df4_dT2, df4_dT3, df4_dT4, df4_dT5;
    double df5_dT1, df5_dT2, df5_dT3, df5_dT4, df5_dT5;

    double Tmp_1_2, Tmp_2_3, Tmp_3_4, Tmp_4_5;
    double deltaT_1_2, deltaT_2_3, deltaT_3_4, deltaT_4_5;

    double hg_0 = hc_e;
    double hg_1;
    double hg_2;
    double hg_3;
    double hg_4;
   
    double hc;
    double hg;
    double Tg;
    double Tsor;
    double Htp;
    double vitesse;
    double Tgjv;
    
    VectorXd v; //vecteur pour stocker les résultat de la fonction calculVitesse()
    MatrixXd qthConstantMatrix;
    MatrixXd qthConstantVector;
    VectorXd qthVector;

    double qth_a1, qth_a2, qth_a3, qth_a4, qth_a5;
    double qc_a1, qc_a2, qc_a3, qc_a4, qc_a5;
    
    // Newton-Raphson iteration loop
    for (int iter = 0; iter < 100; ++iter) {
        if(nombre_couches == 2){
            v = calculVitesse({T1, T2}, Ti, hauteur, largeur, s2);
            hc = v[0];
            Tg = v[1];

            Tmp_1_2 = (T1 + T2) / 2.0;
            deltaT_1_2 = fabs(T2 - T1);

            hg_1 = calculHg(Tmp_1_2, deltaT_1_2, s2, {1, 0, 0, 0}, 90)[0];
            hg_2 = hc_i;
        
            // Recalculate thermophysical properties
            qthConstantMatrix = createQthMatrix({ems[0], ems[1]}, {ems_prime[0], ems_prime[1]}, {tau_th[0], tau_th[1]});
            qthConstantVector = createQthVector({ems[0], ems[1]}, {ems_prime[0], ems_prime[1]}, {tau_th[0], tau_th[1]},
                                                        {qth_0, qth_1, qth_2}, {qth_0_prime, qth_1_prime, qth_2_prime}, {T1, T2});
            qthVector = qthConstantMatrix.inverse() * qthConstantVector;

            qth_0 = _StefanBoltzmann * pow(Te,4);
            qth_1 = qthVector[0];
            qth_0_prime = qthVector[1]; 
            qth_2 = qthVector[2];
            qth_2_prime = _StefanBoltzmann * pow(Ti, 4); 
            qth_1_prime = qthVector[3];

            qth_a1 = ems[0] * qth_0 + ems_prime[0] * qth_1_prime  - (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 4);
            qth_a2 = ems[1] * qth_1 + ems_prime[1] * qth_2_prime  - (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 4);
            
            if(espace == 0){
                qc_a1 = hg_0 * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (Ti - T2);
                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
            }

            else if (espace == 1){
                qc_a1 = hc_e * (Te - T1) + hc * (Tg - T1);
                qc_a2 = hc * (Tg - T2) + hc_i * (Ti - T2);
                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
            }

            VectorXd equation(2);
            equation[0] = bilan_1;
            equation[1] = bilan_2;

            // Update Jacobian matrix
            df1_dT1 = (espace==1) ? -4.0 * (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hc_e - hc : -4.0 * (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hg_0 - hg_1;
            df1_dT2 = (espace==1) ? 0 : hg_1;
            df2_dT1 = (espace==1) ? 0 : hg_1;
            df2_dT2 = (espace==1) ? -4.0 * (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hc - hc_i : -4.0 * (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hg_1 - hg_2;

            MatrixXd jacobian(2, 2);
            jacobian(0, 0) = df1_dT1;
            jacobian(0, 1) = df1_dT2;
            jacobian(1, 0) = df2_dT1;
            jacobian(1, 1) = df2_dT2;

            VectorXd delta = jacobian.colPivHouseholderQr().solve(-equation);  // small adjustment to see effect 

            // Check convergence
            if (fabs(bilan_1) < 0.1 && fabs(bilan_2) < 0.1) {
                break;
            } 
            // Update temperatures
            T1 += delta[0];
            T2 += delta[1];
        }

        else if (nombre_couches == 3)
        {
            v = calculVitesse({T1, T2, T3}, Ti, hauteur, largeur, s2);            
            hc = v[0];
            Tg = v[1];

            Tmp_1_2 = (T1 + T2) / 2;
            Tmp_2_3 = (T2 + T3) / 2;
            deltaT_1_2 = fabs(T2 - T1);
            deltaT_2_3 = fabs(T3 - T2);

            hg_0 = hc_e;
            hg_1 = calculHg(Tmp_1_2, deltaT_1_2, s1, {0.1, 0.9, 0, 0}, 90)[0];
            hg_2 = calculHg(Tg, deltaT_2_3, s2, {1, 0, 0, 0}, 90)[0]; 
            hg_3 = hc_i;

            qthConstantMatrix = createQthMatrix({ems[0], ems[1], ems[2]}, {ems_prime[0], ems_prime[1], ems_prime[2]}, {tau_th[0], tau_th[1], tau_th[2]});
            qthConstantVector = createQthVector({ems[0], ems[1], ems[2]}, {ems_prime[0], ems_prime[1], ems_prime[2]}, {tau_th[0], tau_th[1], tau_th[2]},
                                                {qth_0, qth_1, qth_2, qth_3}, {qth_0_prime, qth_1_prime, qth_2_prime, qth_3_prime}, {T1, T2, T3});
            qthVector = qthConstantMatrix.inverse() * qthConstantVector;

            qth_0 = _StefanBoltzmann * pow(Te,4);
            qth_1 = qthVector[0];
            qth_0_prime = qthVector[1];
            qth_2 = qthVector[2];
            qth_1_prime = qthVector[3];
            qth_3 = qthVector[4];
            qth_2_prime = qthVector[5];
            qth_3_prime = _StefanBoltzmann * pow(Ti, 4);

            qth_a1 = ems[0] * qth_0 + ems_prime[0] * qth_1_prime  - (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 4);
            qth_a2 = ems[1] * qth_1 + ems_prime[1] * qth_2_prime  - (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 4);
            qth_a3 = ems[2] * qth_2 + ems_prime[2] * qth_3_prime  - (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 4);

            if(espace == 0){ // Espace fermé
                hg_2 = calculHg(Tmp_2_3, deltaT_2_3, s2, {0.1, 0.9, 0, 0}, 90)[0]; 

                qc_a1 = hc_e * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hc_i * (Ti - T3);

                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2 ;
                bilan_3 = alphaE[2] * Es + qth_a3 + qc_a3;
            }
            
            else if(espace == 1){ // Espace ouvert 
                qc_a1 = hc_e * (Te - T1) + hg_1 * (Tmp_1_2 - T1);
                qc_a2 = hg_1 * (Tmp_1_2 - T2) + hc * (Tg - T2);
                qc_a3 = hc * (Tg - T3) + hc_i * (Ti - T3);

                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
                bilan_3 = alphaE[2] * Es + qth_a3 + qc_a3;
            }

            VectorXd equation(3);
            equation[0] = bilan_1;
            equation[1] = bilan_2;
            equation[2] = bilan_3;

            // Update Jacobian matrix
            df1_dT1 = (espace==1) ? -4.0 * (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hc_e - hg_1 : -4.0 * (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hg_0 - hg_1;
            df1_dT2 = (espace==1) ? 0 : hg_1;
            df1_dT3 = 0;
            
            df2_dT1 = (espace==1) ? 0 : hg_1;
            df2_dT2 = (espace==1) ? -4.0 * (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hg_1 - hc : -4.0 * (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hg_1 - hg_2;
            df2_dT3 = 0;
            
            df3_dT1 = 0;
            df3_dT2 = (espace==1) ? 0 : hg_2;
            df3_dT3 = (espace==1) ? -4.0 * (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 3) - hc - hc_i : -4.0 * (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 3) - hg_2 - hg_3;

            MatrixXd jacobian(3, 3);
            jacobian(0, 0) = df1_dT1;
            jacobian(0, 1) = df1_dT2;
            jacobian(0, 2) = df1_dT3;
            jacobian(1, 0) = df2_dT1;
            jacobian(1, 1) = df2_dT2;
            jacobian(1, 2) = df2_dT3;
            jacobian(2, 0) = df3_dT1;
            jacobian(2, 1) = df3_dT2;
            jacobian(2, 2) = df3_dT3;

            VectorXd delta = jacobian.colPivHouseholderQr().solve(-equation);  // small adjustment to see effect 

            // Check convergence
            if (fabs(bilan_1) < 0.1 && fabs(bilan_2) < 0.1 && fabs(bilan_3) < 0.1) {
                break;
            } 
            
            // Update temperatures
            T1 += delta[0];
            T2 += delta[1];
            T3 += delta[2];
        }

        else if (nombre_couches==4)
        {
            v = calculVitesse({T2, T2, T3, T4}, Ti, hauteur, largeur, s3);
            hc = v[0];
            Tg = v[1];
            
            Tgjv = (Tmp_1_2 + Te) / 2;

            Tmp_1_2 = (T1 + T2) / 2;
            Tmp_2_3 = (T2 + T3) / 2;
            Tmp_3_4 = (T3 + T4) / 2;
            deltaT_1_2 = fabs(T2 - T1);
            deltaT_2_3 = fabs(T3 - T2);
            deltaT_3_4 = fabs(T4 - T3);

            hg_1 = calculHg(Tmp_1_2, deltaT_1_2, s1, {0.1, 0.9, 0, 0}, 90)[0];
            hg_2 = calculHg(Tmp_2_3, deltaT_2_3, s2, {0.1, 0.9, 0, 0}, 90)[0];
            hg_3 = calculHg(Tg, deltaT_3_4, s3, {1, 0, 0, 0}, 90)[0];

            qthConstantMatrix = createQthMatrix({ems[0], ems[1], ems[2], ems[3]}, {ems_prime[0], ems_prime[1], ems_prime[2], ems_prime[3]}, {tau_th[0], tau_th[1], tau_th[2], tau_th[3]});
            qthConstantVector = createQthVector({ems[0], ems[1], ems[2], ems[3]}, {ems_prime[0], ems_prime[1], ems_prime[2], ems_prime[3]}, {tau_th[0], tau_th[1], tau_th[2], tau_th[3]},
                                                {qth_0, qth_1, qth_2, qth_3, qth_4}, {qth_0_prime, qth_1_prime, qth_2_prime, qth_3_prime, qth_4_prime}, {T1, T2, T3, T4});
            qthVector = qthConstantMatrix.inverse() * qthConstantVector;

            qth_0 = _StefanBoltzmann * pow(Te, 4);
            qth_1 = qthVector[0];
            qth_0_prime = qthVector[1];
            qth_2 = qthVector[2];
            qth_1_prime = qthVector[3];
            qth_3 = qthVector[4];
            qth_2_prime = qthVector[5];
            qth_4 = qthVector[6];
            qth_3_prime = qthVector[7];
            qth_4_prime = _StefanBoltzmann * pow(Ti, 4);

            qth_a1 = ems[0] * qth_0 + ems_prime[0] * qth_1_prime - (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 4);
            qth_a2 = ems[1] * qth_1 + ems_prime[1] * qth_2_prime - (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 4);
            qth_a3 = ems[2] * qth_2 + ems_prime[2] * qth_3_prime - (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 4);
            qth_a4 = ems[3] * qth_3 + ems_prime[3] * qth_4_prime - (ems[3] + ems_prime[3]) * _StefanBoltzmann * pow(T4, 4);

            if(espace==1){
                qc_a1 = hc_e * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hc   * (Tg - T3);
                qc_a4 = hc   * (Tg - T4) + hc_i * (Ti - T4);

                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
                bilan_3 = alphaE[2] * Es + qth_a3 + qc_a3;
                bilan_4 = alphaE[3] * Es + qth_a4 + qc_a4;
            }

            else if (espace==0){
                hg_3 = calculHg(Tmp_3_4, deltaT_3_4, s3, {0.1, 0.9, 0, 0}, 90)[0];

                qc_a1 = hc_e * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hg_3 * (T4 - T3);
                qc_a4 = hg_3 * (T3 - T4) + hc_i * (Ti - T4);

                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
                bilan_3 = alphaE[2] * Es + qth_a3 + qc_a3;
                bilan_4 = alphaE[3] * Es + qth_a4 + qc_a4;
            }

            VectorXd equation(4);
            equation[0] = bilan_1;
            equation[1] = bilan_2;
            equation[2] = bilan_3;
            equation[3] = bilan_4;

            df1_dT1 = - 4 * (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hc_e - hg_1;
            df1_dT2 = hg_1;
            df1_dT3 = 0;
            df1_dT4 = 0;

            df2_dT1 = hg_1;
            df2_dT2 = -4 * (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hg_1 - hg_2;;
            df2_dT3 = hg_2;
            df2_dT4 = 0;

            df3_dT1 = 0;
            df3_dT2 = hg_2;
            df3_dT3 = (espace == 1) ? -4 * (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 3) - hg_2 - hc : -4 * (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 3) - hg_2 - hg_3;
            df3_dT4 = (espace == 1) ? 0 : hg_3;

            df4_dT1 = 0;
            df4_dT2 = 0;
            df4_dT3 = (espace == 1) ? 0 : hg_3;
            df4_dT4 = (espace == 1) ? -4 * (ems[3] + ems_prime[3]) * _StefanBoltzmann * pow(T4, 3) - hc - hc_i : -4 * (ems[3] + ems_prime[3]) * _StefanBoltzmann * pow(T4, 3) - hg_3 - hc_i;

            MatrixXd jacobian(4, 4);
            jacobian(0, 0) = df1_dT1;
            jacobian(0, 1) = df1_dT2;
            jacobian(0, 2) = df1_dT3;
            jacobian(0, 3) = df1_dT4;

            jacobian(1, 0) = df2_dT1;
            jacobian(1, 1) = df2_dT2;
            jacobian(1, 2) = df2_dT3;
            jacobian(1, 3) = df2_dT4;

            jacobian(2, 0) = df3_dT1;
            jacobian(2, 1) = df3_dT2;
            jacobian(2, 2) = df3_dT3;
            jacobian(2, 3) = df3_dT4;
            
            jacobian(3, 0) = df4_dT1;
            jacobian(3, 1) = df4_dT2;
            jacobian(3, 2) = df4_dT3;
            jacobian(3, 3) = df4_dT4;

            VectorXd delta = jacobian.colPivHouseholderQr().solve(-equation);  // small adjustment to see effect 

            //cout << "Iteration " << iter << ":  T1 = " << T1 << ", T2 = " << T2 << ", T3 = " << T3 << ", T4 = " << T4 <<endl;
            //cout << "Iteration " << iter << ": Bilan 1 = " << bilan_1 << ", Bilan 2 = " << bilan_2 << ", Bilan 3 = " << bilan_3 << ", Bilan 4 = " << bilan_4 << endl; 

            // Check convergence
            if (fabs(bilan_1) < 0.1 && fabs(bilan_2) < 0.1 && fabs(bilan_3) < 0.1 && fabs(bilan_4) < 0.1) {
                break;
            } 
            
            // Update temperatures
            T1 += delta[0];
            T2 += delta[1];
            T3 += delta[2];
            T4 += delta[3];
        }

        else if (nombre_couches==5)
        {
            v = calculVitesse({T2, T2, T3, T4, T5}, Ti, hauteur, largeur, s4);
            hc = v[0];
            Tg = v[1];
            
            //Tgjv = (Tmp_1_2 + Te) / 2;

            Tmp_1_2 = (T1 + T2) / 2;
            Tmp_2_3 = (T2 + T3) / 2;
            Tmp_3_4 = (T3 + T4) / 2;
            Tmp_4_5 = (T4 + T5) / 2;
            deltaT_1_2 = fabs(T2 - T1);
            deltaT_2_3 = fabs(T3 - T2);
            deltaT_3_4 = fabs(T4 - T3);
            deltaT_4_5 = fabs(T5 - T4); 

            hg_1 = calculHg(Tmp_1_2, deltaT_1_2, s1, {0.1, 0.9, 0, 0}, 90)[0];
            hg_2 = calculHg(Tmp_2_3, deltaT_2_3, s2, {0.1, 0.9, 0, 0}, 90)[0];
            hg_3 = calculHg(Tmp_3_4, deltaT_3_4, s3, {0.1, 0.9, 0, 0}, 90)[0];
            hg_4 = calculHg(Tg, deltaT_4_5, s4, {1, 0, 0, 0}, 90)[0];

            qthConstantMatrix = createQthMatrix({ems[0], ems[1], ems[2], ems[3], ems[4]}, {ems_prime[0], ems_prime[1], ems_prime[2], ems_prime[3], ems_prime[4]}, 
                                                {tau_th[0], tau_th[1], tau_th[2], tau_th[3], tau_th[4]});
            qthConstantVector = createQthVector({ems[0], ems[1], ems[2], ems[3], ems[4]}, {ems_prime[0], ems_prime[1], ems_prime[2], ems_prime[3], ems_prime[4]}, 
                                                {tau_th[0], tau_th[1], tau_th[2], tau_th[3], tau_th[4]}, {qth_0, qth_1, qth_2, qth_3, qth_4, qth_5}, 
                                                {qth_0_prime, qth_1_prime, qth_2_prime, qth_3_prime, qth_4_prime, qth_5_prime}, 
                                                {T1, T2, T3, T4, T5});
            qthVector = qthConstantMatrix.inverse() * qthConstantVector;

            qth_0 = _StefanBoltzmann * pow(Te, 4);
            qth_1 = qthVector[0];
            qth_0_prime = qthVector[1];
            qth_2 = qthVector[2];
            qth_1_prime = qthVector[3];
            qth_3 = qthVector[4];
            qth_2_prime = qthVector[5];
            qth_4 = qthVector[6];
            qth_3_prime = qthVector[7];
            qth_5 = qthVector[8];
            qth_4_prime = qthVector[9];            
            qth_5_prime = _StefanBoltzmann * pow(Ti, 4);

            qth_a1 = ems[0] * qth_0 + ems_prime[0] * qth_1_prime - (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 4);
            qth_a2 = ems[1] * qth_1 + ems_prime[1] * qth_2_prime - (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 4);
            qth_a3 = ems[2] * qth_2 + ems_prime[2] * qth_3_prime - (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 4);
            qth_a4 = ems[3] * qth_3 + ems_prime[3] * qth_4_prime - (ems[3] + ems_prime[3]) * _StefanBoltzmann * pow(T4, 4);
            qth_a5 = ems[4] * qth_4 + ems_prime[4] * qth_5_prime - (ems[4] + ems_prime[4]) * _StefanBoltzmann * pow(T5, 4);

            if(espace==1){
                qc_a1 = hc_e * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hg_3 * (T4 - T3);
                qc_a4 = hg_3 * (T3 - T4) + hc   * (T5 - T4);
                qc_a5 = hc   * (Tg - T5) + hc_i * (Ti - T5);

                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
                bilan_3 = alphaE[2] * Es + qth_a3 + qc_a3;
                bilan_4 = alphaE[3] * Es + qth_a4 + qc_a4;
                bilan_5 = alphaE[4] * Es + qth_a5 + qc_a5;
            }

            else if (espace==0){
                hg_4 = calculHg(Tmp_4_5, deltaT_4_5, s4, {0.1, 0.9, 0, 0}, 90)[0];

                qc_a1 = hc_e * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hg_3 * (T4 - T3);
                qc_a4 = hg_3 * (T3 - T4) + hg_4 * (T5 - T4);
                qc_a5 = hg_4 * (T4 - T5) + hc_i * (Ti - T5);

                bilan_1 = alphaE[0] * Es + qth_a1 + qc_a1;
                bilan_2 = alphaE[1] * Es + qth_a2 + qc_a2;
                bilan_3 = alphaE[2] * Es + qth_a3 + qc_a3;
                bilan_4 = alphaE[3] * Es + qth_a4 + qc_a4;
                bilan_5 = alphaE[4] * Es + qth_a5 + qc_a5;
            }

            VectorXd equation(5);
            equation[0] = bilan_1;
            equation[1] = bilan_2;
            equation[2] = bilan_3;
            equation[3] = bilan_4;
            equation[4] = bilan_5;

            df1_dT1 = - 4 * (ems[0] + ems_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hc_e - hg_1;
            df1_dT2 = hg_1;
            df1_dT3 = 0;
            df1_dT4 = 0;
            df1_dT5 = 0;

            df2_dT1 = hg_1;
            df2_dT2 = -4 * (ems[1] + ems_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hg_1 - hg_2;;
            df2_dT3 = hg_2;
            df2_dT4 = 0;
            df2_dT5 = 0;
            
            df3_dT1 = 0;
            df3_dT2 = hg_2;
            df3_dT3 = -4 * (ems[2] + ems_prime[2]) * _StefanBoltzmann * pow(T3, 3) - hg_2 - hg_3;
            df3_dT4 = hg_3;
            df3_dT5 = 0;
            
            df4_dT1 = 0;
            df4_dT2 = 0;
            df4_dT3 = hg_3;
            df4_dT4 = (espace == 1) ? -4 * (ems[3] + ems_prime[3]) * _StefanBoltzmann * pow(T4, 3) - hg_3 - hc : -4 * (ems[3] + ems_prime[3]) * _StefanBoltzmann * pow(T4, 3) - hg_3 - hg_4;
            df4_dT5 = (espace == 1) ? hc : hg_4;

            df5_dT1 = 0;
            df5_dT2 = 0;
            df5_dT3 = 0;
            df5_dT4 = (espace == 1) ? 0 : hg_4;
            df5_dT5 = (espace == 1) ? -4 * (ems[4] + ems_prime[4]) * _StefanBoltzmann * pow(T5, 3) - hc - hc_i : -4 * (ems[4] + ems_prime[4]) * _StefanBoltzmann * pow(T5, 3) - hg_4 - hc_i;

            MatrixXd jacobian(5, 5);
            jacobian(0, 0) = df1_dT1;
            jacobian(0, 1) = df1_dT2;
            jacobian(0, 2) = df1_dT3;
            jacobian(0, 3) = df1_dT4;
            jacobian(0, 4) = df1_dT5;

            jacobian(1, 0) = df2_dT1;
            jacobian(1, 1) = df2_dT2;
            jacobian(1, 2) = df2_dT3;
            jacobian(1, 3) = df2_dT4;
            jacobian(1, 4) = df2_dT5;

            jacobian(2, 0) = df3_dT1;
            jacobian(2, 1) = df3_dT2;
            jacobian(2, 2) = df3_dT3;
            jacobian(2, 3) = df3_dT4;
            jacobian(2, 4) = df3_dT5;
            
            jacobian(3, 0) = df4_dT1;
            jacobian(3, 1) = df4_dT2;
            jacobian(3, 2) = df4_dT3;
            jacobian(3, 3) = df4_dT4;
            jacobian(3, 4) = df4_dT5;

            jacobian(4, 0) = df5_dT1;
            jacobian(4, 1) = df5_dT2;
            jacobian(4, 2) = df5_dT3;
            jacobian(4, 3) = df5_dT4;
            jacobian(4, 4) = df5_dT5;

            VectorXd delta = jacobian.colPivHouseholderQr().solve(-equation);  // small adjustment to see effect 

            cout << "Iteration " << iter << ":  T1 = " << T1 << ", T2 = " << T2 << ", T3 = " << T3 << ", T4 = " << T4 << ", T5 = " << T5 << endl;
            cout << "Iteration " << iter << ": Bilan 1 = " << bilan_1 << ", Bilan 2 = " << bilan_2 << ", Bilan 3 = " << bilan_3 << ", Bilan 4 = " << bilan_4 << ", Bilan 5 = " << bilan_5 << endl; 

            // Check convergence
            if (fabs(bilan_1) < 0.1 && fabs(bilan_2) < 0.1 && fabs(bilan_3) < 0.1 && fabs(bilan_4) < 0.1 && fabs(bilan_5) < 0.1){
                break;
            } 
            
            // Update temperatures
            T1 += delta[0];
            T2 += delta[1];
            T3 += delta[2];
            T4 += delta[3];
            T5 += delta[4];
        }            
    }

    cout << "\nTg = " << Tg << endl;
    //cout << "Tgjv = " << Tgjv << endl;
    cout << "hc = " << hc << endl;
    cout << "hg_1 = " << hg_1 << endl;
    cout << "hg_2 = " << hg_2 << endl;
    cout << "hg_3 = " << hg_3 << endl; 
    cout << "hg_4 = " << hg_4 << endl; 
    
    VectorXd result = VectorXd(11);
    result[0] = T1;
    result[1] = T2;
    result[2] = (nombre_couches >= 3) ? T3 : 0;
    result[3] = (nombre_couches >= 4) ? T4 : 0;
    result[4] = (nombre_couches == 5) ? T5 : 0;
    result[5] = Tg;
    result[6] = hc;
    result[7] = qth_2;
    result[8] = (nombre_couches >= 3) ? qth_3 : 0;
    result[9] = (nombre_couches >= 4) ? qth_4 : 0;
    //result[9] = (nombre_couches >= 4) ? Tgjv : 0;
    result[10] = (nombre_couches ==5) ? qth_5 : 0;

    return result;
}

vector<double> calculFacteurSolaire(double Es, double transmission, double Te, double Ti, double hc_e, double hc_i, vector<double> ems, vector<double> ems_prime, 
                        vector<double> alphaE, vector<double> tau_th, double hauteur, double largeur, vector<double> s, int type_espace, int position_store)
{   
    int nombre_couches = ems.size();
    vector<double> temperatureInitialAvecFlux, temperatureInitialSansFlux;
    double g_th, g_c, g_v, g_tot;

    double qth_0 = _StefanBoltzmann * pow(Te,4);
    double qth_n_prime = _StefanBoltzmann * pow(Ti,4);

    /*------------------------------------------------------- resultat avec flux -------------------------------------------------------*/
    MatrixXd temperatureMatrixAvecFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, ems, ems_prime, s);
    VectorXd constantVectorTemperatureAvecFlux = createTemperatureVector(Es, hc_e, hc_i, Te, Ti, alphaE, s);
    VectorXd temperatureVectorAvecFlux = temperatureMatrixAvecFlux.inverse() * constantVectorTemperatureAvecFlux;

    for (int i = 0; i < nombre_couches; ++i) {
        temperatureInitialAvecFlux.push_back(temperatureVectorAvecFlux[i]);
    }

    VectorXd bilanResultAvecFlux = bilan(Es, temperatureInitialAvecFlux, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, s, type_espace);    
    
    double T1_avecFlux = bilanResultAvecFlux[0];
    double T2_avecFlux = bilanResultAvecFlux[1];
    double T3_avecFlux = (nombre_couches >= 3) ? bilanResultAvecFlux[2] : 0;
    double T4_avecFlux = (nombre_couches >= 4) ? bilanResultAvecFlux[3] : 0;
    double T5_avecFlux = (nombre_couches == 5) ? bilanResultAvecFlux[4] : 0;
    double Tg_avecFlux = bilanResultAvecFlux[5];
    double hc_avecFlux = bilanResultAvecFlux[6];
    double qth_2_avecFlux = bilanResultAvecFlux[7];
    double qth_3_avecFlux = bilanResultAvecFlux[8];
    double qth_4_avecFlux = bilanResultAvecFlux[9];
    double qth_5_avecFlux = bilanResultAvecFlux[10];

    cout << "\nT1(Es) = " << T1_avecFlux - 273.15 << endl;
    cout << "T2(Es) = " << T2_avecFlux - 273.15 << endl;
    cout << "T3(Es) = " << T3_avecFlux - 273.15 << endl;
    cout << "T4(Es) = " << T4_avecFlux - 273.15 << endl;
    cout << "T5(Es) = " << T5_avecFlux - 273.15 << endl;

    double qth_i_avecFlux = (nombre_couches == 2) ? (qth_2_avecFlux - qth_n_prime) : 
                            (nombre_couches == 3) ? (qth_3_avecFlux - qth_n_prime) :
                            (nombre_couches == 4) ? (qth_4_avecFlux - qth_n_prime) : 
                            (nombre_couches == 5) ? (qth_5_avecFlux - qth_n_prime) : 0;
    
    cout << "\nqth_i(Es) = " << qth_i_avecFlux << endl;

    double qc_i_avecFlux = (nombre_couches == 2) ? (hc_i * (Ti - T2_avecFlux)) : 
                           (nombre_couches == 3) ? (hc_i * (Ti - T3_avecFlux)) :
                           (nombre_couches == 4) ? (hc_i * (Ti - T4_avecFlux)) :
                           (nombre_couches == 5) ? (hc_i * (Ti - T5_avecFlux)) : 0;
    
    cout << "qc_i(Es) = " << qc_i_avecFlux << endl;

    //double Tgjv = bilanResultAvecFlux[9];
    
    double q_v_avecFlux = (nombre_couches == 2) ? (hc_avecFlux * (Tg_avecFlux - T1_avecFlux) + hc_avecFlux * (Tg_avecFlux - T2_avecFlux)) : 
                          (nombre_couches == 3) ? (hc_avecFlux * (Tg_avecFlux - T2_avecFlux) + hc_avecFlux * (Tg_avecFlux - T3_avecFlux)) :
                          (nombre_couches == 4) ? (hc_avecFlux * (Tg_avecFlux - T3_avecFlux) + hc_avecFlux * (Tg_avecFlux - T4_avecFlux)) :
                          (nombre_couches == 5) ? (hc_avecFlux * (Tg_avecFlux - T4_avecFlux) + hc_avecFlux * (Tg_avecFlux - T5_avecFlux)) : 0;

    cout << "qv(Es) = " << q_v_avecFlux << endl;
    
    /* --------------------------------- Résultat sans flux --------------------------------------*/
    MatrixXd temperatureMatrixSansFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, ems, ems_prime, s);
    VectorXd constantVectorTemperatureSansFlux = createTemperatureVector(0, hc_e, hc_i, Te, Ti, alphaE, s);
    VectorXd temperatureVectorSansFlux = temperatureMatrixSansFlux.inverse() * constantVectorTemperatureSansFlux;

    for (int i = 0; i < nombre_couches; ++i){
        temperatureInitialSansFlux.push_back(temperatureVectorSansFlux[i]);
    }

    VectorXd bilanResultSansFlux = bilan(0, temperatureInitialSansFlux, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, s, type_espace);    

    double T1_sansFlux = bilanResultSansFlux[0];
    double T2_sansFlux = bilanResultSansFlux[1];
    double T3_sansFlux = (nombre_couches >= 3) ? bilanResultSansFlux[2] : 0;
    double T4_sansFlux = (nombre_couches >= 4) ? bilanResultSansFlux[3] : 0;
    double T5_sansFlux = (nombre_couches == 5) ? bilanResultSansFlux[4] : 0;
    double Tg_sansFlux = bilanResultSansFlux[5];
    double hc_sansFlux = bilanResultSansFlux[6];
    double qth_2_sansFlux = bilanResultSansFlux[7];
    double qth_3_sansFlux = bilanResultSansFlux[8];
    double qth_4_sansFlux = bilanResultSansFlux[9];
    double qth_5_sansFlux = bilanResultSansFlux[10];

    cout << "\nT1(0) = " << T1_sansFlux << endl;
    cout << "T2(0) = " << T2_sansFlux << endl;
    cout << "T3(0) = " << T3_sansFlux << endl;
    cout << "T4(0) = " << T4_sansFlux << endl;
    cout << "T5(0) = " << T5_sansFlux << endl;

    double qth_i_sansFlux = (nombre_couches == 2) ? (qth_2_sansFlux - qth_n_prime) : 
                            (nombre_couches == 3) ? (qth_3_sansFlux - qth_n_prime) :
                            (nombre_couches == 4) ? (qth_4_sansFlux - qth_n_prime) :
                            (nombre_couches == 5) ? (qth_5_sansFlux - qth_n_prime) : 0;
    
    cout << "\nqth_i(0) = " << qth_i_sansFlux << endl;

    double qc_i_sansFlux = (nombre_couches == 2) ? (hc_i * (Ti - T2_sansFlux)) : 
                           (nombre_couches == 3) ? (hc_i * (Ti - T3_sansFlux)) :
                           (nombre_couches == 4) ? (hc_i * (Ti - T4_sansFlux)) : 
                           (nombre_couches == 5) ? (hc_i * (Ti - T5_sansFlux)) : 0;
    
    cout << "qc_i(0) = " << qc_i_sansFlux << endl;

    //double Tgjv_sansFlux = bilanResultSansFlux[9];
    
    double q_v_sansFlux = (nombre_couches == 2) ? (hc_sansFlux * (Tg_sansFlux - T1_sansFlux) + hc_sansFlux * (Tg_sansFlux - T2_sansFlux)) : 
                          (nombre_couches == 3) ? (hc_sansFlux * (Tg_sansFlux - T2_sansFlux) + hc_sansFlux * (Tg_sansFlux - T3_sansFlux)) :
                          (nombre_couches == 4) ? (hc_sansFlux * (Tg_sansFlux - T3_sansFlux) + hc_sansFlux * (Tg_sansFlux - T4_sansFlux)) :
                          (nombre_couches == 5) ? (hc_sansFlux * (Tg_sansFlux - T4_sansFlux) + hc_sansFlux * (Tg_sansFlux - T5_sansFlux)) : 0;

    cout << "qv(0) = " << q_v_sansFlux << endl;

    g_th = fabs(qth_i_avecFlux - qth_i_sansFlux) / Es;
    g_c = fabs(qc_i_avecFlux - qc_i_sansFlux) / Es;
    g_v = (type_espace == 1) ? fabs(q_v_avecFlux - q_v_sansFlux) / Es : 0;
    g_tot = g_tot = g_th + g_c + g_v + transmission;

    vector<double> result; 
    result.push_back(g_th);
    result.push_back(g_c);
    result.push_back(g_v);
    result.push_back(g_tot); 

    return result;
}


int main(){
    // Example usage of the createThermalMatrix function
    double Te = 298.15;
    double Ti = 298.15;
    double hc_e = 8;
    double hc_i = 2.5;
    double Es = 500;
    double transmission = 0.03;

    double alphaE1 = 0.3914;
    double alphaE2 = 0.0285;
    double alphaE3 = 0.0578;  
    double alphaE4 = 0.0314;
    double alphaE5 = 0.1251;
    
    double ems1 = 0.837657;
    double ems1_prime = 0.011946; 
    double tau_th1 = 0;   
    double ems2 = 0.837657;
    double ems2_prime = 0.837657;
    double tau_th2 = 0; 
    double ems3 = 0.03552156; 
    double ems3_prime = 0.837657;
    double tau_th3 = 0; 
    double ems4 = 0.837657;
    double ems4_prime = 0.837657;
    double tau_th4 = 0;
    double ems5 = 0.84672;
    double ems5_prime = 0.84672;
    double tau_th5 = 0.03; 

    /* double alphaE1 = 0.3291;
    double alphaE2 = 0.0109;
    double alphaE3 = 0.0091;  
    double alphaE4 = 0.0225; 
    
    double ems1 = 0.837657;
    double ems1_prime = 0.011946; 
    double tau_th1 = 0;   
    double ems2 = 0.837657;
    double ems2_prime = 0.837657;
    double tau_th2 = 0; 
    double ems3 = 0.837657; 
    double ems3_prime = 0.837657;
    double tau_th3 = 0; 
    double ems4 = 0.03552156;
    double ems4_prime = 0.837657;
    double tau_th4 = 0; */
        
    double hauteur = 1.5;
    double largeur = 1.5;
    double s1 = 0.016;
    double s2 = 0.016;
    double s3 = 0.016;
    double s4 = 0.050;

    double T1, T2, T3, T4, T5;

    vector<double> resultat;

    double g_th, g_c, g_v, g_tot;

    vector<double> ems = {ems1, ems2, ems3, ems4, ems5};
    vector<double> ems_prime = {ems1_prime, ems2_prime, ems3_prime, ems4_prime, ems5_prime};
    vector<double> alphaE = {alphaE1, alphaE2, alphaE3, alphaE4, alphaE5};  
    vector<double> tau_th = {tau_th1, tau_th2, tau_th3, tau_th4, tau_th5}; 
    vector<double> s = {s1, s2, s3, s4};               
    
    int type_espace = 1;   
    int position_store  = 1;

    VectorXd bilanResult = bilan(Es, {T1, T2, T3, T4, T5}, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, s, 1);
    resultat = calculFacteurSolaire(Es, transmission, Te, Ti, hc_e, hc_i, ems, ems_prime, alphaE, tau_th, hauteur, largeur, s, type_espace, position_store);

    g_th = resultat[0];
    g_c = resultat[1];
    g_v = resultat[2];
    g_tot = resultat[3];

    g_th = round(g_th * 1000) / 1000;
    g_c = round(g_c * 1000) / 1000;
    g_v = round(g_v * 1000) / 1000;
    g_tot = round(g_tot * 1000) / 1000;

    cout << "\ng_th = " << g_th * 100 << " %" << endl;
    cout << "g_c = " << g_c * 100 << " %" << endl;
    cout << "g_v = " << g_v * 100 << " %" << endl;
    cout << "g_tot = " << g_tot * 100 << " %" << endl;

    return 0;
}