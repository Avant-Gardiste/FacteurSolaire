#include <Eigen/Dense>
#include <iostream>
#include <vector>
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

MatrixXd CreateTemperatureMatrix(double Te, double Ti, double hc_e, double hc_i, vector<double>ems, vector<double>ems_prime, vector<double>s)
{
    int nombre_couches = ems.size();
    MatrixXd temperatureMatrix = MatrixXd::Zero(nombre_couches, nombre_couches);

    double deltaT;
    double he, hi;
    double hr_1_2, hr_2_1, hr_2_3, hr_3_2;

    he = hc_e + 4.0*_StefanBoltzmann*pow(Te,3);
    hi = hc_i + 4.0*_StefanBoltzmann*pow(Ti, 3);
    deltaT = 15 / (nombre_couches-1);

    vector<double> hg(nombre_couches, 0.0);  // Vector to store hg values for each layer
    for (int i = 0; i < nombre_couches; ++i) {
        if (i < s.size()) {  // Ensure s has enough elements to avoid accessing out of bounds
            hg[i] = calculHg(_T0, deltaT, s[i], {1, 0, 0, 0}, 90)[0];
        }
    }

    if(nombre_couches==2)
    {   
        hr_1_2 = 4.0 * (1 / (1/ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_1 = 4.0 * (1 / (1/ems_prime[1] + 1/ems[0] - 1)) * _StefanBoltzmann * pow(_T0,3); 
        temperatureMatrix(0,0) = -(he+hr_1_2+hg[0]);
        temperatureMatrix(1,0) = hr_2_1;
        temperatureMatrix(0,1) = hr_1_2;
        temperatureMatrix(1,1) = -(hr_2_1+hi+hg[0]);
    }

    else if (nombre_couches==3)
    {   
        hr_1_2 = 4.0 * (1 / (1/ems_prime[0] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_2_3 = 4.0 * (1 / (1/ems_prime[1] + 1/ems[2] - 1)) * _StefanBoltzmann * pow(_T0,3);
        hr_3_2 = 4.0 * (1 / (1/ems_prime[2] + 1/ems[1] - 1)) * _StefanBoltzmann * pow(_T0,3);

        temperatureMatrix(0,0) = -(he+hr_1_2+hg[0]);
        temperatureMatrix(1,0) = hr_1_2;
        temperatureMatrix(2,0) = 0.0;

        temperatureMatrix(0,1) = hr_1_2;
        temperatureMatrix(1,1) = -(hr_1_2+hr_2_3+hg[0]+hg[1]);
        temperatureMatrix(2,1) = hr_3_2;

        temperatureMatrix(0,2) = 0.0;
        temperatureMatrix(1,2) = hr_2_3;
        temperatureMatrix(2,2) = -(hi+hr_3_2+hg[1]);
    }

    /* 
    else if(nombre_couches == 4){
        
    }

    else if(nombre_couches == 5){

    }  
    */

    return temperatureMatrix;

}

VectorXd createTemperatureVector(double Es, double hc_e, double hc_i, double Te, double Ti, vector<double> alphaE, vector<double> s)
{   
    int nombre_couches = alphaE.size();
    VectorXd temperatureVect(nombre_couches);

    double deltaT;
    double he, hi;
    double hr_1_2, hr_2_3, hr_3_2;
    

    he = hc_e + 4.0*_StefanBoltzmann*pow(Te,3);
    hi = hc_i + 4.0*_StefanBoltzmann*pow(Ti, 3);
    deltaT = 15 / (nombre_couches-1);

    vector<double> hg(nombre_couches, 0.0);  // Vector to store hg values for each layer
    for (int i = 0; i < nombre_couches; ++i) {
        if (i < s.size()) {  // Ensure s has enough elements to avoid accessing out of bounds
            hg[i] = calculHg(_T0, deltaT, s[i], {1, 0, 0, 0}, 90)[0];
        }
    }

    if(nombre_couches==2){
        temperatureVect << -alphaE[0]*Es - he*Te - hg[0]*_T0,
                        -alphaE[1]*Es - hi*Ti - hg[0]*_T0;
    }
    else if(nombre_couches==3){
        temperatureVect << -alphaE[0]*Es - he*Te - hg[0]*_T0,
                        -alphaE[1]*Es - hg[0]*_T0 - hg[1]*_T0,
                        -alphaE[2]*Es - hi*Ti - hg[1]*_T0;
    }
    
    /* 
    else if(nombre_couches==4){

    }
    
    else if(nombre_couches==5){

    }
     
    */
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
    
    /* double top, bottom, lateral;
    double As, At, Ab, Ah, Al, Ar;
    double A1eq, A2eq;
    double Z1, Z2;

    double tau_th_effective = 0;
    for(double value : tau_th){
        if(value != 0){
            tau_th_effective = value;
            break;
        }
    }

    top = top / 1000.0;
    bottom = bottom / 1000.0;
    lateral = lateral / 1000.0;

    As = largeur * s;
    At = largeur * top;
    Ab = largeur * bottom;
    Ah = tau_th_effective * largeur * hauteur;
    Al = hauteur * lateral;
    Ar = hauteur * lateral;

    A1eq = fmin(As, (Ab + (Al + Ar + Ah)/4));
    A2eq = fmin(As, (At + (Al + Ar + Ah)/4));

    (A1eq < 0.00001f) ? (A1eq = 0.00001f) : false;
    (A2eq < 0.00001f) ? (A2eq = 0.00001f) : false; 

    Z1  = pow(As/0.6f/A1eq-1, 2.0f);
    Z2  = pow(As/0.6f/A2eq-1, 2.0f); */

    double Z1 = 0.44444444; //3.1331942;
    double Z2 = 0.44444444; //5.3928875;
        
    double T1 = T[0];
    double T2 = T[1];
    double T3 = (nombre_couches == 3) ? T[2] : 0;

    double Tgi;
    double rho, muy, lambda, c, hg;

    double Tmp_1_2, Tmp_2_3;
    double deltaT_1_2;
    double deltaT_2_3;

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

    double T1 = temperatureVector[0];
    double T2 = temperatureVector[1];
    double T3 = (nombre_couches == 3) ? temperatureVector[2] : 0;

    double qth_0;
    double qth_1, qth_2, qth_3;
    double qth_0_prime, qth_1_prime;
    double qth_2_prime, qth_3_prime;

    double bilan_1, bilan_2, bilan_3, df1_dT1, df1_dT2, df1_dT3, df2_dT1, df2_dT2, df2_dT3, df3_dT1, df3_dT2, df3_dT3;

    double Tmp_1_2, Tmp_2_3;
    double deltaT_1_2, deltaT_2_3;

    double hg_0;
    double hg_1;
    double hg_2;
    double hg_3;
   
    double hc;
    double Tg;
    
    VectorXd v; //vecteur pour stocker les résultat de la fonction calculVitesse()
    MatrixXd qthConstantMatrix;
    MatrixXd qthConstantVector;
    VectorXd qthVector;

    double qth_a1, qth_a2, qth_a3;
    double qc_a1, qc_a2, qc_a3;

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
                qc_a1 = hg_0 * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hg_3 * (Ti - T3);

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
        
    }
    
    VectorXd result = VectorXd(7);
    result[0] = T1;
    result[1] = T2;
    result[2] = (nombre_couches == 3) ? T3 : 0;
    result[3] = Tg;
    result[4] = hc;
    result[5] = qth_2;
    result[6] = (nombre_couches == 3) ? qth_3 : 0;
    return result;
}

vector<double> calculFacteurSolaire(double Es, double transmission, double Te, double Ti, double hc_e, double hc_i, vector<double> ems, vector<double> ems_prime, 
                                    vector<double> alphaE, vector<double> tau_th, double hauteur, double largeur, vector<double> s, int type_espace, int position_store)
{   
    int nombre_couches = ems.size();
    vector<double> temperatureInitialAvecFlux, temperatureInitialSansFlux;
    double g_th, g_c, g_v, g_tot;

    double s_1 = s[0];
    double s_2 = s[1];

    double qth_0 = _StefanBoltzmann * pow(Te,4);
    double qth_n_prime = _StefanBoltzmann * pow(Ti,4);

    /*------------------------------------------------------- resultat avec flux -------------------------------------------------------*/
    MatrixXd temperatureMatrixAvecFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, ems, ems_prime, {s_1, s_2});
    VectorXd constantVectorTemperatureAvecFlux = createTemperatureVector(Es, hc_e, hc_i, Te, Ti, alphaE, {s_1, s_2});
    VectorXd temperatureVectorAvecFlux = temperatureMatrixAvecFlux.inverse() * constantVectorTemperatureAvecFlux;

    for (int i = 0; i < nombre_couches; ++i) {
        temperatureInitialAvecFlux.push_back(temperatureVectorAvecFlux[i]);
    }

    VectorXd bilanResultAvecFlux = bilan(Es, temperatureInitialAvecFlux, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, {s_1, s_2}, type_espace);    

    double T1_avecFlux = bilanResultAvecFlux[0];
    double T2_avecFlux = bilanResultAvecFlux[1];
    double T3_avecFlux = (nombre_couches == 3) ? bilanResultAvecFlux[2] : 0;
    double Tg_avecFlux = bilanResultAvecFlux[3];
    double hc_avecFlux = bilanResultAvecFlux[4];

    cout << "\nT1 = " << T1_avecFlux << endl;
    cout << "T2 = " << T2_avecFlux << endl;
    cout << "T3 = " << T3_avecFlux << endl;
    
    double qth_2_avecFlux = bilanResultAvecFlux[5];
    double qth_3_avecFlux = (nombre_couches == 3) ? bilanResultAvecFlux[6] : 0;

    double qth_i_avecFlux = (temperatureInitialAvecFlux.size() == 2) ? (qth_2_avecFlux - qth_n_prime) : (qth_3_avecFlux - qth_n_prime);
    cout << "\nqth_i = " << qth_i_avecFlux << endl;

    double qc_i_avecFlux = (temperatureInitialAvecFlux.size() == 2) ? (hc_i * (Ti - T2_avecFlux)) : (hc_i * (Ti - T3_avecFlux));
    cout << "qc_i = " << qc_i_avecFlux << endl;
    
    double q_v_avecFlux = (temperatureInitialAvecFlux.size() == 2) ? (hc_avecFlux * (Tg_avecFlux - T1_avecFlux) + hc_avecFlux * (Tg_avecFlux - T2_avecFlux)) : (hc_avecFlux * (T2_avecFlux - Tg_avecFlux) + hc_avecFlux * (T3_avecFlux - Tg_avecFlux));
    cout << "\nqv = " << q_v_avecFlux << endl;
    /*------------------------------------------------------- resultat sans flux -------------------------------------------------------*/
    MatrixXd temperatureMatrixSansFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, ems, ems_prime, {s_1, s_2});
    VectorXd constantVectorTemperatureSansFlux = createTemperatureVector(0, hc_e, hc_i, Te, Ti, alphaE, {s_1, s_2});
    VectorXd temperatureVectorSansFlux = temperatureMatrixSansFlux.inverse() * constantVectorTemperatureSansFlux;

    for (int i = 0; i < nombre_couches; ++i) {
        temperatureInitialSansFlux.push_back(temperatureVectorSansFlux[i]);
    }

    VectorXd bilanResultSansFlux = bilan(0, temperatureInitialSansFlux, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, {s_1, s_2}, type_espace);    

    double T1_sansFlux = bilanResultSansFlux[0];
    double T2_sansFlux = bilanResultSansFlux[1];
    double T3_sansFlux = (nombre_couches == 3) ? bilanResultSansFlux[2] : 0;
    double Tg_sansFlux = bilanResultSansFlux[3];
    double hc_sansFlux = bilanResultSansFlux[4];
    
    double qth_2_sansFlux = bilanResultSansFlux[5];
    double qth_3_sansFlux = (nombre_couches == 3) ? bilanResultSansFlux[6] : 0;

    double qth_i_sansFlux = (temperatureInitialSansFlux.size() == 2) ? (qth_2_sansFlux - qth_n_prime) : (qth_3_sansFlux - qth_n_prime);

    double qc_i_sansFlux = (temperatureInitialSansFlux.size() == 2) ? (hc_i * (Ti - T2_sansFlux)) : (hc_i * (Ti - T3_sansFlux));
    
    double q_v_sansFlux = (temperatureInitialSansFlux.size() == 2) ? (hc_sansFlux * (Tg_sansFlux - T1_sansFlux) + hc_sansFlux * (Tg_sansFlux - T2_sansFlux)) : (hc_sansFlux * (T2_sansFlux - Tg_sansFlux) + hc_sansFlux * (T3_sansFlux - Tg_sansFlux));
    
    g_th = fabs(qth_i_avecFlux - qth_i_sansFlux) / Es;
    g_c = fabs(qc_i_avecFlux - qc_i_sansFlux) / Es;
    g_v = (position_store == 1) ? fabs(q_v_avecFlux - q_v_sansFlux) / Es : 0;
    g_tot = g_tot = g_th + g_c + g_v + transmission;

    vector<double> result; 
    result.push_back(g_th);
    result.push_back(g_c);
    result.push_back(g_v);
    result.push_back(g_tot);

    return result;
}

int main() {
    // Example usage 
    double Te = 298.15;
    double Ti = 298.15;
    double hc_e = 8;
    double hc_i = 2.5;
    double Es = 500;
    double transmission = 0.006;
    
    //double alphaE1 = 0.4887;
    //double alphaE2 = 0.0116;
    //double alphaE3 = 0.0107;   

    /* double alphaE1 = 0.0839;
    double alphaE2 = 0.1055;
    double alphaE3 = 0.5941;    */

    double alphaE1 = 0.0979;
    double alphaE2 = 0.1488;
    double alphaE3 = 0.5941;    
    

    /* double ems1 = 0.837657;
    double ems1_prime = 0.837657; 
    double tau_th1 = 0.0;   
    double ems2 = 0.837657;
    double ems2_prime = 0.837657;
    double tau_th2 = 0.03; 
    double ems3 = 0.84672; 
    double ems3_prime = 0.84672;
    double tau_th3 = 0.04;  */
 
    double ems1 = 0.837;
    double ems1_prime = 0.837; 
    double tau_th1 = 0.0;   
    double ems2 = 0.837;
    double ems2_prime = 0.837;
    double tau_th2 = 0.0; 
    double ems3 = 0.837; 
    double ems3_prime = 0.837;
    double tau_th3 = 0.02;   
        
    double hauteur = 2;
    double largeur = 1;
    double s1 = 0.018;
    double s2 = 0.100;

    vector<double> resultat;

    double g_th, g_c, g_v, g_tot;

    vector<double> ems = {ems1, ems2, ems3};
    vector<double> ems_prime = {ems1_prime, ems2_prime, ems3_prime};
    vector<double> alphaE = {alphaE1, alphaE2, alphaE3};  
    vector<double> tau_th = {tau_th1, tau_th2, tau_th3}; 
    vector<double> s = {s1, s2};               
    
    int type_espace = 1;   
    int position_store  = 1;

    resultat = calculFacteurSolaire(Es, transmission, Te, Ti, hc_e, hc_i, ems, ems_prime, alphaE, tau_th, hauteur, largeur, s, type_espace, position_store);

    g_th = resultat[0];
    g_c = resultat[1];
    g_v = resultat[2];
    g_tot = resultat[3];

    g_th = round(g_th * 1000.0) / 1000.0;
    g_c = round(g_c * 1000.0) / 1000.0;
    g_v = round(g_v * 1000.0) / 1000.0; 
    g_tot = round(g_tot * 1000.0) / 1000.0;  

    cout << "---------------------------- Final Result ----------------------------" << endl;
    cout << "g_th = " << g_th * 100.0 << " %" << endl;
    cout << "g_c = " << g_c * 100.0 << " %" << endl;
    cout << "g_v = " << g_v * 100.0 << " %" << endl;
    cout << "g_tot = " << g_tot * 100.0  << " %" <<endl;  

    return 0;
}
