#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "constantes.h"

using namespace Eigen;
using namespace std;

enum Gaz{AIR, ARGON, KRYPTON, XENON, NUM_GAZ};

/* Fonction pour calculer le coefficient "hg" en fonction des propriétés des gazs */
vector<double> calculHg(double T, double dT, double s, const std::vector<double>& fractions, int inclinaison){
    if (fractions.size() != NUM_GAZ) {
        throw std::invalid_argument("Nombre incorrect de fractions de gaz fournies");
    }

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

/* Fonction pour calculer les température du système initiales avec lesquelles on commence le calcul */
MatrixXd CreateTemperatureMatrix(double Te, double Ti, double hc_e, double hc_i, vector<double>ems, vector<double>ems_prime, vector<double>s){
    int nombre_couches = ems.size();
    MatrixXd temperatureMatrix = MatrixXd::Zero(nombre_couches, nombre_couches);

    double deltaT;
    double he, hi;
    double hr_1_2, hr_2_1, hr_2_3, hr_3_2;
    double hg_1, hg_2;
    double s_1, s_2;

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
    
    return temperatureMatrix;
}

VectorXd createTemperatureVector(double Es, double hc_e, double hc_i, double Te, double Ti, vector<double> alphaE, vector<double> s){
    int nombre_couches = alphaE.size();
    VectorXd temperatureVect(nombre_couches);

    double deltaT;
    double he, hi;
    double hr_1_2, hr_2_3, hr_3_2;
    double hg_1, hg_2;
    double s_1 = s[0];
    double s_2 = s[1];

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
    
    return temperatureVect;
}

// Fonction pour créer la matrice des q_th qui va être utilisé ultérieurement pour résoudre le système linéaire afin de déterminer les q_th intiiales
MatrixXd createQthMatrix(vector<double> ems, vector<double> ems_prime, vector<double> tau_th){
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

// Fonction pour créer le vecteurs des constantes qui va être utilisé ultérieurement pour résoudre le système linéaire afin de déterminer les q_th intiiales
VectorXd createQthVector(vector<double> ems, vector<double> ems_prime, vector<double> tau_th,
                            vector<double> q_th, vector<double> q_th_prime, vector<double> temperature){
        
    int nombre_couches = temperature.size();
    VectorXd constantVector;
    //if (ems.size() != ems_prime.size() || ems.size() != ems_prime.size() || ems.size() != tau_th.size()
    //    || ems.size() != q_th.size() || ems.size() != q_th_prime.size() || ems.size() || temperature.size()) {
    //    throw std::invalid_argument("Size Error !");
    //}
    constantVector.resize(temperature.size() * 2);

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

VectorXd calculVitesse(vector<double> T, double Ti, double hauteur, double largeur, double s){
    int nombre_couches = T.size();
    
    double Z1 = 3.1331942;
    double Z2 = 5.3928875;
        
    double T1 = T[0];
    double T2 = T[1];
    double T3 = (nombre_couches == 3) ? T[2] : 0;

    //double T0 = 293.15;
    //double Ti = 293.15;
    double Tgi;
    double rho, muy, lambda, c, hg;

    double Tmp_1_2, Tmp_2_3;
    double deltaT_1_2;
    double deltaT_2_3;

    /* 
    top = top / 1000.0;
    bottom = bottom / 1000.0;
    lateral = lateral / 1000.0;

    As = largeur * s;
    At = largeur * top;
    Ab = largeur * bottom;
    Ah = params.tau_th * largeur * hauteur;
    Al = hauteur * lateral;
    Ar = hauteur * lateral;

    A1eq = fmin(As, (Ab + (Al + Ar + Ah)/4));
    A2eq = fmin(As, (At + (Al + Ar + Ah)/4));

    (A1eq < 0.00001f) ? (A1eq = 0.00001f) : false;
    (A2eq < 0.00001f) ? (A2eq = 0.00001f) : false; 

    Z1  = pow(As/0.6f/A1eq-1, 2.0f);
    Z2  = pow(As/0.6f/A2eq-1, 2.0f);
    */

    if(nombre_couches == 2){
        Tmp_1_2 = (T1 + T2) / 2.0;
        Tgi = (Tmp_1_2 + Ti) / 2.0;
        deltaT_1_2 = fabs(T2 - T1);
        //s = 0.050;
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
        //s = 0.050;
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
            hg = calculHg(Tg, deltaT_1_2, s, {1, 0, 0, 0}, 90)[0]; //check !!!!!!
        }

        else if (nombre_couches == 3)
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

    VectorXd result(6);
    result[0] = hc;
    result[1] = Htp;
    result[2] = Tsor;
    result[3] = Tg;
    result[4] = vitesse;
    result[5] = hg;

    return result;
}

VectorXd bilan(double Es, vector<double> temperatureVector,vector<double> epsilon, vector<double> epsilon_prime, vector<double> alphaE, 
            vector<double> tau_th, double hc_e, double hc_i, double Te, double Ti, double hauteur, double largeur, vector<double> s, int espace){

    int nombre_couches = temperatureVector.size();

    double s1 = s[0];
    double s2 = s[1];

    double T1 = temperatureVector[0];
    double T2 = temperatureVector[1];
    double T3 = (nombre_couches == 3) ? temperatureVector[2] : 0;
    //double Tmp = (T1 + T2) / 2.0;
    double T0 = 293.15;

    double qth_0;
    double qth_1, qth_2, qth_3;
    double qth_0_prime, qth_1_prime;
    double qth_2_prime, qth_3_prime;

    double bilan_1, bilan_2, bilan_3, df1_dT1, df1_dT2, df1_dT3, df2_dT1, df2_dT2, df2_dT3, df3_dT1, df3_dT2, df3_dT3;

    double Tmp_1_2, Tmp_2_3;
    double deltaT_1_2, deltaT_2_3;

    double hg_0 = hc_e;
    double hg_1;
    double hg_2;
    double hg_3;
   
    double hc;
    double hg;
    double Tg;
    double Tsor;
    double Htp;
    double vitesse;
    
    VectorXd v; //vecteur pour stocker les résultat de la fonction calculVitesse()
    MatrixXd thermalMatrix;
    MatrixXd constantVectorQth;
    VectorXd qthVectors;

    double qth_a1, qth_a2, qth_a3;
    double qc_a1, qc_a2, qc_a3;

    // Newton-Raphson iteration loop
    for (int iter = 0; iter < 10; ++iter) {
        if(nombre_couches == 2){
            v = calculVitesse({T1, T2}, Ti, hauteur, largeur, s2);
            hc = v[0];
            Htp = v[1];
            Tsor = v[2];
            Tg = v[3];
            vitesse = v[4];
            hg = v[5];
            
            Tmp_1_2 = (T1 + T2) / 2.0;
            deltaT_1_2 = fabs(T2 - T1);

            hg_1 = calculHg(Tmp_1_2, deltaT_1_2, s2, {1, 0, 0, 0}, 90)[0];
            hg_2 = hc_i;
        
            // Recalculate thermophysical properties
            thermalMatrix = createQthMatrix({epsilon[0], epsilon[1]}, {epsilon_prime[0], epsilon_prime[1]}, {tau_th[0], tau_th[1]});
            constantVectorQth = createQthVector({epsilon[0], epsilon[1]}, {epsilon_prime[0], epsilon_prime[1]}, {tau_th[0], tau_th[1]},
                                                        {qth_0, qth_1, qth_2}, {qth_0_prime, qth_1_prime, qth_2_prime}, {T1, T2});
            qthVectors = thermalMatrix.inverse() * constantVectorQth;

            qth_0 = _StefanBoltzmann * pow(Te,4);
            qth_1 = qthVectors[0];
            qth_0_prime = qthVectors[1]; 
            qth_2 = qthVectors[2];
            qth_2_prime = _StefanBoltzmann * pow(Ti, 4); 
            qth_1_prime = qthVectors[3];

            qth_a1 = epsilon[0] * qth_0 + epsilon_prime[0] * qth_1_prime  - (epsilon[0] + epsilon_prime[0]) * _StefanBoltzmann * pow(T1, 4);
            qth_a2 = epsilon[1] * qth_1 + epsilon_prime[1] * qth_2_prime  - (epsilon[1] + epsilon_prime[1]) * _StefanBoltzmann * pow(T2, 4);
            
            if(espace == 0){
                qc_a1 = hg_0 * (T0 - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (Ti - T2);
            }

            else if (espace == 1){
                qc_a1 = hc_e * (Te - T1) + hc * (Tg - T1);
                qc_a2 = hc * (Tg - T2) + hc_i * (Ti - T2);
            }
            // Recalculate bilan equations
            bilan_1 = alphaE[0] * Es + epsilon[0] * qth_0 + epsilon_prime[0] * qth_1_prime  - (epsilon[0] + epsilon_prime[0]) * _StefanBoltzmann * pow(T1, 4) + hc_e * (Te - T1) + hc * (Tg - T1);
            bilan_2 = alphaE[1] * Es + epsilon[1] * qth_1 + epsilon_prime[1] * qth_2_prime  - (epsilon[1] + epsilon_prime[1]) * _StefanBoltzmann * pow(T2, 4) + hc * (Tg - T2) + hc_i * (Ti - T2);

            VectorXd equation(2);
            equation[0] = bilan_1;
            equation[1] = bilan_2;

            // Update Jacobian matrix
            df1_dT1 = -4.0 * (epsilon[0] + epsilon_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hc_e - hc;
            df1_dT2 = 0;
            df2_dT1 = 0;
            df2_dT2 = -4.0 * (epsilon[1] + epsilon_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hc - hc_i;

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
            Htp = v[1];
            Tsor = v[2];
            Tg = v[3];
            vitesse = v[4];
            hg = v[5];

            Tmp_1_2 = (T1 + T2) / 2;
            Tmp_2_3 = (T2 + T3) / 2;
            deltaT_1_2 = fabs(T2 - T1);
            deltaT_2_3 = fabs(T3 - T2);

            hg_1 = calculHg(Tg, deltaT_2_3, s1, {1, 0, 0, 0}, 90)[0];
            hg_2 = calculHg(Tmp_1_2, deltaT_1_2, s2, {0.1, 0.9, 0, 0}, 90)[0];
            hg_3 = hc_i;

            thermalMatrix = createQthMatrix({epsilon[0], epsilon[1], epsilon[2]}, {epsilon_prime[0], epsilon_prime[1], epsilon_prime[2]}, {tau_th[0], tau_th[1], tau_th[2]});
            constantVectorQth = createQthVector({epsilon[0], epsilon[1], epsilon[2]}, {epsilon_prime[0], epsilon_prime[1], epsilon_prime[2]}, {tau_th[0], tau_th[1], tau_th[2]},
                                                {qth_0, qth_1, qth_2, qth_3}, {qth_0_prime, qth_1_prime, qth_2_prime, qth_3_prime}, {T1, T2, T3});
            qthVectors = thermalMatrix.inverse() * constantVectorQth;

            qth_0 = _StefanBoltzmann * pow(Te,4);
            qth_1 = qthVectors[0];
            qth_0_prime = qthVectors[1];
            qth_2 = qthVectors[2];
            qth_1_prime = qthVectors[3];
            qth_3 = qthVectors[4];
            qth_2_prime = qthVectors[5];
            qth_3_prime = _StefanBoltzmann * pow(Ti, 4);

            qth_a1 = epsilon[0] * qth_0 + epsilon_prime[0] * qth_1_prime  - (epsilon[0] + epsilon_prime[0]) * _StefanBoltzmann * pow(T1, 4);
            qth_a2 = epsilon[1] * qth_1 + epsilon_prime[1] * qth_2_prime  - (epsilon[1] + epsilon_prime[1]) * _StefanBoltzmann * pow(T2, 4);
            qth_a3 = epsilon[2] * qth_2 + epsilon_prime[2] * qth_3_prime  - (epsilon[2] + epsilon_prime[2]) * _StefanBoltzmann * pow(T3, 4);

            // CHECK !!!
            if(espace == 0){ // Espace fermé
                qc_a1 = hc_e * (Te - T1) + hg_1 * (T2 - T1);
                qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
                qc_a3 = hg_2 * (T2 - T3) + hc_i * (Ti - T3);
            }
            
            else if(espace == 1){ // Espace ouvert 
                qc_a1 = hc_e * (Te - T1) + hg_2 * (Tmp_1_2 - T1);
                qc_a2 = hg_2 * (Tmp_1_2 - T2) + hc * (Tg - T2);
                qc_a3 = hc * (Tg - T3) + hc_i * (Ti - T3);
            }
        
            // Recalculate bilan equations
            bilan_1 = alphaE[0] * Es + epsilon[0] * qth_0 + epsilon_prime[0] * qth_1_prime  - (epsilon[0] + epsilon_prime[0]) * _StefanBoltzmann * pow(T1, 4) + hc_e * (Te - T1) + hg_1 * (Tmp_1_2 - T1);
            bilan_2 = alphaE[1] * Es + epsilon[1] * qth_1 + epsilon_prime[1] * qth_2_prime  - (epsilon[1] + epsilon_prime[1]) * _StefanBoltzmann * pow(T2, 4) + hg_1 * (Tmp_1_2 - T2) + hc * (Tg - T2);
            bilan_3 = alphaE[2] * Es + epsilon[2] * qth_2 + epsilon_prime[2] * qth_3_prime  - (epsilon[2] + epsilon_prime[2]) * _StefanBoltzmann * pow(T3, 4) + hc * (Tg - T3) + hc_i * (Ti - T3);

            VectorXd equation(3);
            equation[0] = bilan_1;
            equation[1] = bilan_2;
            equation[2] = bilan_3;

            // Update Jacobian matrix
            df1_dT1 = -4.0 * (epsilon[0] + epsilon_prime[0]) * _StefanBoltzmann * pow(T1, 3) - hc_e - hg_1;
            df1_dT2 = hg_1;
            df1_dT3 = 0;
            df2_dT1 = hg_1;
            df2_dT2 = -4.0 * (epsilon[1] + epsilon_prime[1]) * _StefanBoltzmann * pow(T2, 3) - hg_1 - hc;
            df2_dT3 = 0;
            df3_dT1 = 0;
            df3_dT2 = 0;
            df3_dT3 = -4.0 * (epsilon[2] + epsilon_prime[2]) * _StefanBoltzmann * pow(T3, 3) - hc - hc_i;

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
    
    VectorXd result = (nombre_couches == 2) ? VectorXd(13) : VectorXd(15);
    result[0] = T1;
    result[1] = T2;
    result[2] = (nombre_couches == 3) ? T3 : 0;
    result[3] = hc;
    result[4] = Htp;
    result[5] = Tsor;
    result[6] = Tg;
    result[7] = vitesse;
    result[8] = hg;
    result[9] = qth_1;
    result[10] = qth_0_prime;
    result[11] = qth_2;
    result[12] = qth_1_prime;
    if (nombre_couches == 3) {
        result[13] = qth_3;
        result[14] = qth_2_prime;
    }
    return result;
}

vector<double> calculFacteurSolaire(double Es, double transmission, double Te, double Ti, double hc_e, double hc_i, vector<double> ems, vector<double> ems_prime, 
                                    vector<double> alphaE, vector<double> tau_th, double hauteur, double largeur, vector<double> s, int type_espace, int position_store){
    
    int nombre_couches = ems.size();
    vector<double> temperatureInitialAvecFlux, temperatureInitialSansFlux;
    double g_th, g_c, g_v, g_tot;

    double s_1 = s[0];
    double s_2 = s[1];

    double qth_0 = _StefanBoltzmann * pow(Te,4);
    double qth_n_prime = _StefanBoltzmann * pow(Ti,4);

    /*                                  resultat avec flux                                  */
    MatrixXd temperatureMatrixAvecFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, ems, ems_prime, {s_1, s_2});
    VectorXd constantVectorTemperatureAvecFlux = createTemperatureVector(Es, hc_e, hc_i, Te, Ti, alphaE, {s_1, s_2});
    VectorXd temperatureVectorAvecFlux = temperatureMatrixAvecFlux.inverse() * constantVectorTemperatureAvecFlux;

    for (int i = 0; i < nombre_couches; ++i) {
        temperatureInitialAvecFlux.push_back(temperatureVectorAvecFlux[i]);
    }

    /* double T1_initial = temperatureVector[0];
    double T2_initial = temperatureVector[1];
    double T3_initial = (nombre_couches == 3) ? temperatureVector[2] : 0; */

    VectorXd bilanResultAvecFlux = bilan(Es, temperatureInitialAvecFlux, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, {s_1, s_2}, type_espace);    

    double T1_avecFlux = bilanResultAvecFlux[0];
    double T2_avecFlux = bilanResultAvecFlux[1];
    double T3_avecFlux = (nombre_couches == 3) ? bilanResultAvecFlux[2] : 0;
    double hc_avecFlux = bilanResultAvecFlux[3];
    double Htp_avecFlux = bilanResultAvecFlux[4];
    double Tsor_avecFlux = bilanResultAvecFlux[5];
    double Tg_avecFlux = bilanResultAvecFlux[6];
    double vitesseAvecFlux = bilanResultAvecFlux[7];
    double hg_avecFlux = bilanResultAvecFlux[8];

    double qth_1_avecFlux = bilanResultAvecFlux[9];
    double qth_0_prime_avecFlux = bilanResultAvecFlux[10];
    double qth_2_avecFlux = bilanResultAvecFlux[11];
    double qth_1_prime_avecFlux = bilanResultAvecFlux[12];

    double qth_e_avecFlux = qth_0_prime_avecFlux - qth_0;
    double qth_i_avecFlux = qth_2_avecFlux - qth_n_prime;

    double qc_e_avecFlux = hc_e * (T1_avecFlux - Te);
    double qc_i_avecFlux = hc_i * (Ti - T2_avecFlux);

    //double q_v_Es = (type_espace == 1 && position_store ==  1) ? hc_avecFlux * (Tg_avecFlux - T1_avecFlux) + hc_avecFlux * (Tg_avecFlux - T2_avecFlux) : 0;
    double q_v_avecFlux = hc_avecFlux * (Tg_avecFlux - T1_avecFlux) + hc_avecFlux * (Tg_avecFlux - T2_avecFlux);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /*                                  resultat sans flux                                  */
    MatrixXd temperatureMatrixSansFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, ems, ems_prime, {s_1, s_2});
    VectorXd constantVectorTemperatureSansFlux = createTemperatureVector(0, hc_e, hc_i, Te, Ti, alphaE, {s_1, s_2});
    VectorXd temperatureVectorSansFlux = temperatureMatrixSansFlux.inverse() * constantVectorTemperatureSansFlux;

    /* double T1_initial_sansFlux = temperatureVectorSansFlux[0];
    double T2_initial_sansFlux = temperatureVectorSansFlux[1];
    double T3_initial_sansFlux = (nombre_couches == 3) ? temperatureVectorSansFlux[2] : 0; */

    for (int i = 0; i < nombre_couches; ++i) {
        temperatureInitialSansFlux.push_back(temperatureVectorSansFlux[i]);
    }

    VectorXd bilanResultSansFlux = bilan(0, temperatureInitialSansFlux, ems, ems_prime, alphaE, tau_th, hc_e, hc_i, Te, Ti, hauteur, largeur, {s_1, s_2}, type_espace);    

    double T1_sansFlux = bilanResultSansFlux[0];
    double T2_sansFlux = bilanResultSansFlux[1];
    double T3_sansFlux = (nombre_couches == 3) ? bilanResultSansFlux[2] : 0;
    double hc_sansFlux = bilanResultSansFlux[3];
    double Htp_sansFlux = bilanResultSansFlux[4];
    double Tsor_sansFlux = bilanResultSansFlux[5];
    double Tg_sansFlux = bilanResultSansFlux[6];
    double vitesse_sansFlux = bilanResultSansFlux[7];
    double hg_sansFlux = bilanResultSansFlux[8];

    double qth_1_sansFlux = bilanResultSansFlux[9];
    double qth_0_prime_sansFlux = bilanResultSansFlux[10];
    double qth_2_sansFlux = bilanResultSansFlux[11];
    double qth_1_prime_sansFlux = bilanResultSansFlux[12];

    double qth_e_sansFlux = qth_0_prime_sansFlux - qth_0;
    double qth_i_sansFlux = qth_2_sansFlux - qth_n_prime;

    double qc_e_sansFlux = hc_e * (T1_sansFlux - Te);
    double qc_i_sansFlux = hc_i * (Ti - T2_sansFlux);

    //double q_v_sansFlux = (type_espace == 1 && position_store ==  1) ? hc_sansFlux * (Tg_sansFlux - T1_sansFlux) + hc_sansFlux * (Tg_sansFlux - T2_sansFlux) : 0;
    double q_v_sansFlux = hc_sansFlux * (Tg_sansFlux - T1_sansFlux) + hc_sansFlux * (Tg_sansFlux - T2_sansFlux);

    g_th = fabs(qth_i_avecFlux - qth_i_sansFlux) / Es;
    g_c = fabs(qc_i_avecFlux - qc_i_sansFlux) / Es;
    g_v = fabs(q_v_avecFlux - q_v_sansFlux) / Es;
    //g_v = (type_espace == 1) ? fabs(q_v_Es - q_v_sansFlux) / Es : 0;
    g_tot = g_tot = g_th + g_c + g_v + transmission;

    vector<double> result; 
    result.push_back(g_th);
    result.push_back(g_c);
    result.push_back(g_v);
    result.push_back(g_tot);

    return result;
}


int main() {
    // Example usage of the createThermalMatrix function
    double Te = 298.15;
    double Ti = 298.15;
    double hc_e = 8;
    double hc_i = 2.5;
    double alphaE1 = 0.0839;
    double alphaE2 = 0.1055;
    double alphaE3 = 0.0107;    
    double Es = 500;
    double transmission = 0.166;

    double epsilon1 = 0.837657;
    double epsilon1_prime = 0.837657; 
    double tau_th1 = 0.0;   
    double epsilon2 = 0.837657;
    double epsilon2_prime = 0.837657;
    double tau_th2 = 0.03; 
    double epsilon3 = 0.84672; 
    double epsilon3_prime = 0.84672;
    double tau_th3 = 0.04;  
        
    double hauteur = 1.4;
    double largeur = 2.4;
    double s1 = 0.016;
    double s2 = 0.050;

    vector<double> resultat;

    double g_th, g_c, g_v, g_tot;

    vector<double> ems = {epsilon1, epsilon2};
    vector<double> ems_prime = {epsilon1_prime, epsilon2_prime};
    vector<double> alphaE = {alphaE1, alphaE2};  
    vector<double> tau_th = {tau_th1, tau_th2}; 
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

    cout << "\n---------------------------- Final Result ----------------------------" << endl;
    cout << "g_th = " << g_th * 100.0 << " %" << endl;
    cout << "g_c = " << g_c * 100.0 << " %" << endl;
    cout << "g_v = " << g_v * 100.0 << " %" << endl;
    cout << "g_tot = " << g_tot * 100.0  << " %" <<endl;  

    return 0;
}