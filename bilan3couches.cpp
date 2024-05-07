#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

enum Gaz{AIR, ARGON, KRYPTON, XENON, NUM_GAZ};

vector<double> calculHg(double T, double dT, double s, const std::vector<double>& fractions) {
    if (fractions.size() != NUM_GAZ) {
        throw std::invalid_argument("Nombre incorrect de fractions de gaz fournies");
    }

    // Propriétés de base pour chaque gaz
    std::vector<double> rho(NUM_GAZ), muy(NUM_GAZ), lambda(NUM_GAZ), c(NUM_GAZ);
    rho[AIR] = 1.189 - 0.0044 * (T - 293.15);
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
    double final_rho = 0, final_muy = 0, final_lambda = 0, final_c = 0;
    for (int i = 0; i < NUM_GAZ; i++) {
        final_rho += rho[i] * fractions[i];
        final_muy += muy[i] * fractions[i];
        final_lambda += lambda[i] * fractions[i];
        final_c += c[i] * fractions[i];
    }

    // Calcul de Grashof, Prandtl, et Nusselt
    double Gr = 9.81 * pow(s, 3) * dT * pow(final_rho, 2) / (T * pow(final_muy, 2));
    double Pr = final_muy * final_c / final_lambda;
    double Nu = fmax(0.035 * pow(Gr * Pr, 0.38), 1);  // Coefficient ajusté pour une meilleure corrélation

    // Calcul de hg
    double hg = Nu * final_lambda / s;

    std::vector<double> results;
    results.push_back(hg);   // Index 0
    results.push_back(final_rho);  // Index 1
    results.push_back(final_muy);  // Index 2
    results.push_back(final_c);

    return results;
    //return hg;
}

MatrixXd CreateTemperatureMatrix(double Te, double Ti, double hc_e, double hc_i, vector<double> ems, vector<double> ems_prime){
    double sigma = 5.67e-8;
    double T0 = 293.15;
    double deltaT = 2.5;
    double A = 0.035;
    double n = 0.38;

    double he, hi;
    double hr_1_2, hr_2_3, hr_3_2;
    double hg_1, hg_2;

    //double Tm = (T0 + Te) / 2;
    double s_1 = 0.018;
    double s_2 = 0.100;


    hg_1 = calculHg(T0, deltaT, s_1, {1, 0, 0, 0})[0];
    cout << "hg_1 = " << hg_1 << endl;

    hg_2 = calculHg(T0, deltaT, s_2, {1, 0, 0, 0})[0];
    cout << "hg_2 = " << hg_2 << endl;

    he = hc_e + 4.0*sigma*pow(Te,3);
    hi = hc_i + 4.0*sigma*pow(Ti, 3);
    cout << "\nhe = " << he << endl;
    cout << "hi = " << hi << endl;

    hr_1_2 = 4.0 * (1 / (1 / ems_prime[0] + 1/ems[1] - 1)) * sigma * pow(T0,3);
    hr_2_3 = 4.0 * (1 / (1 / ems_prime[1] + 1/ems[2] - 1)) * sigma * pow(T0,3);
    hr_3_2 = 4.0 * (1 / (1 / ems_prime[2] + 1/ems[1] - 1)) * sigma * pow(T0,3);

    cout << "\nhr_1_2 = " << hr_1_2 << endl;
    cout << "hr_2_3 = " << hr_2_3 << endl;
    cout << "hr_3_2 = " << hr_3_2 << endl;

    MatrixXd temperatureMatrix = MatrixXd::Zero(3, 3);
    
    //int size = ems.size() * 2;
    //MatrixXd thermalMatrix = MatrixXd::Zero(size, size);
    //thermalMatrix.diagonal().setOnes();

    temperatureMatrix(0,0) = -(he+hr_1_2+hg_1);
    temperatureMatrix(1,0) = hr_1_2;
    temperatureMatrix(2,0) = 0.0;

    temperatureMatrix(0,1) = hr_1_2;
    temperatureMatrix(1,1) = -(hr_1_2+hr_2_3+hg_1+hg_2);
    temperatureMatrix(2,1) = hr_3_2;

    temperatureMatrix(0,2) = 0.0;
    temperatureMatrix(1,2) = hr_2_3;
    temperatureMatrix(2,2) = -(hi+hr_3_2+hg_2);


    return temperatureMatrix;
}

VectorXd createTemperatureVector(double Es, double hc_e, double hc_i, double Te, double Ti, vector<double> alphaE){
    VectorXd temperatureVect(3);
    double sigma = 5.67e-8;
    double T0 = 293.15;
    double deltaT = 2.5;
   
    double he, hi;
    double hr_1_2, hr_2_3, hr_3_2;
    double hg_1, hg_2;
    double s_1 = 0.018;
    double s_2 = 0.100;

    hg_1 = calculHg(T0, deltaT, s_1, {1, 0, 0, 0})[0];
    cout << "hg_1 = " << hg_1 << endl;

    hg_2 = calculHg(T0, deltaT, s_2, {1, 0, 0, 0})[0];
    cout << "hg_2 = " << hg_2 << endl;

    he = hc_e + 4.0*sigma*pow(Te,3);
    hi = hc_i + 4.0*sigma*pow(Ti, 3);

    temperatureVect << -alphaE[0]*Es - he*Te - hg_1*T0,
                    -alphaE[1]*Es - hg_1*T0 - hg_2*T0,
                    -alphaE[2]*Es - hi*Ti - hg_2*T0;
    
    return temperatureVect;
}

MatrixXd createQthMatrix(vector<double> ems, vector<double> ems_prime, vector<double> tau_th){
        if(ems.size() != ems_prime.size() || ems.size() != tau_th.size()){
            throw std::invalid_argument("Size Error");
        }

        int size = ems.size() * 2;
        MatrixXd thermalMatrix = MatrixXd::Zero(size, size);
        thermalMatrix.diagonal().setOnes();

        for(int i=0; i<ems.size(); ++i){
            if(ems.size()==2)
            {
                thermalMatrix(0, 3) = -(1 - ems_prime[0] - tau_th[0]);
                thermalMatrix(1, 3) = -tau_th[0];
                thermalMatrix(2, 0) = -tau_th[1];
                thermalMatrix(3, 0) = -(1 - ems[1] - tau_th[1]);
            }
            else if(ems.size()==3)
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
            else if (ems.size()==4)
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
            else if (ems.size()==5)
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

VectorXd createQthVector(vector<double> ems, vector<double> ems_prime, vector<double> tau_th,
                            vector<double> q_th, vector<double> q_th_prime, vector<double> temperature){
        
    VectorXd constantVector;
    double stefanBoltzmann = 5.67e-8;
    
    //if (ems.size() != ems_prime.size() || ems.size() != ems_prime.size() || ems.size() != tau_th.size()
    //    || ems.size() != q_th.size() || ems.size() != q_th_prime.size() || ems.size() || temperature.size()) {
    //    throw std::invalid_argument("Size Error !");
    //}
    constantVector.resize(temperature.size() * 2);

    for(int i=0; i<temperature.size(); ++i){
        if(temperature.size()==2)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * stefanBoltzmann * pow(temperature[0], 4),
            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * stefanBoltzmann * pow(temperature[0], 4),
            (1 - ems_prime[1] - tau_th[1]) * q_th_prime[2] + ems_prime[1] * stefanBoltzmann * pow(temperature[1], 4),
            tau_th[1] * q_th_prime[2] + ems[1] * stefanBoltzmann * pow(temperature[1], 4);
        }
    
        else if (temperature.size()==3)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * stefanBoltzmann * pow(temperature[0], 4),
            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * stefanBoltzmann * pow(temperature[0], 4),
            ems_prime[1] * stefanBoltzmann * pow(temperature[1], 4),
            ems[1] * stefanBoltzmann * pow(temperature[1], 4),
            (1 - ems_prime[2] - tau_th[2]) * q_th_prime[3] + ems_prime[2] * stefanBoltzmann * pow(temperature[2], 4),
            tau_th[2] * q_th_prime[3] + ems[2] * stefanBoltzmann * pow(temperature[2], 4);
        }

        else if (temperature.size()==4)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * stefanBoltzmann * pow(temperature[0], 4),
            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * stefanBoltzmann * pow(temperature[0], 4),
            ems_prime[1] * stefanBoltzmann * pow(temperature[1], 4),
            ems[1] * stefanBoltzmann * pow(temperature[1], 4),
            ems_prime[2] * stefanBoltzmann * pow(temperature[2], 4),
            ems[2] * stefanBoltzmann * pow(temperature[2], 4),
            (1 - ems_prime[3] - tau_th[3]) * q_th_prime[4] + ems_prime[3] * stefanBoltzmann * pow(temperature[3], 4),
            tau_th[3] * q_th_prime[4] + ems[3] * stefanBoltzmann * pow(temperature[3], 4);
        }

        else if (temperature.size()==5)
        {
            constantVector << tau_th[0] * q_th[0] + ems_prime[0] * stefanBoltzmann * pow(temperature[0], 4),
            (1 - ems[0] - tau_th[0]) * q_th[0] + ems[1] * stefanBoltzmann * pow(temperature[0], 4),
            ems_prime[1] * stefanBoltzmann * pow(temperature[1], 4),
            ems[1] * stefanBoltzmann * pow(temperature[1], 4),
            ems_prime[2] * stefanBoltzmann * pow(temperature[2], 4),
            ems[2] * stefanBoltzmann * pow(temperature[2], 4),
            ems_prime[3] * stefanBoltzmann * pow(temperature[3], 4),
            ems[3] * stefanBoltzmann * pow(temperature[3], 4),
            (1 - ems_prime[4] - tau_th[4]) * q_th_prime[5] + ems_prime[4] * stefanBoltzmann * pow(temperature[4], 4),
            tau_th[4] * q_th_prime[5] + ems[4] * stefanBoltzmann * pow(temperature[4], 4);
        }    
    }
    
    return constantVector;
}

VectorXd calculVitesseStable(double T1, double T2, double T3, double hauteur, double largeur){
    //double As, At, Ab, Ah, Al, Ar;
    //double A1eq, A2eq;
    double Z1 = 0.44444444;
    double Z2 = 0.44444444;
    double g = 9.81;
    //double hg = 1.2081;
    double hc;

    double T0 = 293.15;
    double Ti = 25 + 273.15;
    double Tmp_2_3 = (T2 + T3) / 2.0;
    double Tgi = (Tmp_2_3 + Ti) / 2.0;
    double deltaT_2_3 = fabs(T3 - T2);
    
    double Tg;
    double Tsor;

    double hg;
    double rho;
    double muy;
    double c;
    double s = 0.100;

    double vitesse = 1.0;  // Initialisation de la vitesse
    double tol = 0.00001;  // Tolérance de convergence pour la vitesse
    double vitesse_precedente;
    bool continueLoop = true;
    int max_iter = 100;   // Pour éviter la boucle infinie
    int iter = 0;

    double DPB, DPHP, DPZ, DPT;
    double A, B, C;
    double Htp;

    /* 
    top = top / 1000.0;
    bottom = bottom / 1000.0;
    lateral = lateral / 1000.0;

    As = largeur * s;
    At = largeur * top;
    Ab = largeur * bottom;
    Ah = params.tau_th3 * largeur * hauteur;
    Al = hauteur * lateral;
    Ar = hauteur * lateral;

    A1eq = fmin(As, (Ab + (Al + Ar + Ah)/4));
    A2eq = fmin(As, (At + (Al + Ar + Ah)/4));

    double ZZ1  = pow(As/0.6f/A1eq-1, 2.0f);
    double ZZ2  = pow(As/0.6f/A2eq-1, 2.0f);

    (A1eq < 0.00001f) ? (A1eq = 0.00001f) : false;
    (A2eq < 0.00001f) ? (A2eq = 0.00001f) : false; 
    */

    rho = calculHg(T0, deltaT_2_3, s, {1, 0, 0, 0})[1];
    muy = calculHg(T0, deltaT_2_3, s, {1, 0, 0, 0})[2];
    c = calculHg(T0, deltaT_2_3, s, {1, 0, 0, 0})[3];
    hg = calculHg(T0, deltaT_2_3, s, {1, 0, 0, 0})[0];
    cout << "hg_1 = " << hg << endl;

    DPB = (rho * pow(vitesse, 2)) / 2;
    DPHP = (12 * muy * hauteur * vitesse) / pow(s, 2);
    DPZ = rho * vitesse * ((Z1 + Z2) / 2);

    while (continueLoop && iter < max_iter) {
        vitesse_precedente = vitesse;
        // Mise à jour des pertes de pression et des autres paramètres en fonction de la nouvelle vitesse
        DPT = 1.189 * 293 * 9.81 * hauteur * fabs(Tgi - Ti) / (Tgi * Ti);

        A = DPB + DPZ;
        B = DPHP;
        C = -DPT;
        vitesse = (-B+pow(B*B - 4*A*C,0.5))/2/A; 

        // Calcul de hc, Htp, Tsor, et Tg
        hc = 2 * hg + 4 * vitesse;
        Htp = rho * c * s * vitesse / (2 * hc);
        Tsor = Tmp_2_3 - (Tmp_2_3 - Ti) * exp(-hauteur / Htp);
        Tg = Tmp_2_3 - (Htp / hauteur) * (Tsor - Ti);
        hg = calculHg(Tg, deltaT_2_3, s, {1, 0, 0, 0})[0];
    
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

    // Affichage des résultats finaux
    //cout << "\nVitesse stabilisée : " << vitesse << " après " << iter << " itérations." << endl;
    return result;
}

VectorXd bilan(double Es, double hc_e, double hc_i, double Te, double Ti, vector<double> temperatureVector, vector<double> epsilon, 
                vector<double> epsilon_prime, vector<double> alphaE, vector<double> tau_th, double hauteur, double largeur){
    
    double A = 0.035;
    double n = 0.38;
    double T1 = temperatureVector[0];
    double T2 = temperatureVector[1];
    double T3 = temperatureVector[2];
    double Tmp_1_2 = (T1 + T2) / 2;
    double Tmp_2_3 = (T2 + T3) / 2;
    double deltaT_1_2 = fabs(T2 - T1);
    double deltaT_2_3 = fabs(T3 - T2);

    double epsilon1 = epsilon[0];
    double epsilon1_prime = epsilon_prime[0];
    double epsilon2 = epsilon[1];
    double epsilon2_prime = epsilon_prime[1];
    double epsilon3 = epsilon[2];
    double epsilon3_prime = epsilon_prime[2];

    double alphaE1 = alphaE[0];
    double alphaE2 = alphaE[1];
    double alphaE3 = alphaE[2];

    double tau_th1 = tau_th[0];
    double tau_th2 = tau_th[1];
    double tau_th3 = tau_th[2];
    
    //double Tm = (T1 + T2) / 2.0;
    double T0 = 293.15;
    double deltaT = 2.5; 
    double sigma = 5.67e-8;

    double s_1 = 0.018;

    double qth_0 = sigma * pow(Te,4);
    double qth_1, qth_2, qth_3;
    double qth_0_prime, qth_1_prime, qth_2_prime;
    double qth_3_prime = sigma * pow(Ti, 4);

    double bilan_1, bilan_2, bilan_3, df1_dT1, df1_dT2, df1_dT3, df2_dT1, df2_dT2, df2_dT3, df3_dT1, df3_dT2, df3_dT3;

    double hg_0 = hc_e;
    double hg_1;
    double hg_3 = hc_i;
   
    double hc;
    double hg;
    double Tg;
    double Tsor;
    double Htp;
    double vitesse;
    VectorXd v;

    MatrixXd thermalMatrix;
    MatrixXd constantVectorQth;
    VectorXd qthVectors;
    
    // Newton-Raphson iteration loop
    for (int iter = 0; iter < 100; ++iter) {
        // Recalculate thermophysical properties
        v = calculVitesseStable(T1, T2, T3, hauteur, largeur);
        //v =  calculVitesseStable(T1, T2, T3, hauteur, largeur, top, bottom, lateral, params)
        
        hc = v[0];
        Htp = v[1];
        Tsor = v[2];
        Tg = v[3];
        vitesse = v[4];
        hg = v[5];

        Tmp_1_2 = (T1+T2) / 2.0;
        Tmp_2_3 = (T2+T3) / 2.0;
        deltaT_1_2 = fabs(T2 - T1);
        deltaT_2_3 = fabs(T3 - T2);

        hg_1 = calculHg(Tmp_1_2, deltaT_1_2, s_1, {0.1, 0.9, 0, 0})[0];
        
        thermalMatrix = createQthMatrix({epsilon1, epsilon2, epsilon3}, {epsilon1_prime, epsilon2_prime, epsilon3_prime}, {tau_th1, tau_th2, tau_th3});
        constantVectorQth = createQthVector({epsilon1, epsilon2, epsilon3}, {epsilon1_prime, epsilon2_prime, epsilon3_prime}, {tau_th1, tau_th2, tau_th3},
                                            {qth_0, qth_1, qth_2, qth_3}, {qth_0_prime, qth_1_prime, qth_2_prime, qth_3_prime}, {T1, T2, T3});
        qthVectors = thermalMatrix.inverse() * constantVectorQth;

        //cout << "\n Qth Vectors : \n" << qthVectors << endl;
        qth_1 = qthVectors[0];
        qth_0_prime = qthVectors[1];
        qth_2 = qthVectors[2];
        qth_1_prime = qthVectors[3];
        qth_3 = qthVectors[4];
        qth_2_prime = qthVectors[5];

        double qth_a1 = epsilon1 * qth_0 + epsilon1_prime * qth_1_prime  - (epsilon1 + epsilon1_prime) * sigma * pow(T1, 4);
        double qth_a2 = epsilon2 * qth_1 + epsilon2_prime * qth_2_prime  - (epsilon2 + epsilon2_prime) * sigma * pow(T2, 4);
        double qth_a3 = epsilon3 * qth_2 + epsilon3_prime * qth_3_prime  - (epsilon3 + epsilon3_prime) * sigma * pow(T3, 4);
        //double qc_a1 = hg_0 * (T0 - T1) + hg_1 * (T2 - T1);
        //double qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
        double qc_a1 = hc_e * (Te - T1) + hg_1 * (Tmp_1_2 - T1);
        double qc_a2 = hg_1 * (Tmp_1_2 - T2) + hc * (Tg - T2);
        double qc_a3 = hc * (Tg - T3) + hc_i * (Ti - T3);

        // Recalculate bilan equations
        bilan_1 = alphaE1 * Es + epsilon1 * qth_0 + epsilon1_prime * qth_1_prime  - (epsilon1 + epsilon1_prime) * sigma * pow(T1, 4) + hc_e * (Te - T1) + hg_1 * (Tmp_1_2 - T1);
        bilan_2 = alphaE2 * Es + epsilon2 * qth_1 + epsilon2_prime * qth_2_prime  - (epsilon2 + epsilon2_prime) * sigma * pow(T2, 4) + hg_1 * (Tmp_1_2 - T2) + hc * (Tg - T2);
        bilan_3 = alphaE3 * Es + epsilon3 * qth_2 + epsilon3_prime * qth_3_prime  - (epsilon3 + epsilon3_prime) * sigma * pow(T3, 4) + hc * (Tg - T3) + hc_i * (Ti - T3);

        VectorXd equation(3);
        equation[0] = bilan_1;
        equation[1] = bilan_2;
        equation[2] = bilan_3;

        // Update Jacobian matrix
        df1_dT1 = -4.0 * (epsilon1 + epsilon1_prime) * sigma * pow(T1, 3) - hc_e - hg_1;
        df1_dT2 = hg_1;
        df1_dT3 = 0;
        df2_dT1 = hg_1;
        df2_dT2 = -4.0 * (epsilon2 + epsilon2_prime) * sigma * pow(T2, 3) - hg_1 - hc;
        df2_dT3 = 0;
        df3_dT1 = 0;
        df3_dT2 = 0;
        df3_dT3 = -4.0 * (epsilon3 + epsilon3_prime) * sigma * pow(T3, 3) - hc - hc_i;

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
        
        cout << "Iteration " << iter << ":  T1 = " << T1 << ", T2 = " << T2 << ", T3 = " << T3 <<endl;
        cout << "Iteration " << iter << ":  Tmp_1_2 = " << Tmp_1_2 << ", Tmp_2_3 = " << Tmp_2_3 <<endl;
        cout << "Iteration " << iter << ":  deltaT_1_2 = " << deltaT_1_2 << ", deltaT_2_3 = " << deltaT_2_3 <<endl;
        cout << "Iteration " << iter << ":  qc_a1 = " << qc_a1 << ", qc_a2 = " << qc_a2 << ", qc_a2 = " << qc_a3 << endl;
        cout << "Iteration " << iter << ": Bilan 1 = " << bilan_1 << ", Bilan 2 = " << bilan_2 << ", Bilan 3 = " << bilan_3 << endl; 
        
        // Check convergence
        if (fabs(bilan_1) < 0.1 && fabs(bilan_2) < 0.1) {
            break;
        } 
        // Update temperatures
        T1 += delta[0];
        T2 += delta[1];
        T3 += delta[2];
    }

    cout << "\n Qth Vectors : \n" << qthVectors << endl;

    VectorXd result(15);
    result[0] = T1;
    result[1] = T2;
    result[2] = T3;
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
    result[13] = qth_3;
    result[14] = qth_2_prime;

    return result;
}

struct ThermalParameters {
    double Te, Ti;
    double hc_e, hc_i;
    double alphaE1, alphaE2, alphaE3;
    double Es, transmission;
    double epsilon1, epsilon1_prime, tau_th1;
    double epsilon2, epsilon2_prime, tau_th2;
    double epsilon3, epsilon3_prime, tau_th3;
};

struct GlazingPropreties
{
    double hauteur, largeur;
    double top, bottom, lateral;
};


int main() {
    // Example usage of the createThermalMatrix function
    const double sigma = 5.67e-8;
    /* double Te = 298.15;
    double Ti = 298.15;
    double hc_e = 8.0;
    double hc_i = 2.5;
    double alphaE1 = 0.0979;
    double alphaE2 = 0.1488;
    double alphaE3 = 0.5941;    
    double Es = 500.0;
    double transmission = 0.006;

    double epsilon1 = 0.837;
    double epsilon1_prime = 0.837; 
    double tau_th1 = 0.0;   
    double epsilon2 = 0.837;
    double epsilon2_prime = 0.837;
    double tau_th2 = 0.0; 
    double epsilon3 = 0.837; 
    double epsilon3_prime = 0.837;
    double tau_th3 = 0.02;  */

    ThermalParameters params;

    // Get input from the user
    std::cout << "Enter Te: ";
    std::cin >> params.Te;
    std::cout << "Enter Ti: ";
    std::cin >> params.Ti;
    std::cout << "Enter hc_e: ";
    std::cin >> params.hc_e;
    std::cout << "Enter hc_i: ";
    std::cin >> params.hc_i;
    std::cout << "Enter Es: ";
    std::cin >> params.Es;
    std::cout << "Enter transmission: ";
    std::cin >> params.transmission;
    std::cout << "Enter alphaE1: ";
    std::cin >> params.alphaE1;
    std::cout << "Enter alphaE2: ";
    std::cin >> params.alphaE2;
    std::cout << "Enter alphaE3: ";
    std::cin >> params.alphaE3;
    std::cout << "Enter ems 1: ";
    std::cin >> params.epsilon1;
    std::cout << "Enter ems'1: ";
    std::cin >> params.epsilon1_prime;
    std::cout << "Enter ems 2: ";
    std::cin >> params.epsilon2;
    std::cout << "Enter ems'2: ";
    std::cin >> params.epsilon2_prime;
    std::cout << "Enter ems 3: ";
    std::cin >> params.epsilon3;
    std::cout << "Enter ems'3: ";
    std::cin >> params.epsilon3_prime;
    std::cout << "Enter tau_th1: ";
    std::cin >> params.tau_th1;
    std::cout << "Enter tau_th2: ";
    std::cin >> params.tau_th2;
    std::cout << "Enter tau_th3: ";
    std::cin >> params.tau_th3;

    double Te = params.Te;
    double Ti = params.Ti;
    double hc_e = params.hc_e;
    double hc_i = params.hc_i;
    double alphaE1 = params.alphaE1;
    double alphaE2 = params.alphaE2;
    double alphaE3 = params.alphaE3;    
    double Es = params.Es;
    double transmission = params.transmission;

    double epsilon1 = params.epsilon1;
    double epsilon1_prime = params.epsilon1_prime; 
    double tau_th1 = params.tau_th1;   
    double epsilon2 = params.epsilon2;
    double epsilon2_prime = params.epsilon2_prime;
    double tau_th2 = params.tau_th2; 
    double epsilon3 = params.epsilon3; 
    double epsilon3_prime = params.epsilon3_prime;
    double tau_th3 = params.tau_th3; 

    double qth_0 = sigma * pow(Te,4);
    double qth_0_prime;
    double qth_1;
    double qth_1_prime;
    double qth_2;
    double qth_2_prime;
    double qth_3;
    double qth_3_prime = sigma * pow(Ti,4);
    
    double qth_e, qth_i;
    double qth_a1, qth_a2, qth_a3;
    double qc_a1, qc_a2, qc_a3;

    double qc_e, qc_i;

    double hauteur = 1.4;
    double largeur = 2.4;

    cout << "\n------------------------Resultat Avec Flux------------------------------------" << endl;
    double T1_initial_avecFlux, T2_initial_avecFlux, T3_initial_avecFlux;

    MatrixXd temperatureMatrix = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, {epsilon1, epsilon2, epsilon3}, {epsilon1_prime, epsilon2_prime, epsilon3_prime});
    VectorXd constantVectorTemperature = createTemperatureVector(Es, hc_e, hc_i, Te, Ti, {alphaE1, alphaE2, alphaE3});
    VectorXd temperatureVector = temperatureMatrix.inverse() * constantVectorTemperature;

    T1_initial_avecFlux = temperatureVector[0];
    T2_initial_avecFlux = temperatureVector[1];
    T3_initial_avecFlux = temperatureVector[2];

    VectorXd bilanResult = bilan(Es, hc_e, hc_i, Te, Ti, {T1_initial_avecFlux, T2_initial_avecFlux}, {epsilon1, epsilon2, epsilon3}, 
                                {epsilon1_prime, epsilon2_prime, epsilon3_prime}, {alphaE1, alphaE2, alphaE3}, {tau_th1, tau_th2, tau_th3}, hauteur, largeur);    
    
    double T1_avecFlux = bilanResult[0];
    double T2_avecFlux = bilanResult[1];
    double T3_avecFlux = bilanResult[2];
    double hc_avecFlux = bilanResult[3];
    double Htp_avecFlux = bilanResult[4];
    double Tsor_avecFlux = bilanResult[5];
    double Tg_avecFlux = bilanResult[6];
    double vitesseAvecFlux = bilanResult[7];
    double hg_avecFlux = bilanResult[8];
    
    cout << "\nBilan Results (Avec Flux): " << endl;
    cout << "T1 = " << T1_avecFlux << " K" << endl;
    cout << "T2 = " << T2_avecFlux << " K" << endl;
    cout << "T3 = " << T3_avecFlux << "K" << endl;
    cout << "hc = " << hc_avecFlux << endl;
    cout << "Htp = " << Htp_avecFlux << endl;
    cout << "Tsor = " << Tsor_avecFlux << " K" << endl;
    cout << "Tg = " << Tg_avecFlux << " K" << endl;
    cout << "Vitesse = " << vitesseAvecFlux << endl;
    cout << "hg = " << hg_avecFlux << endl; 

    qth_1 = bilanResult[9];
    qth_0_prime = bilanResult[10];
    qth_2 = bilanResult[11];
    qth_1_prime = bilanResult[12];
    qth_3 =  bilanResult[13];
    qth_2_prime =  bilanResult[14];

    cout << "\nqth_0 (Es) = " << qth_0 << endl;
    cout << "qth_0_prime (Es) = " << qth_0_prime << endl;
    cout << "qth_1 (Es) = " << qth_1 << endl;
    cout << "qth_1_prime (Es) = " << qth_1_prime << endl;
    cout << "qth_2 (Es) = " << qth_2 << endl;
    cout << "qth_2_prime (Es) = " << qth_2_prime << endl;
    cout << "qth_3 (Es) = " << qth_3 << endl;
    cout << "qth_3_prime (Es) = " << qth_3_prime << endl;

    qth_e = qth_0_prime - qth_0;
    qth_i = qth_3 - qth_3_prime;
    
    cout << "\nqth_e (Es) = " << qth_e << endl;
    cout << "qth_i (Es) = " << qth_i << endl;

    qc_e = hc_e * (T1_avecFlux - Te);
    qc_i = hc_i * (Ti - T3_avecFlux);

    cout << "\nqc_e (Es) = " << qc_e << endl;
    cout << "qc_i (Es) = " << qc_i << endl;

    double q_v_Es = hc_avecFlux * (T2_avecFlux - Tg_avecFlux) + hc_avecFlux * (T3_avecFlux - Tg_avecFlux);

    cout << "\nq_v (Es) = " << q_v_Es << endl;

    cout << "\n------------------------Resultat Sans Flux------------------------------------" << endl;
    double T1_initial_sansFlux, T2_initial_sansFlux, T3_initial_sansFlux;

    MatrixXd temperatureMatrixSansFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, {epsilon1, epsilon2, epsilon3}, {epsilon1_prime, epsilon2_prime, epsilon3_prime});
    VectorXd constantVectorTemperatureSansFlux = createTemperatureVector(0, hc_e, hc_i, Te, Ti, {alphaE1, alphaE2, alphaE3});
    VectorXd temperatureVectorSansFlux = temperatureMatrixSansFlux.inverse() * constantVectorTemperatureSansFlux;

    T1_initial_sansFlux = temperatureVectorSansFlux[0];
    T2_initial_sansFlux = temperatureVectorSansFlux[1];
    T3_initial_sansFlux = temperatureVectorSansFlux[2];

    VectorXd bilanResultsansFlux = bilan(0, hc_e, hc_i, Te, Ti, {T1_initial_sansFlux, T2_initial_sansFlux}, {epsilon1, epsilon2, epsilon3}, 
                                {epsilon1_prime, epsilon2_prime, epsilon3_prime}, {alphaE1, alphaE2, alphaE3}, {tau_th1, tau_th2, tau_th3}, hauteur, largeur); 

    double T1_sansFlux = bilanResultsansFlux[0];
    double T2_sansFlux = bilanResultsansFlux[1];
    double T3_sansFlux = bilanResultsansFlux[2];
    double hc_sansFlux = bilanResultsansFlux[3];
    double Htp_sansFlux = bilanResultsansFlux[4];
    double Tsor_sansFlux = bilanResultsansFlux[5];
    double Tg_sansFlux = bilanResultsansFlux[6];
    double vitesseSansFlux = bilanResultsansFlux[7];
    double hg_sansFlux = bilanResultsansFlux[8];
    
    cout << "\nBilan Results (Sans Flux): " << endl;
    cout << "T1 = " << T1_sansFlux << " K" << endl;
    cout << "T2 = " << T2_sansFlux << " K" << endl;
    cout << "T3 = " << T3_sansFlux << " K" << endl;
    cout << "hc = " << hc_sansFlux << endl;
    cout << "Htp = " << Htp_sansFlux << endl;
    cout << "Tsor = " << Tsor_sansFlux << " K" << endl;
    cout << "Tg = " << Tg_sansFlux << " K" << endl;
    cout << "Vitesse = " << vitesseSansFlux << endl;
    cout << "hg = " << hg_sansFlux << endl; 

    double qth_1_sansFlux = bilanResultsansFlux[9];
    double qth_0_prime_sansFlux = bilanResultsansFlux[10];
    double qth_2_sansFlux = bilanResultsansFlux[11];
    double qth_1_prime_sansFlux = bilanResultsansFlux[12];
    double qth_3_sansFlux =  bilanResultsansFlux[13];
    double qth_2_prime_sansFlux =  bilanResultsansFlux[14];

    double qth_3_prime_sansFlux = sigma * pow(Ti,4);
    double qth_0_sansFlux = sigma * pow(Te,4);

    cout << "\nqth_0 (0) = " << qth_0_sansFlux << endl;
    cout << "qth_0_prime (0) = " << qth_0_prime_sansFlux << endl;
    cout << "qth_1 (0) = " << qth_1_sansFlux << endl;
    cout << "qth_1_prime (0) = " << qth_1_prime_sansFlux << endl;
    cout << "qth_2 (0) = " << qth_2_sansFlux << endl;
    cout << "qth_2_prime (0) = " << qth_2_prime_sansFlux << endl;
    cout << "qth_3 (0) = " << qth_3_sansFlux << endl;
    cout << "qth_3_prime (0) = " << qth_3_prime_sansFlux << endl;

    double qth_e_sansFlux = qth_0_prime_sansFlux - qth_0_prime_sansFlux;
    double qth_i_sansFlux = qth_3_sansFlux - qth_3_prime_sansFlux;
    
    cout << "\nqth_e (0) = " << qth_e_sansFlux << endl;
    cout << "qth_i (0) = " << qth_i_sansFlux << endl;

    double qc_e_sansFlux = hc_e * (T1_sansFlux - Te);
    double qc_i_sansFlux = hc_i * (Ti - T3_sansFlux);

    cout << "\nqc_e (0) = " << qc_e_sansFlux << endl;
    cout << "qc_i (0) = " << qc_i_sansFlux << endl;

    double q_v_sansFlux = hc_sansFlux * (T2_sansFlux - Tg_sansFlux) + hc_sansFlux * (T3_sansFlux - Tg_sansFlux);

    cout << "\nq_v (0) = " << q_v_sansFlux << endl;

 
    double g_th, g_c, g_v, g_tot;

    g_th = fabs(qth_i - qth_i_sansFlux) / Es;
    g_c = fabs(qc_i - qc_i_sansFlux) / Es;
    g_v = fabs(q_v_Es - q_v_sansFlux) / Es;

    //g_th = std::round(g_th * 1000.0) / 1000.0;
    //g_c = std::round(g_c * 1000.0) / 1000.0;
    //g_v = std::round(g_v * 1000.0) / 1000.0; 

    g_tot = g_th + g_c + g_v + transmission;
    //g_tot = std::round(g_tot * 1000.0) / 1000.0; 

    cout << "\n---------------------------- Final Result ----------------------------" << endl;
    cout << "g_th = " << g_th * 100.0 << " %" << endl;
    cout << "g_c = " << g_c * 100.0 << " %" << endl;
    cout << "g_v = " << g_v * 100.0 << " %" << endl;
    cout << "g_tot = " << g_tot * 100.0  << " %" <<endl; 
 
    return 0;
}
