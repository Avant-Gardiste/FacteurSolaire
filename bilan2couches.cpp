#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace Eigen;
using namespace std;

MatrixXd CreateTemperatureMatrix(double Te, double Ti, double hc_e, double hc_i, double AlphaE1, double AlphaE2){
    double sigma = 5.67e-8;
    double T0 = 293.15;
    double deltaT = 5.0;
    double A = 0.035;
    double n = 0.38;

    double he, hi;
    double hr_1, hr_2;
    double hg_1;

    double rho = 1.189;
    double muy = 0.00001811;
    double lambda = 0.02576;
    double c_1 = 1008.0;
    double s_1 = 0.050;
    double Gr, Pr, Nu;

    double Tm = (T0 + Te) / 2;

    double rho_1 = 1.189 - 0.0044*(Tm - 293.15); // 20°C = 293K (Annexe F - 52022-3)
    double muy_1 = 1.811e-5 + 0.005e-5*(Tm - 293.15);
    double lambda_1 = 2.576e-2 + 0.008e-2*(Tm - 293.15);

    Gr = (9.81*pow(s_1, 3)*deltaT*pow(rho_1, 2)) / (Tm*pow(muy_1,2));
    Pr = (muy_1*c_1) / lambda_1;
    Nu = fmax(A*pow(Gr*Pr, n), 1);

    hg_1 = Nu * (lambda_1/s_1);

    //cout << "Nu_1 = " << Nu << endl;
    //cout << "hg_1 = " << hg_1 << endl;

    he = hc_e + 4.0*sigma*pow(Te,3);
    hi = hc_i + 4.0*sigma*pow(Ti, 3);
    //cout << "\nhe = " << he << endl;
    //cout << "hi = " << hi << endl;

    hr_1 = 4.0*0.9*sigma*pow(T0,3);
    hr_2 = 4.0*0.9*sigma*pow(T0,3);
    //cout << "\nhr_1 = " << hr_1 << endl;
    //cout << "hr_2 = " << hr_2 << endl;

    MatrixXd temperatureMatrix = MatrixXd::Zero(2, 2);

    temperatureMatrix(0,0) = -(he+hr_1+hg_1);
    temperatureMatrix(1,0) = hr_2;
    temperatureMatrix(0,1) = hr_1;
    temperatureMatrix(1,1) = -(hr_2+hi+hg_1);

    return temperatureMatrix;
}

VectorXd createTemperatureVector(double alphaE1, double alphaE2, double Es, double hc_e, double hc_i, double Te, double Ti) {
    
    VectorXd temperatureVect(2);
    double he, hi;
    double sigma = 5.67e-8;
    double hg_1;
    double Tg1;

    double rho = 1.189;
    double muy = 0.00001811;
    double lambda = 0.02576;
    double c_1 = 1008.0;
    double T0 = 293.15; // 20°C
    double deltaT = 5.0;
    double s_1 = 0.050;

    double A = 0.035;
    double n = 0.38;
    double Gr, Pr, Nu;

    double Tm = (T0 + Te) / 2.0;

    double rho_1 = 1.189 - 0.0044*(Tm - 293.15); // 20°C = 293K (Annexe F - 52022-3)
    double muy_1 = 1.811e-5 + 0.005e-5*(Tm - 293.15);
    double lambda_1 = 2.576e-2 + 0.008e-2*(Tm - 293.15);

    Gr = (9.81*pow(s_1, 3)*deltaT*pow(rho_1, 2)) / (Tm*pow(muy_1,2));
    Pr = (muy_1*c_1) / lambda_1;
    Nu =  fmax(A*pow(Gr*Pr, n), 1);

    hg_1 = Nu * (lambda_1/s_1);
    Tg1 = T0;

    he = hc_e + 4.0*sigma*pow(Te,3);
    hi = hc_i + 4.0*sigma*pow(Ti, 3);
    temperatureVect << -alphaE1*Es - he*Te - hg_1*Tg1,
                    -alphaE2*Es - hi*Ti - hg_1*Tg1;
    
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


VectorXd calculVitesseStable(double T1, double T2, double hauteur, double largeur){
    double Z1 = 3.133;
    double Z2 = 5.392;
    double g = 9.81;
    //double hg = 1.2081;
    double hc;

    double Tmp = (T1 + T2) / 2.0;
    double T0 = 293.15;
    double Ti = 25.0 + 273.15;
    double Tgi = (Tmp + Ti) / 2.0;
    double deltaT = 5.0;
    
    double Tg;
    double Tsor;

    double rho = 1.189 - 0.0044 * (T0 - 293.15);
    double muy = 1.811e-5 + 0.005e-5 * (T0 - 293.15);
    double lambda = 2.576e-2 + 0.008e-2*(T0 - 293.15);
    double c = 1008.0;
    double s = 0.050;
    double Pr = (muy*c) / lambda;
    double Gr = (9.81*pow(s, 3)*deltaT*pow(rho, 2)) / (Tmp*pow(muy,2));
    double Nu =  fmax(0.035*pow(Gr*Pr, 0.38), 1);
    double hg = Nu * lambda / s;


    double vitesse = 1.0;  // Initialisation de la vitesse
    double tol = 0.0000001;     // Tolérance de convergence pour la vitesse
    double vitesse_precedente;
    bool continueLoop = true;
    int max_iter = 100;     // Pour éviter la boucle infinie
    int iter = 0;

    double DPB, DPHP, DPZ, DPT;
    double A, B, C;
    double Htp;

    DPB = (rho * pow(vitesse, 2)) / 2;
    DPHP = (12 * muy * hauteur * vitesse) / pow(s, 2);
    DPZ = rho * vitesse * ((Z1 + Z2) / 2);
    
    while (continueLoop && iter < max_iter) {
        vitesse_precedente = vitesse;

        // Mise à jour des pertes de pression et des autres paramètres en fonction de la nouvelle vitesse
        //DPB = (rho * pow(vitesse, 2)) / 2;
        //DPHP = (12 * muy * hauteur * vitesse) / pow(s, 2);
        //DPZ = rho * vitesse * ((Z1 + Z2) / 2); 
         
        DPT = 1.189 * 293 * 9.81 * hauteur * fabs(Tgi - Ti) / (Tgi * Ti);

        A = DPB + DPZ;
        B = DPHP;
        C = -DPT;

        /* cout << "\nA = " << A << endl;
        cout << "B = " << B << endl;
        cout << "C = " << C << endl; */

        // Recalcul de la vitesse
        //A = 0.5*rho*(1+Z1+Z2);
        //B = 12.0 * muy * hauteur / pow(s,2);
        //C = -DPT;

        vitesse = (-B+pow(B*B - 4*A*C,0.5))/2/A; 

        // Calcul de hc, Htp, Tsor, et Tg
        hc = 2.0 * hg + 4.0 * vitesse;
        Htp = rho * c * s * vitesse / (2.0 * hc);
        Tsor = Tmp - (Tmp - Ti) * exp(-hauteur / Htp);
        Tg = Tmp - (Htp / hauteur) * (Tsor - Ti);
        
        rho = 1.189 - 0.0044 * (Tg - 293.15);
        muy = 1.811e-5 + 0.005e-5 * (Tg - 293.15);
        lambda = 2.576e-2 + 0.008e-2*(Tg - 293.15);

        c = 1008.0;
        s = 0.050;
        Gr = (9.81*pow(s, 3)*deltaT*pow(rho, 2)) / (Tmp*pow(muy,2));
        Pr = (muy*c) / lambda;
        Nu =  fmax(0.035*pow(Gr*Pr, 0.38), 1);
        hg = Nu * lambda / s;


        // Mise à jour de Tgi pour le prochain calcul de DPT
        Tgi = Tg;

        /* cout << "\nhc = " << hc << endl;
        cout << "Htp = " << Htp << endl;
        cout << "Tsor = " << Tsor - 273.15 << endl;
        cout << "Tg = " << Tg - 273.15 << endl;
        cout << "Vitesse = " << vitesse << endl; */
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

VectorXd bilan(vector<double> temperatureVector, double Es, double epsilon1, double epsilon1_prime, double epsilon2, double epsilon2_prime,
        double alphaE1, double alphaE2, double hc_e, double hc_i, double Te, double Ti, double tau_th1, double tau_th2, double hauteur, double largeur){
    // Mettre à jour hg_1

    double A = 0.035;
    double n = 0.38;
    double T1 = temperatureVector[0];
    double T2 = temperatureVector[1];
    double Tmp = (T1 + T2) / 2.0;
    double T0 = 293.15;
    double deltaT = 5.0;
    double T3 = Ti;
    double sigma = 5.67e-8;

    double Gr;
    double Pr;
    double Nu; 
    
    // Calculer hg_1 en fonction des nouvelles valeurs de T1, T2 et Tm
    double hg_2 = hc_i;
    //double hg_1;
    //double hg_1 = Nu * lambda_1 / s_1;

    double qth_0 = sigma * pow(Te,4);
    double qth_1, qth_2;
    double qth_0_prime, qth_1_prime;
    double qth_2_prime = sigma * pow(Ti, 4);

    double bilan_1, bilan_2, df1_dT1, df1_dT2, df2_dT1, df2_dT2;


    double Tm = (T0 + Te) / 2.0;

    double rho_1 = 1.189 - 0.0044*(Tm - 293.15); // 20°C = 293K (Annexe F - 52022-3)
    double muy_1 = 1.811e-5 + 0.005e-5*(Tm - 293.15);
    double lambda_1 = 2.576e-2 + 0.008e-2*(Tm - 293.15);
    double c_1 = 1008.0;
    double s_1 = 0.050;

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

        v = calculVitesseStable(T1, T2, hauteur, largeur);
        hc = v[0];
        Htp = v[1];
        Tsor = v[2];
        Tg = v[3];
        vitesse = v[4];
        hg = v[5];
        Tmp = (T1 + T2) / 2.0;

        rho_1 = 1.189 - 0.0044*(Tmp - 293.15); // 20°C = 293K (Annexe F - 52022-3)
        muy_1 = 1.811e-5 + 0.005e-5*(Tmp - 293.15);
        lambda_1 = 2.576e-2 + 0.008e-2*(Tmp - 293.15);
        Gr = (9.81*pow(s_1, 3)*deltaT*pow(rho_1, 2)) / (Tmp*pow(muy_1,2));
        Pr = (muy_1*c_1) / lambda_1;
        Nu =  fmax(A*pow(Gr*Pr, n), 1);
        hg_1 = Nu * (lambda_1 / s_1);
    
        // Recalculate thermophysical properties
        thermalMatrix = createQthMatrix({epsilon1, epsilon2}, {epsilon1_prime, epsilon2_prime}, {tau_th1, tau_th2});
        constantVectorQth = createQthVector({epsilon1, epsilon2}, {epsilon1_prime, epsilon2_prime}, {tau_th1, tau_th2},
                                                    {qth_0, qth_1, qth_2}, {qth_0_prime, qth_1_prime, qth_2_prime}, {T1, T2});
        qthVectors = thermalMatrix.inverse() * constantVectorQth;

        qth_1 = qthVectors[0];
        qth_0_prime = qthVectors[1]; 
        qth_2 = qthVectors[2]; 
        qth_1_prime = qthVectors[3];

        double qth_a1 = epsilon1 * qth_0 + epsilon1_prime * qth_1_prime  - (epsilon1 + epsilon1_prime) * sigma * pow(T1, 4);
        double qth_a2 = epsilon2 * qth_1 + epsilon2_prime * qth_2_prime  - (epsilon2 + epsilon2_prime) * sigma * pow(T2, 4);
        //double qc_a1 = hg_0 * (T0 - T1) + hg_1 * (T2 - T1);
        //double qc_a2 = hg_1 * (T1 - T2) + hg_2 * (T3 - T2);
        double qc_a1 = hc_e * (Te - T1) + hc * (Tg - T1);
        double qc_a2 = hc * (Tg - T2) + hc_i * (Ti - T2);

        // Recalculate bilan equations
        bilan_1 = alphaE1 * Es + epsilon1 * qth_0 + epsilon1_prime * qth_1_prime  - (epsilon1 + epsilon1_prime) * sigma * pow(T1, 4) + hc_e * (Te - T1) + hc * (Tg - T1);
        bilan_2 = alphaE2 * Es + epsilon2 * qth_1 + epsilon2_prime * qth_2_prime  - (epsilon2 + epsilon2_prime) * sigma * pow(T2, 4) + hc * (Tg - T2) + hc_i * (Ti - T2);

        VectorXd equation(2);
        equation[0] = bilan_1;
        equation[1] = bilan_2;

        // Update Jacobian matrix
        df1_dT1 = -4.0 * (epsilon1 + epsilon1_prime) * sigma * pow(T1, 3) - hc_e - hc;
        df1_dT2 = 0;
        df2_dT1 = 0;
        df2_dT2 = -4.0 * (epsilon2 + epsilon2_prime) * sigma * pow(T2, 3) - hc - hc_i;

        MatrixXd jacobian(2, 2);
        jacobian(0, 0) = df1_dT1;
        jacobian(0, 1) = df1_dT2;
        jacobian(1, 0) = df2_dT1;
        jacobian(1, 1) = df2_dT2;

        VectorXd delta = jacobian.colPivHouseholderQr().solve(-equation);  // small adjustment to see effect 
        
        cout << "Iteration " << iter << ":  T1 = " << T1 << ", T2 = " << T2 << endl;
        cout << "Iteration " << iter << ":  qc_a1 = " << qc_a1 << ", qc_a2 = " << qc_a2 << endl;
        cout << "Iteration " << iter << ": Bilan 1 = " << bilan_1 << ", Bilan 2 = " << bilan_2 << endl;
        
        // Check convergence
        if (fabs(bilan_1) < 0.00001 && fabs(bilan_2) < 0.00001) {
            break;
        } 
        // Update temperatures
        T1 += delta[0];
        T2 += delta[1];
        //Tm += (T1 + T2) / 2.0;  // Update mean temperature
        //hc += hc_1;
        //Tgi += Tg; 

    }

    cout << "\n Qth Vectors : \n" << qthVectors << endl;
    
    VectorXd result(13);
    result[0] = T1;
    result[1] = T2;
    result[2] = hc;
    result[3] = Htp;
    result[4] = Tsor;
    result[5] = Tg;
    result[6] = vitesse;
    result[7] = hg;
    result[8] = qth_1;
    result[9] = qth_0_prime;
    result[10] = qth_2;
    result[11] = qth_1_prime;
    return result;
}


int main() {
    // Example usage of the createThermalMatrix function
    double sigma = 5.67e-8;
    double Te = 278.15;
    double Ti = 293.15;
    double hc_e = 18.0;
    double hc_i = 3.6;
    double alphaE1 = 0.0839;
    double alphaE2 = 0.1055;
    double Es = 500.0;
    double T0 = 293.15;
    double deltaT = 5.0;

    double A = 0.035;
    double n = 0.38;
    double Gr, Pr, Nu;
    double hg_0, hg_1, hg_2;

    double epsilon1 = 0.837657;
    double epsilon1_prime = 0.837657; 
    double tau_th1 = 0.0;   
    double epsilon2 = 0.837657; 
    double epsilon2_prime = 0.837657;
    double tau_th2 = 0.03;   
    
    double qth_0 = sigma * pow(Te,4);
    double qth_0_prime;
    double qth_1;
    double qth_1_prime;
    double qth_2;
    double qth_2_prime = sigma * pow(Ti,4);

    double qth_e, qth_i;
    double qth_a1, qth_a2;
    double qc_a1, qc_a2;
    
    double hauteur = 2;
    double largeur = 1;

    // condition initial pour T0 et deltaT;
    MatrixXd temperatureMatrix = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, alphaE1, alphaE2);
    VectorXd constantVectorTemperature = createTemperatureVector(alphaE1, alphaE2, Es, hc_e, hc_i, Te, Ti);
    VectorXd temperatureVector = temperatureMatrix.inverse() * constantVectorTemperature;

    double T1_initial = temperatureVector[0];
    double T2_initial = temperatureVector[1];
    double T3 = Ti;

    cout << "\nMatrice temperature  : \n" << temperatureVector<< endl;
    cout << "\nT1 Initiale  : " << T1_initial << endl;
    cout << "T2 Initiale  : " << T2_initial << endl;

    VectorXd bilanResult = bilan({T1_initial, T2_initial}, Es, epsilon1, epsilon1_prime, epsilon2, epsilon2_prime, alphaE1, alphaE2, hc_e, hc_i, Te, Ti, tau_th1, tau_th2, hauteur, largeur);    

    double T1_avecFlux = bilanResult[0];
    double T2_avecFlux = bilanResult[1];
    double hc_avecFlux = bilanResult[2];
    double Htp_avecFlux = bilanResult[3];
    double Tsor_avecFlux = bilanResult[4];
    double Tg_avecFlux = bilanResult[5];
    double vitesseAvecFlux = bilanResult[6];
    double hg_avecFlux = bilanResult[7];

    cout << "\nBilan Results (Avec Flux): " << endl;
    cout << "T1 = " << T1_avecFlux << endl;
    cout << "T2 = " << T2_avecFlux << endl;
    cout << "hc = " << hc_avecFlux << endl;
    cout << "Htp = " << Htp_avecFlux << endl;
    cout << "Tsor = " << Tsor_avecFlux << " C" << endl;
    cout << "Tg = " << Tg_avecFlux << " C" << endl;
    cout << "Vitesse = " << vitesseAvecFlux << endl;
    cout << "hg = " << hg_avecFlux << endl; 

    qth_1 = bilanResult[8];
    qth_0_prime = bilanResult[9];
    qth_2 = bilanResult[10];
    qth_1_prime = bilanResult[11];

    cout << "\nqth_0 (Es) = " << qth_0 << endl;
    cout << "qth_0_prime (Es) = " << qth_0_prime << endl;
    cout << "qth_1 (Es) = " << qth_1 << endl;
    cout << "qth_1_prime (Es) = " << qth_1_prime << endl;
    cout << "qth_2 (Es) = " << qth_2 << endl;
    cout << "qth_2_prime (Es) = " << qth_2_prime << endl;

    qth_e = qth_0_prime - qth_0;
    qth_i = qth_2 - qth_2_prime;
    
    cout << "\nqth_e (Es) = " << qth_e << endl;
    cout << "qth_i (Es) = " << qth_i << endl;

    //qc_e = hc_e * (T1_avecFlux - Te);
    ///qc_i = hc_i * (Ti - T3_avecFlux);
    double qc_e_Es = hc_e * (T1_avecFlux - Te);
    double qc_i_Es = hc_i * (Ti - T2_avecFlux);

    cout << "\nqc_e (Es) = " << qc_e_Es << endl;
    cout << "qc_i (Es) = " << qc_i_Es << endl;

    /* cout <<"Tg avec flux : " << Tg_avecFlux << endl;
    cout << "T2 avec flux : " << T2_avecFlux << endl;
    cout << "hc avec flux :" << hc_avecFlux << endl; */
    
    double q_v_Es = hc_avecFlux * (Tg_avecFlux - T1_avecFlux) + hc_avecFlux * (Tg_avecFlux - T2_avecFlux);

    cout << "\nq_v (Es) = " << q_v_Es << endl;   
    

    cout << "\n------------------------Resultat Sans Flux------------------------------------" << endl;
    
    MatrixXd temperatureMatrixsansFlux = CreateTemperatureMatrix(Te, Ti, hc_e, hc_i, alphaE1, alphaE2);
    VectorXd constantVectorTemperaturesansFlux = createTemperatureVector(alphaE1, alphaE2, 0.0, hc_e, hc_i, Te, Ti);
    VectorXd temperatureVectorsansFlux = temperatureMatrixsansFlux.inverse() * constantVectorTemperaturesansFlux;

    double T1_initial_sansflux = temperatureVectorsansFlux[0];
    double T2_initial_sansflux = temperatureVectorsansFlux[1];
    //double T3 = Ti;

    cout << "\nMatrice temperature  : \n" << temperatureVector<< endl;
    cout << "\nT1 Initiale  : " << T1_initial << endl;
    cout << "T2 Initiale  : " << T2_initial << endl;

    VectorXd bilanResultSansFlux = bilan({T1_initial, T2_initial}, 0.0, epsilon1, epsilon1_prime, epsilon2, epsilon2_prime, alphaE1, alphaE2, hc_e, hc_i, Te, Ti, tau_th1, tau_th2, hauteur, largeur);    

    double T1_sansFlux = bilanResultSansFlux[0];
    double T2_sansFlux = bilanResultSansFlux[1];
    double hc_sansFlux = bilanResultSansFlux[2];
    double Htp_sansFlux = bilanResultSansFlux[3];
    double Tsor_sansFlux = bilanResultSansFlux[4];
    double Tg_sansFlux = bilanResultSansFlux[5];
    double vitesseSansFlux = bilanResultSansFlux[6];
    double hg_sansFlux = bilanResultSansFlux[7];

    cout << "\nBilan Results (Avec Flux): " << endl;
    cout << "T1 = " << T1_sansFlux << endl;
    cout << "T2 = " << T2_sansFlux << endl;
    cout << "hc = " << hc_sansFlux << endl;
    cout << "Htp = " << Htp_sansFlux << endl;
    cout << "Tsor = " << Tsor_sansFlux << " C" << endl;
    cout << "Tg = " << Tg_sansFlux << " C" << endl;
    cout << "Vitesse = " << vitesseSansFlux << endl;
    cout << "hg = " << hg_sansFlux << endl; 

    double qth_1_sansFlux = bilanResultSansFlux[8];
    double qth_0_prime_sansFlux = bilanResultSansFlux[9];
    double qth_2_sansFlux = bilanResultSansFlux[10];
    double qth_1_prime_sansFlux = bilanResultSansFlux[11];

    /* cout << "\nqth_0 (Es) = " << qth_0 << endl;
    cout << "qth_0_prime (Es) = " << qth_0_prime << endl;
    cout << "qth_1 (Es) = " << qth_1 << endl;
    cout << "qth_1_prime (Es) = " << qth_1_prime << endl;
    cout << "qth_2 (Es) = " << qth_2 << endl;
    cout << "qth_2_prime (Es) = " << qth_2_prime << endl; */

    double qth_e_sansFlux = qth_0_prime_sansFlux - qth_0;
    double qth_i_sansFlux = qth_2_sansFlux - qth_2_prime;
    
    cout << "\nqth_e (0) = " << qth_e_sansFlux << endl;
    cout << "qth_i (0) = " << qth_i_sansFlux << endl;

    //qc_e = hc_e * (T1_avecFlux - Te);
    ///qc_i = hc_i * (Ti - T3_avecFlux);
    double qc_e_sansFlux = hc_e * (T1_sansFlux - Te);
    double qc_i_sansFlux = hc_i * (Ti - T2_sansFlux);

    cout << "\nqc_e (0) = " << qc_e_sansFlux << endl;
    cout << "qc_i (0) = " << qc_i_sansFlux << endl;

    /* cout <<"Tg avec flux : " << Tg_avecFlux << endl;
    cout << "T2 avec flux : " << T2_avecFlux << endl;
    cout << "hc avec flux :" << hc_avecFlux << endl; */
    
    double q_v_sansFlux = hc_sansFlux * (Tg_sansFlux - T1_sansFlux) + hc_sansFlux * (Tg_sansFlux - T2_sansFlux);

    cout << "\nq_v (0) = " << q_v_sansFlux << endl;   
    
 
    double g_th, g_c, g_v, g_tot;

    g_th = fabs(qth_i - qth_i_sansFlux) / Es;
    g_c = fabs(qc_i_Es - qc_i_sansFlux) / Es;
    g_v = fabs(q_v_Es - q_v_sansFlux) / Es;

    //g_th = std::round(g_th * 1000.0) / 1000.0;
    //g_c = std::round(g_c * 1000.0) / 1000.0;
    //g_v = std::round(g_v * 1000.0) / 1000.0; 

    g_tot = g_th + g_c + g_v + 0.166;
    //g_tot = std::round(g_tot * 1000.0) / 1000.0; 

    cout << "\n---------------------------- Final Result ----------------------------" << endl;
    cout << "g_th = " << g_th * 100.0 << " %" << endl;
    cout << "g_c = " << g_c * 100.0 << " %" << endl;
    cout << "g_v = " << g_v * 100.0 << " %" << endl;
    cout << "g_tot = " << g_tot * 100.0  << " %" <<endl;  

    return 0;
}
