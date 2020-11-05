#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <cmath>

#include "integral.h"

/// This is srtuct of columns;
/// You enter them in main() and take an ionosphere potential;
struct Column {
    double area;
    std::function<double(double)> conductivity;
    std::function<double(double)> current;
};

///At the moment it's useless
struct RToR {
    virtual ~RToR() {}
    virtual double operator()(double z) = 0;
};

///It calculates a geometric altitude from a geopotential altitude [km]
/// |DEF| Geometric height is an altitude above mean sea level.
/// |DEF| Geopotential is an adjustment to geometric height that
///       accounts for the variation of gravity with latitude and altitude.
double geom_from_gp(double H){
    double earth_radius = 6369.0; ///[km]
    return H * earth_radius / (H + earth_radius);
}
///altitude is from 0 km to 70 km
///temperature and pressure can be calculated
class StandardAtmosphere{
private:
    double T_0 = 288.15;
    double p_0 = 101'325;
    double R = 8.31432;
    double M = 28.9644;
    double g = 9.80655;
    std::vector <double> z = {0.0,
                             11.0,
                             20.0,
                             32.0,
                             47.0,
                             51.0,
                             70.0};
    std::vector <double> T = {T_0}, p = {p_0};
    std::vector <double> gamma = {0.0,
                                 -6.5,
                                 0.0,
                                 1.0,
                                 2.8,
                                 0.0,
                                 -2.8};
public:
    StandardAtmosphere(){
        T.resize(6);
        p.resize(6);
        std::cout << T[0] << std::endl;
        for(size_t n = 1; n < 7; ++n){
            T[n] = T[n-1] + gamma[n]*(z[n] - z[n-1]);
            std::cout << T[n] << std::endl;
            if (gamma[n] == 0.0){
                p[n] = p[n-1]*exp(-g*M*(z[n] - z[n-1])/(R*T[n-1]));
            } else {
                p[n] = p[n-1]*pow((1.0 + gamma[n]*(z[n] - z[n-1])/T[n-1]), -g*M/(gamma[n]*R));
            }
        }
    };
    double temperature(double zeta){
        size_t n = 1;
        for (size_t k = 1; k < 7; ++k) {
            if (zeta >= z[k-1] and z[k] >= zeta){n = k;}
        }
        return T[n-1] + gamma[n]*(zeta - z[n-1]);
    }
    double pressure(double alt){
        size_t n = 1;
        for (size_t k = 1; k < 7; ++k) {
            if (alt >= z[k-1] and z[k] >= alt){n = k;}
        }
        if (gamma[n] == 0.0){
            return p[n-1]*exp(-g*M*(alt - z[n-1])/(R*T[n-1]));
        } else {
            return p[n-1]*pow((1.0 + gamma[n]*(alt - z[n-1])/T[n-1]), -g*M/(gamma[n]*R));
        }
    }
};

/// explicit conductivity functions
class Conductivities {
protected:
    static constexpr double sigma_0 = 5.0e-14, H_0 = 6'000.0;
public:
    static double exp_conductivity(double z){
             return sigma_0*exp(z/H_0);
    };
    static double const_conductivity(double z){
        return sigma_0;
    }
};

/// explicit current functions
class Currents {
protected:
    static constexpr double j_0 = 1e-9;
public:
    static double step_current(double z){

        if (z>5'000.0 and z<10'000.0){
            return j_0;
        } else{
            return 0.0;
        }
    }
    static double zero_current(double z){
        return 0.0;
    }
};

struct Tank{
    std::vector<double> a;
    std::vector<double> b;
};


/// It is the central class, here you can find any parameters of the model
class GECModel {
protected:
    double H = 70'000;
    double V = 0;
    //std::vector<std::function<double(double)>> phi;
    std::vector<Tank> phi;
    size_t phi_points = 100;
    bool isVCalculated = false;
    bool isPhiCalculated = false;
    std::vector<Column> model;
    unsigned number_of_points = 100'001;
    /// This is a function that calculates an ionosphere potention
    double calc_ip() {
        std::vector<double> area_by_int(model.size());
        std::vector<double> int_of_cur_by_cond(model.size());
        for (unsigned i = 0; i < model.size(); ++i){
            area_by_int[i] = model[i].area / integrate_Simpson(
                        [this,i](double z) {return 1 / model[i].conductivity(z);},
                        0.0, H, number_of_points);
            int_of_cur_by_cond[i] = integrate_Simpson(
                        [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                        0.0, H, number_of_points);
        }
        double up = 0.0, down = 0.0;
        for (unsigned i = 0; i < model.size(); ++i){
            up += area_by_int[i] * int_of_cur_by_cond[i];
            down += area_by_int[i];
        }
        return up / down;
    }
    /// This is a function that calculates dependence of the potention on the varying altitude
    void calc_phi(){
        size_t new_number_of_points = 2.0 * unsigned(number_of_points / phi_points) + 1;
        phi.reserve(model.size());
        double I1 = 0.0, I2 = 0.0, C1, C2;
        for (unsigned i = 0; i < model.size(); ++i){
            std::vector<double> potential(phi_points), zeta(phi_points);
            C1 = integrate_Simpson(
                        [this,i](double z) {return 1 / model[i].conductivity(z);},
                        0.0, H, number_of_points);
            C2 = integrate_Simpson(
                        [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                        0.0, H, number_of_points);
            zeta[0] = 0.0;
            potential[0] = 0.0;
            for (size_t j = 1; j < phi_points; ++j){
                zeta[j] = j * H / phi_points;
                I1 = integrate_Simpson(
                            [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                            0.0, zeta[j], new_number_of_points);
                I2 = integrate_Simpson(
                            [this,i](double z) {return 1 / model[i].conductivity(z);},
                            0.0, zeta[j], new_number_of_points);
                potential[j] = I1 - (I2 / C1) * (C2 - getIP());
            }
            phi[i] = {zeta, potential};
        }
    }
public:
    GECModel(std::vector<Column> aModel): model(std::move(aModel)) {
    }
    GECModel() {
    }
    double getIP() {
        if (not isVCalculated) {
            V = calc_ip();
            isVCalculated = true;
        }
        return V;
    }
    double getPhi(size_t i, size_t j) {
        if (not isPhiCalculated) {
            calc_phi();
            isPhiCalculated = true;
        }
        return phi[i].b[j];
    }
    double getPhiX(size_t i, size_t j) {
        if (not isPhiCalculated) {
            calc_phi();
            isPhiCalculated = true;
        }
        return phi[i].a[j];
    }
    size_t getPhiPoints(){return phi_points;}
    double getH(){return H;}
};

///This is the considered model
class ConcreteGECModel: public GECModel, private Conductivities, private Currents {
public:
    ConcreteGECModel(): GECModel() {
        model.reserve(2);
        model.push_back({1.0, Conductivities::exp_conductivity, Currents::step_current});
        model.push_back({1.0, Conductivities::exp_conductivity, Currents::zero_current});
    }
};



class ConductivitySMZ15{
private:
    /// constants
    std::vector<double> A = {49.0 * M_PI / 180.0,
                             49.0 * M_PI / 180.0,
                             49.0 * M_PI / 180.0,
                             57.0 * M_PI / 180.0,
                             63.0 * M_PI / 180.0,
                             63.0 * M_PI / 180.0,
                             63.0 * M_PI / 180.0};
    std::vector<double> a = {8.0 * M_PI / 180.0,
                             8.0 * M_PI / 180.0,
                             8.0 * M_PI / 180.0,
                             10.0 * M_PI / 180.0,
                             10.0 * M_PI / 180.0,
                             9.0 * M_PI / 180.0};
    std::vector<double> H = {0.0, 5.0, 10.0, 16.0, 23.0, 37.0};                      ///[km]
    std::vector<double> B = {0.0, 0.0, 5.0, 21.0, 14.0, 0.0};                        ///[km]
    std::vector<double> b = {0.0, 0.0, 0.0, 14.0, 9.0, 0.0};                         ///[km]
    std::vector<double> C = {1.4E6, 15.0E6, 64.0E6, 93.0E6, 74.0E6, 45.0E6};        ///[m^(-3)*c^(-1)]
    std::vector<double> c = {0.1E6, 0.5E6, 2.0E6, 5.0E6, 3.0E6, 2.0E6};             ///[m^(-3)*c^(-1)]
    std::vector<double> D = {0.8E6, 10.0E6, 236.0E6, 402.0E6, 421.0E6, 450.0E6};    ///[m^(-3)*c^(-1)]
    std::vector<double> d = {0.1E6, 2.5E6, 83.0E6, 225.0E6, 236.0E6, 243.0E6};      ///[m^(-3)*c^(-1)]
    //...
public:
    double Lymbda(double lymbda, double lymbda_0){
        if(std::abs(lymbda) < lymbda_0){
            return std::abs(lymbda);
        } else {
            return lymbda_0;
        }
    }

    /// ion-pair production rate
    /// z is the altitude
    double q(double z, double lymbda, double xi, double p, double T){
        //...
    }



};




int main() {

    /*ConcreteGECModel m;
    std::ofstream fout("plots/potential_2_columns.txt");
    for (size_t i = 0; i < m.getPhiPoints(); ++i){
        fout << m.getPhiX(1,i) << "\t" << m.getPhi(0,i) << "\t" <<  m.getPhi(1,i) << std::endl;
    }
    fout << 70'000 << "\t" << m.getIP() << "\t" << m.getIP() << std::endl;
    fout.close();*/

    std::ofstream fout("plots/my_atm.txt");
    if (fout.is_open() == false){
        std::cout << "Impossible to find a file" << std::endl;
        return 1;
    }
    StandardAtmosphere atm;
    for (double z = 0.0; z <= 70.0; z += 0.1){
        fout << z << "\t" << atm.pressure(geom_from_gp(z)) << "\t" << atm.temperature(geom_from_gp(z)) << "\n";
    }
    fout.close();

    return 0;
}
