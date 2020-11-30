#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <string>
#include <cmath>

#include "integral.h"
#include "StdAtm.h"
#include "sigma.h"

/// This is srtuct of columns;
/// You enter them in main() and take an ionosphere potential;
struct Column {
    double area;
    std::function<double(double)> conductivity;
    std::function<double(double)> current;
};

/// explicit conductivity functions
class Conductivities: private SMZ15, private TZ06, private ZT07{
protected:
    static constexpr double sigma_0 = 5.0e-14;
    static constexpr double H_0 = 6.0;                      /// [km]
    static constexpr double e_0 = 1.602176634e-19;          /// [C]
public:
    ///z in kilometres
    static double exp_conductivity(double z){
             return sigma_0*exp(z/H_0);
    };
    static double const_conductivity(double z){
        return sigma_0;
    }
    double conductivity(double z, double lymbda, double xi){
        return e_0 * mu_sum(z) * sqrt(q(z, lymbda, xi) / alpha(z));
    }
};

/// explicit current functions
class Currents{
protected:
    static constexpr double j_0 = 1e-9;
public:
    static double step_current(double z){

        if(z>=5.0 and z<10.0){
            return j_0;
        } else{
            return 0.0;
        }
    }
    static double zero_current(double z){
        return 0.0;
    }
    static double simple_geo_current(double lat, double z){
        if(std::abs(lat) <= 10){
            return step_current(z);
        } else return zero_current(z);
    }
};

/// It is the central class, here you can find any parameters of the model
class GECModel{
protected:
    static constexpr double H = 70.0;                   ///in km
    static constexpr unsigned int_points = 10'001;
    static constexpr unsigned int_pot_points = 9;
    double V = 0.0;
    bool isVCalculated = false;
    std::vector<Column> model;
    static constexpr double pot_step = 1.0;       ///in km

    /// This is a function that calculates an ionosphere potention
    double calc_ip() {
        std::vector<double> area_by_int(model.size());
        std::vector<double> int_of_cur_by_cond(model.size());
        for(unsigned i = 0; i < model.size(); ++i){
            area_by_int[i] = model[i].area / integrate_Simpson(
                        [this,i](double z) {return 1 / model[i].conductivity(z);},
                        0.0, H, int_points);
            int_of_cur_by_cond[i] = integrate_Simpson(
                        [this,i](double z) {return model[i].current(z) / model[i].conductivity(z);},
                        0.0, H, int_points);
        }
        double up = 0.0, down = 0.0;
        for(unsigned i = 0; i < model.size(); ++i){
            up += area_by_int[i] * int_of_cur_by_cond[i];
            down += area_by_int[i];
        }
        return up / down;
    }

    /// This is a function that calculates the potention on the i * pot_step
    std::vector<double> calc_pot(unsigned column_num){
        unsigned N = ceil(H / pot_step);
        std::vector<double> vec(N);
        double I1 = 0.0, I2 = 0.0, C1, C2;
        C1 = integrate_Simpson([this, column_num](double z) {return 1.0 / model[column_num].conductivity(z);}, 0.0, H, int_points);
        C2 = integrate_Simpson([this, column_num](double z) {return model[column_num].current(z) / model[column_num].conductivity(z);}, 0.0, H, int_points);
        std::vector<double> h(N);
        for(unsigned n = 0; n < N; ++n){
            h[n] = n * pot_step;
        }
        for(unsigned n = 1; n < N; ++n){
            I1 += integrate_Simpson([this, column_num](double z) {return model[column_num].current(z) / model[column_num].conductivity(z);},
                                    h[n-1], h[n], int_pot_points);
            I2 += integrate_Simpson([this, column_num](double z) {return 1.0 / model[column_num].conductivity(z);},
                                    h[n-1], h[n], int_pot_points);
            vec[n] = I1 - I2 * (C2 - getIP()) / C1;
        }
        return vec;
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
    void getPot(std::string filename, unsigned column_num){
        unsigned N = ceil(H / pot_step);
        std::vector<double> vec;
        vec = calc_pot(column_num);
        std::ofstream fout(filename);
        if(fout.is_open() == false){
            std::cout << "Impossible to find a file" << std::endl;
            exit(-1);
        }
        for(unsigned n = 1; n < N; ++n){
            fout << n * pot_step << "\t" << vec[n] << std::endl;
        }
        fout << 70 << "\t" << getIP() << std::endl;
        fout.close();
    }
};

class SimpliestGECModel: public GECModel, private Conductivities, private Currents {
public:
    SimpliestGECModel(): GECModel() {
        ///This is the simplest parametrization
        model.reserve(2);
        model.push_back({1.0, [this](double z){return Conductivities::conductivity(z, 1.0, 1.0);}, Currents::step_current});
        model.push_back({1.0, [this](double z){return Conductivities::conductivity(z, 1.0, 1.0);}, Currents::zero_current});
    }
};

class Geomodel: public GECModel, private Conductivities, private Currents {
private:
    double earth_radius2 = 40408473.9788; //km or m?
public:
    Geomodel(double delta_lat, double delta_lon): GECModel(){
        unsigned N = std::ceil(180.0 / delta_lat);
        unsigned M = std::ceil(360.0 / delta_lon);
        model.reserve(N*M);
        ///может быть проблема на концах
        for(unsigned n = 0; n < N; ++n){
            double lat_n = -90.0 + delta_lat * (0.5 + n);
            for(unsigned m = 0; m < M; ++m){
                model.push_back({delta_lon * M_PI / 180.0 * earth_radius2
                                 * 2.0 * sin(M_PI / 180 * lat_n) * sin(M_PI / 180 * delta_lat),
                                [this, lat_n](double z){return Conductivities::conductivity(z, lat_n, 1.0);},
                                [lat_n](double z){return Currents::simple_geo_current(lat_n, z);}});
            }
        }
    }
};

int main(){

    /// The simplest parametrization (2 columns, considered conductivity, step and zero currents)
    SimpliestGECModel m;
    m.getPot("plots/potential_2_columns.txt", 0);

    return 0;
}
