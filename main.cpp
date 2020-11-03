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
        for(size_t n = 1; n < 7; ++n){
            T[n] = T[n-1] + gamma[n]*(z[n] - z[n-1]);
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
    double pressure(double zeta){
        size_t n = 1;
        for (size_t k = 1; k < 7; ++k) {
            if (zeta >= z[k-1] and z[k] >= zeta){n = k;}
        }
        if (gamma[n] == 0.0){
            return p[n-1]*exp(-g*M*(zeta - z[n-1])/(R*T[n-1]));
        } else {
            return p[n-1]*pow((1.0 + gamma[n]*(zeta - z[n-1])/T[n-1]), -g*M/(gamma[n]*R));
        }
    }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
/// It takes data from file(-s) and approximates them;
/// The least square method return y(x) = a*x + b;
struct LinearLSM : public RToR {
    std::vector<double> x, y;
    LinearLSM(const char *filename) {                  /// Constructor makes a vector of points x_i and a vector of f(x_i);
        std::ifstream fin(filename);
        for (size_t i = 0; fin; ++i){
            fin >> x[i] >> y[i];
        }
        fin.close();
    }
    double operator()(double z) override {                      /// Calculation of the coefficients;
        double sumxx = 0.0;
        double sumxy = 0.0;
        double sumx = 0.0;
        double sumy = 0.0;
        for(unsigned i = 0; i < x.size(); ++i){
            sumxx += x[i] * x[i];
            sumxy += x[i] * y[i];
            sumx += x[i];
            sumy += y[i];
        }
        unsigned n = size(x);
        double a = (n * sumxy - sumx * sumy)/(n * sumxx - sumxx * sumxx);
        double b = (sumy - a * sumx)/n;
        return a * z + b;
    }
};

/// It takes data from file(-s) and approximates them;
/// The least square method return y(x) = exp(A * x + B);
struct ExponentialLSM : public RToR {
    std::vector<double> x, y;
    ExponentialLSM(const char *filename) {                     /// Constructor makes a vector of points x_i and a vector of f(x_i);
        std::ifstream fin(filename);
        for (size_t i = 0; fin; ++i){
            fin >> x[i] >> y[i];
            y[i] = log(y[i]);
        }
        fin.close();
    }
    double operator()(double z) override {                      /// Calculation of the coefficients like in the case of linear function;
        double sumxx = 0.0;
        double sumxy = 0.0;
        double sumx = 0.0;
        double sumy = 0.0;
        for(unsigned i = 0; i < x.size(); ++i){
            sumxx += x[i] * x[i];
            sumxy += x[i] * y[i];
            sumx += x[i];
            sumy += y[i];
        }
        unsigned n = size(x);
        double A = (n * sumxy - sumx * sumy)/                   /// Conversion of the coefficients for the considered case
                (n * sumxx - sumxx * sumxx);                    /// of exponentional function;
        double B = (sumy - A * sumx)/n;
        return exp(A * z + B);
    }
};

/// It takes data from file(-s) and approximates them by the cubic spline method;
/// It's important not to use this method with quick-change function
class CubicSpline : public RToR{
private:
    struct SplineSet{
        double a;
        double b;
        double c;
        double d;
        double x;
    };
    std::vector<SplineSet> coef;
    std::vector<SplineSet> spline(std::vector <double> &x, std::vector <double> &y){
        size_t n = x.size()-1;
        std::vector <double> a;
        a.insert(a.begin(), y.begin(), y.end());
        std::vector <double> b(n);
        std::vector <double> d(n);
        std::vector <double> h;

        for(size_t i = 0; i < n; ++i){
            h.push_back(x[i+1]-x[i]);
        }
        std::vector <double> alpha;
        alpha.push_back(0);
        for(size_t i = 1; i < n; ++i){
            alpha.push_back( 3*(a[i+1]-a[i])/h[i] - 3*(a[i]-a[i-1])/h[i-1]  );
        }
        std::vector <double> c(n+1);
        std::vector <double> l(n+1);
        std::vector <double> mu(n+1);
        std::vector <double> z(n+1);
        l[0] = 1;
        mu[0] = 0;
        z[0] = 0;
        for(size_t i = 1; i < n; ++i){
            l[i] = 2 *(x[i+1]-x[i-1])-h[i-1]*mu[i-1];
            mu[i] = h[i]/l[i];
            z[i] = (alpha[i]-h[i-1]*z[i-1])/l[i];
        }
        l[n] = 1;
        z[n] = 0;
        c[n] = 0;
        for(size_t j = n-1; j > -1; --j){
           // int j = i;
            c[j] = z [j] - mu[j] * c[j+1];
            b[j] = (a[j+1]-a[j])/h[j]-h[j]*(c[j+1]+2*c[j])/3;
            d[j] = (c[j+1]-c[j])/3/h[j];
        }
        std::vector<SplineSet> output_set(n);
        for(size_t i = 0; i < n; ++i){
            output_set[i].a = a[i];
            output_set[i].b = b[i];
            output_set[i].c = c[i];
            output_set[i].d = d[i];
            output_set[i].x = x[i];
        }
        return output_set;
    }
public:
    CubicSpline(std::vector <double> &x, std::vector <double> &y){
        coef = spline(x, y);
    }
    double operator()(double z){
        size_t i = 1;
        while(i < coef.size() and z > coef[i].x){
            ++i;
        }
        SplineSet set = coef[i-1];
        return set.a*(z-set.x)*(z-set.x)*(z-set.x) +
               set.b*(z-set.x)*(z-set.x) + set.c*(z-set.x) + set.d;
    }
};
/// It is the simplest current function which is identical zero
struct ZeroCurrent : public RToR {
    double operator()(double) override {
        return 0;
    }
};*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
            /*CubicSpline approx(zeta, potential);
            phi[i] = approx;*/
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

int main() {

    /*ConcreteGECModel m;
    std::ofstream fout("/home/mrkozelberg/Desktop/Work/Global_Electric_Circuit/plots/potential_2_columns.txt");
    for (size_t i = 0; i < m.getPhiPoints(); ++i){
        fout << m.getPhiX(1,i) << "\t" << m.getPhi(0,i) << "\t" <<  m.getPhi(1,i) << std::endl;
    }
    fout << 70'000 << "\t" << m.getIP() << "\t" << m.getIP() << std::endl;
    fout.close();*/

    std::ofstream fout("/home/mrkozelberg/Desktop/Work/Global_Electric_Circuit/plots/my_atm.txt");
    StandardAtmosphere atm;
    for (double z = 0.0; z <= 70.0; ++z){
        fout << z << "\t" << atm.pressure(z) << "\t" << atm.temperature(z) << "\n";
    }
    fout.close();

    return 0;
}
