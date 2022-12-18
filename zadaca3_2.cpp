#include <iostream>
#include <cmath>
#include <limits>
#include <vector>

template <typename FunType>
std::pair<double, bool> Limit(FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20){
    if(eps <= 0 || nmax <3 || nmax > 30)
        throw std::domain_error("Invalid parameters");
    if(std::fabs(h) < eps)
        h = 0.001 * std::max(1. , x0);

    double y_old = std::numeric_limits<double>::infinity();
    std::vector<double> y(nmax);
    double p;

//Nevillov algoritam za Richardsonovu ekstrapolaciju ka granici
    for(int i = 0;i < nmax; i++){
        y[i] = f(x0 + h);
        p=2;
        for(int k = i - 1; k >= 0; k--){
            y[k] = ( p * y[k+1] - y[k] ) / ( p - 1 );
            p*=2;
        }
        if(std::fabs(y[0] - y_old) < eps )
            return std::make_pair(y[0], true );
        y_old = y[0];
        h/=2.;
    }
    return std::make_pair(y[0], false);
}

int main(){
    std::pair<double,bool> lim;

//1
    lim = Limit( [](double x){ return 1 / x; }, 0); //infinty
    std::cout << "Value: " << lim.first  << std::endl;
//2
    lim = Limit( [](double x){ return sin(x) / x; }, 0); //1
    std::cout << "Value: " << lim.first  << std::endl;
//3
    lim = Limit( [](double x){ return ( 1 - cos(x)) / x; }, 0); //0
    std::cout << "Value: " << lim.first  << std::endl;
//4
    lim = Limit( [](double x){ return std::pow( 1 + 1/x, x); }, std::numeric_limits<double>::infinity()); //1
    std::cout << "Value: " << lim.first  << std::endl;
//5
    lim = Limit( [](double x){ return std::exp(x); }, 0); //1
    std::cout << "Value: " << lim.first  << std::endl;
//Ocekivane vrijednosti se dobijaju
    return 0;
}
