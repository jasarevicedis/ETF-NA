
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>

const double eps = 0.0001;
typedef std::pair<double, double> Point;
typedef std::vector<std::pair<double, double>> PointSet;

bool isEqual(double a, double b)  {
    return std::fabs(a-b) < std::numeric_limits<double>::epsilon();
}

class AbstractInterpolator
{
protected:
    PointSet coordinates;
    mutable int interval;
    int Locate(double x) const;
public:
    AbstractInterpolator(const std::vector<std::pair<double, double> > &data);
    virtual double operator()(double x) const = 0;
};

AbstractInterpolator::AbstractInterpolator(const std::vector<std::pair<double, double> > &data){
    interval = 1;
    coordinates = data;
    std::sort(coordinates.begin(), coordinates.end(), [](Point t1, Point t2) {
        return t1.first < t2.first;
    });
    for(int i = 0; i < coordinates.size() - 1; i++) {
        if(isEqual(coordinates[i].first, coordinates[i + 1].first))
            throw std::domain_error("Invalid data set");
    }
}

int AbstractInterpolator::Locate(double x) const
    {
        if(x <= coordinates[0].first) 
            return 0;
        if(x > coordinates[coordinates.size() - 1].first) 
            return coordinates.size();
        if(coordinates[interval - 1].first <= x && x < coordinates[interval].first) 
            return interval;
        Point tmp(x,0);
        auto i = std::lower_bound(coordinates.begin(), coordinates.end(), tmp, [](Point t1, Point t2){
            return t1.first < t2.first;
        });
        interval = i - coordinates.begin();
        return interval;
    }
/*---------------------------------------------------------------------------------------------------------------*/

class LinearInterpolator: public AbstractInterpolator {
    double LinearIntCalculate(Point c1, Point c2, double x) const;
public:
    LinearInterpolator(const PointSet &data);
    double operator()(double x) const override;
};

//pomocna funkcija, racuna vrijednost y za odgovorajuci x na nekom intervalu
double LinearInterpolator::LinearIntCalculate(Point c1, Point c2, double x) const {
    return (( c2.first - x ) / ( c2.first - c1.first )) * c1.second + (( x - c1.first ) / ( c2.first - c1.first )) * c2.second;
}

LinearInterpolator::LinearInterpolator(const PointSet &data): AbstractInterpolator(data) {}

double LinearInterpolator::operator()(double x) const {
    int interval = Locate(x);
//linearna ekstrapolacija
    int n = coordinates.size();
    if(interval <= 1)
        return LinearIntCalculate(coordinates[0], coordinates[1], x);
    else if(interval >= n)
        return LinearIntCalculate(coordinates[n-2], coordinates[n-1], x);
//linearna interpolacija
    return LinearIntCalculate(coordinates[interval-1], coordinates[interval], x);
}


/*-----------------------------------------------------------------------------------------*/

class PolynomialInterpolator: public AbstractInterpolator {
    std::vector<double> polynomCoefficients;
    std::vector<double> dividedDifferences;//omogucavaju brze racunanje koeficijenata nakon dodavanja nove tacke
public:
    PolynomialInterpolator(const PointSet &data);
    double operator()(double x) const override;
    void AddPoint(const Point &p);
    std::vector<double> GetCoefficients() const;
};

PolynomialInterpolator::PolynomialInterpolator(const PointSet &data): AbstractInterpolator(data) {
    polynomCoefficients.resize(coordinates.size());
    int n = coordinates.size();//da ne racunamo u petlji

    for(int i = 0; i <n; i++)
        polynomCoefficients.push_back(coordinates[i].second);

    dividedDifferences.push_back(polynomCoefficients[0]);
    for(int j = 1; j < n; j++){
        for(int i = n; i > j; i--){
            polynomCoefficients[i-1] = ( polynomCoefficients[i-1] - polynomCoefficients[i-2]) / (coordinates[i-1].first - coordinates[i - j -1].first);
            if(i == n)
                dividedDifferences.push_back(polynomCoefficients[i-1]);
        }
    }
   
}
double PolynomialInterpolator::operator()(double x) const  {
    int n = coordinates.size();
//linearno vrijeme izvrsavanja, Hornerov algoritam
    double f = polynomCoefficients[n-1];
    for(int i= n-1; i>0; i--)
        f = f * (x - coordinates[i-1].first) + polynomCoefficients[i-1];
    return f;
}
void PolynomialInterpolator::AddPoint(const Point &p){
    int n= coordinates.size();
    for(int i=0;i < n; i++)
        if(isEqual(coordinates[i].first, p.first))
            throw std::domain_error("Invalid point");
    
    coordinates.push_back(p);
    n = coordinates.size();
    polynomCoefficients.push_back(p.second);
    for(int j=1;j<=n-1;j++){
        if(j==1){
            polynomCoefficients[n-1] = (polynomCoefficients[n-1] - coordinates[n-2].second) / ( coordinates[n-1].first - coordinates[n-j-1].first);
        }
        else {
            double x = dividedDifferences[j-1];
            dividedDifferences[j-1] = polynomCoefficients[n-1];
            polynomCoefficients[n-1] = (polynomCoefficients[n-1] - x) / (coordinates[n-1].first - coordinates[n-j-1].first);
        }
    }
    dividedDifferences.push_back(polynomCoefficients[n-1]);
}

std::vector<double> PolynomialInterpolator::GetCoefficients() const {
    int n = coordinates.size();
    std::vector<double> p(n); //koeficijenti polinoma
    std::vector<double> w(n+1); //master polinom 
    std::vector<double> v(n+1); //vektor u koji kopiramo koeficijente master polinoma 

//algoritam zasnovan na upotrebi master polinoma, str. 12
    w[0] = 1;
    double a;
    for(int i = 1; i <= n; i++){
        w[i] = w[i-1];
        for(int j = i-1; j>= 1; j--)
            w[j] = w[j-1] - coordinates[i-1].first * w[j];
        w[0] *=  -coordinates[i-1].first;
    }
    for(int i =1; i <= n; i++){
        a = 1; //nema potrebe za cuvanjem u vektoru jer nam poslije ne treba
        for(int j =1; j <= n; j++){
            if(j != i)
                a *= coordinates[i-1].first - coordinates[j-1].first; 
        }
        a = coordinates[i-1].second / a;
        for(int j=0; j<= n; j++)//copy
            v[j] = w[j];
        for(int j = n-1; j>= 0; j--){
            v[j] += coordinates[i-1].first * v[j+1];
            p[j] += a * v[j+1];
        }
    }
    return p;
}

/*------------------------------------------------------------------------------------------*/

class PiecewisePolynomialInterpolator: public AbstractInterpolator {
    int k;//red interpolacije
public:
    PiecewisePolynomialInterpolator(const PointSet &data, int order);
    double operator()(double x) const override;
};

PiecewisePolynomialInterpolator::PiecewisePolynomialInterpolator(const PointSet &data, int order): AbstractInterpolator(data) {
    if(order < 1 || order > coordinates.size())
        throw std::domain_error("Invalid order");
    k=order;
}
double PiecewisePolynomialInterpolator::operator()(double x) const  {
    int i =Locate(x);//interval
    int n = coordinates.size();
    int g1=0,g2=0; 
//odredjivanje granica za j, p6 str 17
    if(k%2 == 0){
        g1 = i - k/2;
        g2 = i + k/2;
        if(g1 < 1){
            g1=1;
            g2 = k+1;
        }
        else if(g2 > n){
            g2 = n;
            g1 = g2 - k;
        }
    }
    else {
        g1 = i - (k-1)/2;
        g2 = i + (k+1)/2;
        if(g1<1){
            g1 = 1;
            g2 = k+1;
        }
        else if(g2 > n){
            g2 =n;
            g1 = g2 - k;
        }
    }
//interpoliranje, lagrangova formula s promjenjenim granicama za i i j petlju, p6 str 5
    double s(0);
    for(i = g1; i<= g2; i++){//interval nam ne treba vise
        double p = coordinates[i-1].second;
        for(int j = g1; j<= g2; j++)
            if(j!=i)
                p *= (x - coordinates[j-1].first) / (coordinates[i-1].first - coordinates[j-1].first);                
        s += p;
    }
    return s;
}

/*------------------------------------------------------------------------------------------*/
class SplineInterpolator: public AbstractInterpolator {
    std::vector<double> r;//dovoljni su nam r koeficijenti
public:
    SplineInterpolator(const PointSet &data);
    double operator()(double x) const override;
};

SplineInterpolator::SplineInterpolator(const PointSet &data): AbstractInterpolator(data) {
    int n = coordinates.size();
    r.resize(n);
    std::vector<double> a(n-1);//ne treba nam nakon nalaska r
//konstrukcija spline-a p7 str 5
    for(int i = 1; i < n-1; i++){
        a[i] = 2 * (coordinates[i+1].first - coordinates[i-1].first);
        r[i] = 3 * ( ( coordinates[i+1].second - coordinates[i].second ) / (coordinates[i+1].first - coordinates[i].first ) - ( coordinates[i].second - coordinates[i-1].second) / (coordinates[i].first - coordinates[i-1].first) );
    }
    for(int i = 1; i < n-2; i++){
        double mi = ( coordinates[i+1].first - coordinates[i].first) / a[i];
        a[i+1] -= mi * ( coordinates[i+1].first - coordinates[i].first);
        r[i+1] -= mi * r[i];
    }
    r[n-2] /= a[n-2];
    for(int i = n-3; i>0; i--){
        r[i] = ( r[i] - ( coordinates[i+1].first - coordinates[i].first) * r[i+1] ) / a[i];
    }
}
double SplineInterpolator::operator()(double x) const  {
    int i = Locate(x); // interval
    int n = coordinates.size();
    if(i == 0)
        i = 1;
    else if(i == n)
        i --;
//racunanje potrebnih parametara
    double delta_x = coordinates[i].first - coordinates[i-1].first;
    double t = x - coordinates[i-1].first;
    double s = ( r[i] - r[i-1]) / (3 * delta_x);
    double q = ( coordinates[i].second - coordinates[i-1].second) / delta_x - delta_x * (r[i] + 2*r[i-1] ) / 3.;
//interpoliranje
    return coordinates[i-1].second + t * ( q + t *(r[i-1] + t*s));
}

/*-------------------------------------------------------------------------------------------------------------------*/

class BarycentricInterpolator: public AbstractInterpolator {
    std::vector<double> weights;
    double barycentricOrder;
public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order);
    double operator()(double x) const override;
    std::vector<double> GetWeights() const;
};

BarycentricInterpolator::BarycentricInterpolator(const PointSet &data, int order): AbstractInterpolator(data) {
    int n = coordinates.size();
    if(order > n || order < 0)
        throw std::domain_error("Invalid order");
//za potrebe koeficijenata    
    barycentricOrder = order;
    weights.resize(coordinates.size());
    int m,l;
    double p;
//odredjivanje tezisnih koeficijenata
    for(int i = 0; i < n; i++){
        weights[i] = 0;
        if( i - barycentricOrder > 1)
            m = i - barycentricOrder;
        else
            m=1;

        if(i < n - barycentricOrder)
            l=i;
        else
            l = n - barycentricOrder - 1;

        for(int k = m-1; k < l + 1; k++){
            p=1;
            for(int j = k; j < k + barycentricOrder; j++)
                if(j != i)
                    p /= (coordinates[i].first - coordinates[j].first);
            if(k%2 == 1)
                p *= (-1);

        }
        weights[i] += p;
    }
}

double BarycentricInterpolator::operator()(double x) const {
    //(slučaj kad je x = x_i za neko 𝑖 = 1,2,…,𝑛 tretira se kao poseban slučaj
    double p(0);
    double q(0);
    double u;
    int n = coordinates.size();

    for(int i = 0; i < n; i++){
        if(isEqual(x,coordinates[i].first))
            return coordinates[i].second;
        u = weights[i] / ( x - coordinates[i].first);
        p +=  u * coordinates[i].first;
        q += u;
    }
    return p / q;
}
std::vector<double> BarycentricInterpolator::GetWeights() const {
    return weights;
}

/*------------------------------------------------------------------------------------------------*/
class TrigonometricInterpolator: public AbstractInterpolator {
public:
    TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data);
    double operator()(double x) const override;
};

TrigonometricInterpolator::TrigonometricInterpolator(const PointSet &data): AbstractInterpolator(data) {
    if( !isEqual(coordinates[0].second,coordinates[coordinates.size() - 1].second ))
        throw std::domain_error("Function is not periodic");
    
}
double TrigonometricInterpolator::operator()(double x) const {
    return 3;
}





int main(){
    try {
        LinearInterpolator li({{1,1}, {1,1}});
    }catch(...){
        std::cout << "Greska: Iste tacke u data setu!"<< std::endl;
    }

    try {
        PolynomialInterpolator pi({{1,1}, {2,2}});
        pi.AddPoint({1,3});
    }catch(...){
        std::cout << "Greska: Iste tacke u data setu!"<< std::endl;
    }

    try {
        PiecewisePolynomialInterpolator pi({{1,1}, {2,2}}, -1);
    }catch(...){
        std::cout << "Greska: Pogresan red!"<< std::endl;//slicno za baricentrincu
    }


    LinearInterpolator li({{1,1}, {2,2}, {3,3}});
    std::cout << li(1.5)<< " " << li(2)<< " " << li(2.5) << std::endl;

    PointSet s1({{1,1}, { 2, 4}, {3, 9}, {4,16}, {5,25}, {6,36}, {7,49}});
    PointSet s2({{1,1}, {3,5}, {5,2}, {1.5, 8}});

    PolynomialInterpolator pi(s1);
    std::cout << pi(1.5)<< " " << pi(2)<< " " << pi(2.5) << std::endl;
    pi = PolynomialInterpolator(s2);
    std::cout << pi(1.5)<< " " << pi(2)<< " " << pi(2.5) << std::endl;


    PiecewisePolynomialInterpolator pc(s1, 3);
    std::cout << pc(1.5)<< " " << pc(2)<< " " << pc(2.5) << std::endl;
    pc = PiecewisePolynomialInterpolator(s2,3);
    std::cout << pc(1.5)<< " " << pc(2)<< " " << pc(2.5) << std::endl;

    SplineInterpolator sl(s1);
    std::cout << sl(1.5)<< " " << sl(2)<< " " << sl(2.5) << std::endl;
    sl = SplineInterpolator(s2);
    std::cout << sl(1.5)<< " " << sl(2)<< " " << sl(2.5) << std::endl;

    BarycentricInterpolator bi(s1, 2);
    std::cout << bi(1.5)<< " " << bi(2)<< " " << bi(2.5) << std::endl;
    bi = BarycentricInterpolator(s2,2);
    std::cout << bi(1.5)<< " " << bi(2)<< " " << bi(2.5) << std::endl;


    return 0;
}
