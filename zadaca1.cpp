#include <exception>
#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>

class Vector {
    std::vector<double> elements;
public:
    explicit Vector(int n){
        if(n<=0)
            throw std::range_error("Bad dimension");
        elements.resize(n);
    }
    Vector(std::initializer_list<double> l){
        if(l.size()==0)
            throw std::range_error("Bad dimension");
        elements.resize(l.size());
        auto it=l.begin();
        for(int i=0;i<l.size();i++){
            elements.at(i)=*it;
            it++;
        }
    }
    int NElems() const {
        return elements.size();
    }
    double &operator[](int i){
        return elements[i];
    }
    double operator[](int i) const {
        return elements[i];
    }
    double &operator()(int i) {
        if(i > elements.size())
            throw std::range_error("Invalid index");
        return elements[i-1];
    }
    double operator()(int i) const{
        if(i > elements.size())
            throw std::range_error("Invalid index");
        return elements[i-1];
    }
    double Norm() const {
        double sum=0;
        for(int i=0;i<elements.size();i++)
            sum+=std::pow(elements[i],2);
        return std::sqrt(sum);
    }
    friend double VectorNorm(const Vector &v){
        double sum=0;
        for(int i=0;i<v.elements.size();i++)
            sum+=v.elements.at(i)*v.elements.at(i);
        return std::sqrt(sum);
    }
    double GetEpsilon() const {
        return 10 * this->Norm() * std::numeric_limits<double>::epsilon();
    }
    void Print(char separator = '\n', double eps = -1) const {
        if(eps < 0)
            eps = GetEpsilon();
        for(auto it = elements.begin(); it !=elements.end();it++){
            double value = *it;
            if(fabs(value) < eps)
                value = 0;
            std::cout << value;
            if(it!=elements.end()-1)
                std::cout << separator;
            if(it==elements.end() && separator=='\n')
                std::cout << std::endl;
        }
    }
    friend void PrintVector(const Vector &v, char separator = '\n',
    double eps = -1){
        if(eps < 0)
            eps = v.GetEpsilon();
        for(auto it = v.elements.begin(); it !=v.elements.end();it++){
            double value = *it;
            if(fabs(value) < eps)
                value = 0;
            std::cout << value;
            if(it!=v.elements.end()-1)
                std::cout << separator;
            if(it==v.elements.end() && separator=='\n')
                std::cout << std::endl;
        }
    }
    friend Vector operator +(const Vector &v1, const Vector &v2){
        Vector v(v1);
        return v+=v2;
    }
    Vector &operator +=(const Vector &v){
        if(elements.size() != v.elements.size())
            throw std::domain_error("Incompatible formats");  
        for(int i=0;i<v.elements.size();i++)
            elements[i] += v.elements[i];
        return *this;
    }
    friend Vector operator -(const Vector &v1, const Vector &v2){
        Vector v(v1);
        return v-=v2;
    }
    Vector &operator -=(const Vector &v){
        if(elements.size() != v.elements.size())
            throw std::domain_error("Incompatible formats"); 
        for(int i=0;i<v.elements.size();i++)
            elements[i] -= v.elements[i];
        return *this;
    }
    friend Vector operator *(double s, const Vector &v){
        Vector tmp(v);
        return tmp*=s;
    }
    friend Vector operator *(const Vector &v, double s){
        Vector proizvod{};
        proizvod.elements.resize(v.elements.size());
        for(int i=0;i< v.elements.size();i++)
            proizvod.elements[i]= s * v.elements[i];
        return proizvod;
    }
    Vector &operator *=(double s){
        for(int i=0;i<elements.size();i++)
            elements[i]= elements[i]*s;
        return *this;
    }
    friend double operator *(const Vector &v1, const Vector &v2){
        if(v1.elements.size() != v2.elements.size())
            throw std::domain_error("Incompatible formats"); 
        double proizvod(0);
        for(int i=0;i< v1.elements.size();i++)
            proizvod+=v1.elements[i]*v2.elements[i];
        return proizvod;
    }
    friend Vector operator /(const Vector &v, double s){
       auto tmp(v);
       return tmp/=s;
    }
    Vector &operator /=(double s){
        if(s==0)
            throw std::domain_error("Division by zero");
        for(int i=0;i<elements.size();i++)
            elements[i]/=s;
        return *this;
    }
};

class Matrix {
    std::vector<std::vector<double>> elements;
    std::vector<std::vector<double>> KopirajMatrix(){
        std::vector<std::vector<double>> kopija;
        kopija.resize(elements.size());
        for(auto i=0;i< elements.size();i++)
            kopija.at(i).resize(elements.at(i).size());
        for(int i=0;i<elements.size();i++)
            for(int j=0;j<elements.at(i).size();j++)
                kopija[i][j]= elements[i][j];
        return kopija;
    }
public:
    Matrix(int m, int n){
        if(m<=0 || n<=0)
            throw std::range_error("Bad dimension");
        this->elements = std::vector<std::vector<double>>(m,std::vector<double>(n));
    }
    /*Vector becoming 1-row matrix*/
    Matrix(const Vector &v): elements(v.NElems(),std::vector<double>(1)) {
        for(int i=0;i< elements.size();i++)
            elements[i][0]=v[i];
    }
    Matrix(std::initializer_list<std::vector<double>> l){
        if(l.size()==0 || (*l.begin()).size() == 0 )
            throw std::range_error("Bad dimension");
        int vel = (*l.begin()).size();
        for(auto it = l.begin();it != l.end();it++)
            if(vel != (*it).size())
                throw std::logic_error("Bad matrix");
        for(auto it = l.begin();it != l.end();it++)
            elements.push_back(*it);
    }
    int NRows() const {
        return elements.size();
    }
    int NCols() const {
        return elements.at(0).size();
    }
    double *operator[](int i){
        return &elements[i][0];
    }
    const double *operator[](int i) const{
        return &elements[i][0];
    }
    double &operator()(int i, int j){
        if(i>elements.size() || j > elements[0].size())
            throw std::range_error("Invalid index");
        return elements[i-1][j-1];
    }
    double operator()(int i, int j) const {
        if(i>elements.size() || j > elements[0].size())
            throw std::range_error("Invalid index");
        return elements[i-1][j-1];
    }
    double Norm() const {
        double sum(0);
        for(int i=0;i<elements.size();i++)
            for(int j=0;j<elements.at(i).size();j++)
                sum+= elements[i][j];
        return std::sqrt(sum);
    }
    friend double MatrixNorm(const Matrix &m){
        double sum(0);
        for(int i=0;i<m.elements.size();i++)
            for(int j=0;j<m.elements.at(i).size();j++)
                sum+= m.elements[i][j];
        return std::sqrt(sum);
    }
    double GetEpsilon() const {
        return 10 * this->Norm() * std::numeric_limits<double>::epsilon();
    }
    void Print(int width = 10, double eps = -1) const{
        for(int i=0;i< NRows();i++){
            for(int j=0;j<NCols();j++){
                auto tmp = elements[i][j];
                if(std::fabs(tmp) < eps)
                    tmp=0;
                std::cout << std::setw(width) << tmp;
            }
            std::cout << std::endl;
        }
    }
    friend void PrintMatrix(const Matrix &m, int width = 10, double eps = -1){
        for(int i=0;i< m.NRows();i++){
            for(int j=0;j<m.NCols();j++){
                auto tmp = m.elements[i][j];
                if(std::fabs(tmp) < eps)
                    tmp=0;
                std::cout << std::setw(width) << tmp;
            }
            std::cout << std::endl;
        }
    }
    friend Matrix operator +(const Matrix &m1, const Matrix &m2){
        Matrix tmp(m1);
        return tmp+=m2;
    }
    Matrix &operator +=(const Matrix &m){
        if(NRows() != m.NRows() || NCols() != m.NCols())
            throw std::domain_error("Incompatible formats");
        for(int i=0;i<elements.size();i++)
            for(int j=0;j<elements[0].size();j++)
                elements[i][j] += m[i][j];
        return *this;
    }
    friend Matrix operator -(const Matrix &m1, const Matrix &m2){
        Matrix tmp(m1);
        return tmp-=m2;
    }
    Matrix &operator -=(const Matrix &m){
        if(NRows() != m.NRows() || NCols() != m.NCols())
            throw std::domain_error("Incompatible formats");
        for(int i=0;i<elements.size();i++)
            for(int j=0;j<elements[0].size();j++)
                elements[i][j] -= m[i][j];
        return *this;
    }
    
    friend Matrix operator *(double s, const Matrix &m){
        Matrix proizvod{int(m.elements.size()),int(m.elements[0].size())};
        for(int i=0;i<m.elements.size();i++)
            for(int j=0;j<m.elements[0].size();j++)
                proizvod[i][j] = m[i][j] * s;
        return proizvod;
    }
    friend Matrix operator *(const Matrix &m, double s){
        Matrix proizvod{int(m.elements.size()),int(m.elements[0].size())};
        for(int i=0;i<m.elements.size();i++)
            for(int j=0;j<m.elements[0].size();j++)
                proizvod[i][j] = m[i][j] * s;
        return proizvod;
    }
    Matrix &operator *=(double s){
        
        for(int i=0;i<elements.size();i++)
            for(int j=0;j<elements[0].size();j++)
                elements[i][j] *=  s;
        return *this;
    }
    friend Matrix operator *(const Matrix &m1, const Matrix &m2){
        Matrix tmp(m1);
        return tmp*=m2;
    }
    Matrix &operator *=(const Matrix &m){
        if(NCols() != m.NRows())
            throw std::domain_error("Incompatible formats");
        Matrix tmp(NRows(),m.NCols());
        for(int i=0;i<NRows();i++)
            for(int j=0;j<m.NCols();j++)
                for(int k=0;k<NCols();k++)
                    tmp[i][j]+=this->elements[i][k]*m[k][j];
        this->elements = tmp.elements;
        return *this;
    }
    friend Vector operator *(const Matrix &m, const Vector &v){
        auto tmp(m);
        tmp*=v;
        Vector tmp2(tmp.NRows());
        for(int i=0;i<tmp.NRows();i++)
            tmp2[i]=tmp[i][0];
        return tmp2;
    }
    friend Matrix Transpose(const Matrix &m){
        /*
        Matrix tmp{int(m.elements[0].size()), int(m.elements.size())};
        for(int i=0;i<tmp.elements.size();i++)
            for(int j=0;j<tmp.elements[0].size();j++)
                tmp.elements[i][j]=m.elements[j][i];
        return tmp;*/
        Matrix tmp(m);
        tmp.Transpose();
        return tmp;
    }
    void Transpose(){
        if(NCols()==NRows())
            for(int i=0;i<NRows()-1;i++)
                for(int j=i+1;j<NCols();j++)
                    std::swap(elements[i][j],elements[j][i]);
        else {
            std::vector<std::vector<double>>  tmp(NCols(), std::vector<double>(NRows()));
            for(int i=0;i<NRows();i++)
                for(int j=0;j<NCols();j++)
                    tmp[j][i]=elements[i][j];
            elements=tmp;
        }
    }
};
int main(){
    /*
    Vector v1(5);
    Vector v2{1,2,3,4,5};

    try{
        Vector v3{-100};
    }catch(std::exception& e){
        std::cout << e.what() << std::endl;
    }
    try{
        Vector v3{};
    }catch(std::exception& e){
        std::cout << e.what() << std::endl;
    }
    Vector v3(v2);
    v1.Print('.');
    std::cout << std::endl;
    v2.Print('.');
    std::cout << std::endl;
    for(int i=0;i<v2.NElems();i++)
        std::cout << v1[i] << " ";
    std::cout << std::endl;
    for(int i=1;i<=v2.NElems();i++)
        std::cout << v2(i) << " ";
    std::cout << std::endl;
    

    Matrix mat1(2,3);
    Matrix mat2{{1,2,3},{4,5,6}};
    Matrix mat3(Vector({7,8,9}));
    //konstruktori
    try{
        Matrix{{},{}};
    }catch(std::exception &e){
        std::cout << e.what() << std::endl;
    }
    try{
        Matrix{-100,-100};
    }catch(std::exception &e){
        std::cout << e.what() << std::endl;
    }
    try{
        Matrix{{1,1,1,1,1,1},{}};
    }catch(std::exception &e){
        std::cout << e.what() << std::endl;
    }
    try{
        Matrix{{1,1,1,1,1,1,1,1,1,1},{2}};
    }catch(std::exception &e){
        std::cout << e.what() << std::endl;
    }
    //PRINT
    mat1.Print(7);
    std::cout<<std::endl;
    mat2.Print(7);
    std::cout<<std::endl;
    mat3.Print(7);
    std::cout<<std::endl;
    (mat1+=mat2).Print(7);
    std::cout<<std::endl;
    (mat1 *=10).Print(7);
    std::cout<<std::endl;
    (mat1-=mat2).Print(7);
    std::cout<<std::endl;
    mat2.Transpose();
    mat3.Transpose();
    (mat2+=mat3).Print(7);
    std::cout<<std::endl;*/
    return 0;
}
