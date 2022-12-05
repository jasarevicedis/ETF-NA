#include <exception>
#include <iostream>
#include <stdexcept>
#include <sys/types.h>
#include <vector>
#include <cmath>
#include <limits>
#include <iomanip>



class Vector {
    std::vector<double> elements;
public:
    explicit Vector(int n);
    Vector(std::initializer_list<double> l);
    int NElems() const;
    double Norm() const;
    friend double VectorNorm(const Vector &v);
    double GetEpsilon() const;
    void Print(char separator , double eps) const;
    friend void PrintVector(const Vector &v, char separator,double eps);
    void Chop(double eps);
    bool EqualTo(const Vector &v, double eps) const;

    double &operator[](int i);
    double operator[](int i) const;
    double &operator()(int i);
    double operator()(int i) const;
    friend Vector operator +(const Vector &v1, const Vector &v2);
    Vector &operator +=(const Vector &v);
    friend Vector operator -(const Vector &v1, const Vector &v2);
    Vector &operator -=(const Vector &v);
    friend Vector operator *(double s, const Vector &v);
    friend Vector operator *(const Vector &v, double s);
    Vector &operator *=(double s);
    friend double operator *(const Vector &v1, const Vector &v2);
    friend Vector operator /(const Vector &v, double s);
    Vector &operator /=(double s);
};

/*---------------FUNKCIJE------------------*/
    Vector::Vector(int n){
        if(n<=0)
            throw std::range_error("Bad dimension");
        elements.resize(n);
    }
    Vector::Vector(std::initializer_list<double> l){
        if(l.size()==0)
            throw std::range_error("Bad dimension");
        elements.resize(l.size());
        auto it=l.begin();
        for(int i=0;i<l.size();i++){
            elements.at(i)=*it;
            it++;
        }
    }
    int Vector::NElems() const {
        return elements.size();
    }
    double Vector::Norm() const {
        double sum=0;
        for(int i=0;i<elements.size();i++)
            sum+=std::pow(elements[i],2);
        return std::sqrt(sum);
    }
    double VectorNorm(const Vector &v){
        double sum=0;
        for(int i=0;i<v.elements.size();i++)
            sum+=v.elements.at(i)*v.elements.at(i);
        return std::sqrt(sum);
    }
    double Vector::GetEpsilon() const {
        return 10 * this->Norm() * std::numeric_limits<double>::epsilon();
    }
    void Vector::Print(char separator = '\n', double eps = -1) const {
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
    void PrintVector(const Vector &v, char separator = '\n',double eps = -1){
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
    void Vector::Chop(double eps = -1){
        double cutParametar = eps;
        if(cutParametar < 0)
            cutParametar = GetEpsilon();
        for(auto el: elements){
            if(el<cutParametar)
                el = 0;
        }
    }
    bool Vector::EqualTo(const Vector &v, double eps = -1) const{
        if(elements.size() != v.elements.size())
            return false;
        for(int i=0;i<elements.size();i++){
            if(elements[i] != v.elements[i])
                return false;
        }    
        return true;
    }

    double &Vector::operator[](int i){
        return elements[i];
    }
    double Vector::operator[](int i) const {
        return elements[i];
    }
    double &Vector::operator()(int i) {
        if(i > elements.size())
            throw std::range_error("Invalid index");
        return elements[i-1];
    }
    double Vector::operator()(int i) const{
        if(i > elements.size())
            throw std::range_error("Invalid index");
        return elements[i-1];
    }
    Vector operator +(const Vector &v1, const Vector &v2){
        Vector v(v1);
        return v+=v2;
    }
    Vector &Vector::operator +=(const Vector &v){
        if(elements.size() != v.elements.size())
            throw std::domain_error("Incompatible formats");  
        for(int i=0;i<v.elements.size();i++)
            elements[i] += v.elements[i];
        return *this;
    }
    Vector operator -(const Vector &v1, const Vector &v2){
        Vector v(v1);
        return v-=v2;
    }
    Vector &Vector::operator -=(const Vector &v){
        if(elements.size() != v.elements.size())
            throw std::domain_error("Incompatible formats"); 
        for(int i=0;i<v.elements.size();i++)
            elements[i] -= v.elements[i];
        return *this;
    }
    Vector operator *(double s, const Vector &v){
        Vector tmp(v);
        return tmp*=s;
    }
    Vector operator *(const Vector &v, double s){
        Vector proizvod{};
        proizvod.elements.resize(v.elements.size());
        for(int i=0;i< v.elements.size();i++)
            proizvod.elements[i]= s * v.elements[i];
        return proizvod;
    }
    Vector &Vector::operator *=(double s){
        for(int i=0;i<elements.size();i++)
            elements[i]= elements[i]*s;
        return *this;
    }
    double operator *(const Vector &v1, const Vector &v2){
        if(v1.elements.size() != v2.elements.size())
            throw std::domain_error("Incompatible formats"); 
        double proizvod(0);
        for(int i=0;i< v1.elements.size();i++)
            proizvod+=v1.elements[i]*v2.elements[i];
        return proizvod;
    }
    Vector operator /(const Vector &v, double s){
       auto tmp(v);
       return tmp/=s;
    }
    Vector &Vector::operator /=(double s){
        if(s==0)
            throw std::domain_error("Division by zero");
        for(int i=0;i<elements.size();i++)
            elements[i]/=s;
        return *this;
    }


class Matrix {
    std::vector<std::vector<double>> elements;
    Matrix VectorToColumn(Vector v) {      
            Matrix mat(v.NElems(), 1);
            for(int i = 0; i < v.NElems(); i++)
                mat[i][0] = v[i];
            return mat;
        }


    public:
        Matrix(int m, int n);
        Matrix(const Vector &v);
        Matrix(std::initializer_list< std::vector<double> >l);
        int NRows() const;
        int NCols() const;
        double *operator[] (int i);
        const double *operator[] (int i) const;
        double &operator() (int i, int j);
        double operator ()(int i, int j) const;
        double Norm() const;
        friend double MatrixNorm(const Matrix &m);
        double GetEpsilon() const;
        void Print(int width = 2, double eps = -1) const;
        friend void PrintMatrix(const Matrix &m, int width, double eps);
        friend Matrix operator +(const Matrix &m1, const Matrix &m2);
        Matrix &operator +=(const Matrix &m);
        friend Matrix operator -(const Matrix &m1, const Matrix &m2);
        Matrix &operator -=(const Matrix &m);
        friend Matrix operator *(double s, const Matrix &m);
        friend Matrix operator *(const Matrix &m, double s);
        Matrix &operator *=(double s);
        friend Matrix operator *(const Matrix &m1, const Matrix &m2);
        Matrix &operator *=(const Matrix &m);
        friend Vector operator *(const Matrix &m, const Vector &v);
        friend Matrix Transpose(const Matrix &m);
        void Transpose();
        void Chop(double eps);
        bool EqualTo(const Matrix &m, double eps = -1) const;
        friend Matrix LeftDiv(Matrix m1, Matrix m2);

        friend Vector LeftDiv(Matrix m, Vector v);

        friend Matrix operator/(const Matrix& m, double s);
        Matrix &operator/=(double s);

        friend Matrix operator/(Matrix m1, Matrix m2);
        Matrix &operator/=(Matrix m);

        double Det() const;
        friend double Det(Matrix m);

        void Invert();
        friend Matrix Inverse(Matrix m);

        friend Matrix RREF(Matrix m);
        void ReduceToRREF();

        int Rank() const;
        friend int Rank(Matrix m);
};

/*METODE*/
Matrix::Matrix(int m, int n) {
    if( m <= 0 || n <= 0) throw std::range_error("Bad dimension");
        elements.resize(m);
    for(int i = 0; i < m; i++) {
        elements[i].resize(n);
        std::fill(elements[i].begin(), elements[i].end(), 0);
    }
}
Matrix::Matrix(const Vector &v) {
    elements.resize(1);
    elements[0].resize(v.NElems());
    for(int i = 0; i < v.NElems(); i++)
        elements[0][i] = v[i];
}
Matrix::Matrix(std::initializer_list< std::vector<double> >l) {
    if(l.size() == 0) throw std::range_error("Bad dimension");
    elements.resize(l.size());
    int i = 0;
    int velicina = 0;
    for(const auto &x : l)  {
        if(i == 0) velicina = x.size();
        if(velicina != x.size() || velicina == 0)
            throw std::logic_error("Bad matrix");
        elements[i++] = x;
    }
}
int Matrix::NRows() const {
    return elements.size();
}
int Matrix::NCols() const {
    return elements[0].size();
}
double *Matrix::operator[] (int i) {
    return &elements[i][0];
}
const double *Matrix::operator[] (int i) const {           
    return &elements[i][0];
}
double &Matrix::operator() (int i, int j) {
    if(i < 1 || j < 1 || i > elements.size() || j > elements[0].size())
        throw std::range_error("Invalid index");
    return elements[i - 1][j - 1];
}
double Matrix::operator ()(int i, int j) const {
    if(i < 1 || j < 1 || i > elements.size() || j > elements[0].size())
        throw std::range_error("Invalid index");
    return elements[i - 1][j - 1];
}
double Matrix::GetEpsilon() const {
    return 10 * Norm() * std::numeric_limits<double>::epsilon();
}
void PrintMatrix(const Matrix &m, int width=2, double eps=-1) {
    m.Print(width, eps);
}     
void Matrix::Chop(double eps = -1) {
    if(eps < 0) eps = GetEpsilon();
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
            if(std::fabs(elements[i][j]) < eps) elements[i][j] = 0;
}

/*metode*/

int Rank(Matrix m) {
    return m.Rank();
}

int Matrix::Rank() const {
    Matrix X(NRows(), NCols());
    X += *this;
    int k = -1;
    int l =  -1;
    int n = X.NCols();
    int m = X.NRows();
    std::vector<bool> w;
    w.resize(n);
    int p;

    for(int j = 0; j < n; j++) {
        w.at(j) = false;
    }

        while( k < m && l < n)
        {
            l++;
            k++;
            double v = 0;

            while( v < X.GetEpsilon() && l < n) {
                p = k;
                for(int i = k; i < m; i++) {
                    if(std::fabs(X.elements[i][l]) > v) {
                        v = std::fabs(X.elements[i][l]);
                        p = i;
                    }
                }
                if(v < X.GetEpsilon()) l++;
            }

            if(l < n) {
                w[l] = true;
                if(p != k) { 
                    X.elements[k].swap(X.elements[p]);
                }

                double mi = X.elements[k][l];
                for(int j = l; j < n; j++)
                    X.elements[k][j] /= mi;
                for(int i = 0; i < m; i++) {
                    if(i != k) {
                        mi = X.elements[i][l];
                        for(int j = l; j < n; j++)
                            X.elements[i][j] -= mi * X.elements[k][j];
                    }
                }
            }
        }
    return k;
}

Matrix RREF(Matrix m) {
    Matrix result(m.NRows(), m.NCols());
    result += m;
    result.ReduceToRREF();
    return result;
}

void Matrix::ReduceToRREF() {
    int k = -1;
    int l =  -1;
    int n = NCols();
    int m = NRows();
    std::vector<bool> w;
    w.resize(n);
    int p;

    for(int j = 0; j < n; j++) {
        w.at(j) = false;
    }

        while( k < m && l < n)
        {
            l++;
            k++;
            double v = 0;

            while( v < GetEpsilon() && l < n) {
                p = k;
                for(int i = k; i < m; i++) {
                    if(std::fabs(elements[i][l]) > v) {
                        v = std::fabs(elements[i][l]);
                        p = i;
                    }
                }
                if(v < GetEpsilon()) l++;
            }

            if(l < n) {
                w[l] = true;
                if(p != k) { 
                    elements[k].swap(elements[p]);
                }

                double mi = elements[k][l];
                for(int j = l; j < n; j++)
                    elements[k][j] /= mi;
                for(int i = 0; i < m; i++) {
                    if(i != k) {
                        mi = elements[i][l];
                        for(int j = l; j < n; j++)
                            elements[i][j] -= mi * elements[k][j];
                    }
                }
            }
        }
}

Matrix operator/(Matrix m1, Matrix m2) {
    if(m2.NCols() != m2.NRows()) throw std::domain_error("Divisor matrix is not square");
    else if(m1.NCols() != m2.NCols()) throw std::domain_error("Incompatible formats");
    else if(std::fabs(m2.Det()) < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Divisor matrix is singular");

    Matrix X (m1.NRows(), m1.NCols());

    int n = m2.NRows();
    int m = m1.NRows();

    for(int k = 0; k < n; k++) {

        //pivotizacija
        int p = k;
        for(int i = k + 1; i < n; i++)
            if(std::fabs(m2[k][i]) > std::fabs(m2[k][p])) p = i;

        if(p != k) {
            for (int i = 0; i < m2.NRows(); i++) {
                double t = m2[i][p];
                m2[i][p] = m2[i][k];
                m2[i][k] = t;
            }
        }

        // svodjenje na trougaonu matricu
        double mi;
        for(int i = k + 1; i < n; i++) {
            mi = m2[k][i] / m2[k][k];
            for(int j = k + 1; j < n; j++)
                m2[j][i] -= mi * m2[j][k];

            for(int j = 0; j < m; j++)
                m1[j][i] -= mi * m1[j][k];
        }
    }

     // supstitucija unazad
    double s;
    for(int k = 0; k < m; k++) {
        for(int i = n - 1; i >= 0; i--) {
            s = m1[k][i];
            for(int j = i + 1; j < n; j++)
                s -= m2[j][i] * X[k][j];
            X[k][i] = s / m2[i][i];
        }
    }
    return X;
}

Matrix & Matrix::operator/=(Matrix m) {
    *this = *this / m;
    return *this;
}

Matrix Inverse(Matrix m) {
   Matrix kopija(m.NRows(), m.NCols());
   kopija += m;
   kopija.Invert();
   return kopija;
}

void Matrix::Invert() {

    if(NCols() != NRows()) throw std::domain_error("Matrix is not square");
    if(std::fabs(Det()) < std::numeric_limits<double>::epsilon()) throw std::domain_error("Matrix is singular");

    int n = NRows();
    double mi;
    for(int k = 0; k < n; k++) {
        mi = elements[k][k];

        elements[k][k] = 1;
        for(int j = 0; j < n; j++)
            elements[k][j] /= mi;
        for(int i = 0; i < n; i++) {
            if(i != k) {
                mi = elements[i][k];
                elements[i][k] = 0;
                for(int j = 0; j < n; j++)
                    elements[i][j] -= mi * elements[k][j];
            }
        }
    }
}

double Det(Matrix m) {
    return m.Det();
}

double Matrix::Det() const{
    if(NCols() != NRows()) throw std::domain_error("Matrix is not square");
    double d =1;
    Matrix X(NRows(), NCols());
    X = *this;
    double p;
    for(int k = 0; k < NRows(); k++) {
        p = k;
        for(int i = k + 1; i < NRows(); i++)
            if(fabs(X.elements[i][k]) > fabs(X.elements[p][k])) p = i;

        if(fabs(X.elements[p][k]) < X.GetEpsilon()) return 0;
        if(std::fabs(p - k) > std::numeric_limits<double>::epsilon()) {
            d = (-1) * d;
            X.elements[p].swap(X.elements[k]);
        }
        d *= X.elements[k][k];


        for(int i = k + 1; i < NRows(); i++) {
            double mi = X.elements[i][k] / X.elements[k][k];
            for(int j = k + 1; j < NRows(); j++)
                X.elements[i][j] -= mi * X.elements[k][j];
        }
    }
    return d;
}

Matrix operator/(const Matrix& m, double s) {
    Matrix result(m.NRows(), m.NCols());
    if(s < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Division by zero");
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        result[i][j] = m[i][j] / s;
    return result;
}

Matrix & Matrix::operator/=(double s) {
    if(s < std::numeric_limits<double>::epsilon() ) throw std::domain_error("Division by zero");
     for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
        elements[i][j] /= s;
    return *this;
}

Vector LeftDiv(Matrix m, Vector v) {
    Matrix m2 (v.NElems(),1);
    m2 = m2.VectorToColumn(v);
    m2 = LeftDiv(m,m2);
    Vector result(v.NElems());
    for(int i = 0; i < v.NElems(); i++)
        result[i] = m2[i][0];
    return result;
}


Matrix LeftDiv(Matrix m1, Matrix m2) {
    if(m1.NCols() != m1.NRows()) 
        throw std::domain_error("Divisor matrix is not square");
    else if(m1.NCols() != m2.NRows()) 
        throw std::domain_error("Incompatible formats");

    Matrix X (m2.NRows(), m2.NCols());
    int n = m1.NRows();
    int m = m2.NCols();
    int p; 
    for(int k = 0; k < n; k++) {

        //pivotizacija
        p = k;
        for(int i = k + 1; i < n; i++)
            if(std::fabs(m1[i][k]) > std::fabs(m1[p][k])) p = i;
        if(std::fabs(m1[p][k]) < m1.GetEpsilon()) throw std::domain_error("Divisor matrix is singular");

        if(p != k) {
            m1.elements[p].swap(m1.elements[k]);
            m2.elements[p].swap(m2.elements[k]);
        }

        // svodjenje na trougaonu matrica
        double mi;
        for(int i = k + 1; i < n; i++) {
            mi = m1[i][k] / m1[k][k];
            for(int j = k + 1; j < n; j++)
                m1[i][j] -= mi * m1[k][j];

            for(int j = 0; j < m; j++)
                m2[i][j] -= mi * m2[k][j];
        }

    }

     // supstitucija unazad
    double s;
    for(int k = 0; k < m; k++) {
        for(int i = n - 1; i >= 0; i--) {
            s = m2[i][k];
            for(int j = i + 1; j < n; j++)
                s -= m1[i][j] * X[j][k];
            X[i][k] = s / m1[i][i];
        }
    }
    return X;
}


bool Matrix:: EqualTo(const Matrix &m, double eps) const {
    if(NCols() != m.NCols() || NRows() != m.NRows()) return false;
    if(eps < 0) eps = GetEpsilon();
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
            if(elements[i][j] - m[i][j] > eps) return false;
    return true;
}

void Matrix::Transpose() {

   if(NCols() == NRows()) {
        for(int i = 0; i < NRows(); i++)
            for(int j = i + 1; j < NCols(); j++)
               std::swap(elements[i][j], elements[j][i]);
   }
    else {
        Matrix pomocna(NCols(), NRows());
        for(int i = 0; i < NRows(); i++)
            for(int j = 0; j < NCols(); j++)
                pomocna[j][i] = elements[i][j];
        *this = pomocna;
    }
}

Matrix Transpose(const Matrix &m) {
    Matrix trans_matrix(m.NCols(), m.NRows());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
            trans_matrix[j][i] = m[i][j];
    return trans_matrix;
}

Vector operator *(const Matrix &m, const Vector &v) {
    if(m.NCols() != v.NElems()) throw std::domain_error("Incompatible formats");
    Vector vec (m.NRows());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < v.NElems(); j++)
            vec[i] += m[i][j] * v[j];
    return vec;
}

inline Matrix& Matrix::operator *=(const Matrix &m) {
    *this = *this * m;
    return *this;
}

Matrix operator *(const Matrix &m1, const Matrix &m2) {
    if(m1.NCols() != m2.NRows()) throw std::domain_error("Incompatible formats");
    Matrix proizvod_matrica (m1.NRows(), m2.NCols());
    for(int i = 0; i < m1.NRows(); i++)
        for(int j = 0; j < m2.NCols(); j++)
            for(int k = 0; k < m1.NCols(); k++)
                proizvod_matrica[i][j] += m1[i][k] * m2[k][j];
    return proizvod_matrica;
}

void Matrix::Print(int width, double eps) const {
    double prag = eps;
    if(eps < 0) prag = GetEpsilon();
    for(int i = 0; i < NRows(); i++) {
        for(int j = 0; j < NCols(); j++) {
            if(fabs(elements[i][j]) < prag) std::cout << std::setw(width) << '0';
            else if(elements[i][j] < 0) std::cout << std::setw(width + 1) << elements[i][j];
            else std::cout << std::setw(width) << elements[i][j];
        }
    std::cout << std::endl;
    }
}

double Matrix::Norm() const {
    double sumOfSquares(0);
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
        sumOfSquares += elements[i][j] * elements[i][j];
    sumOfSquares = sqrtl(sumOfSquares);
    return sumOfSquares;
}

double MatrixNorm(const Matrix &m) {
    return m.Norm();
}

Matrix& Matrix::operator *=(double s) {
    for(int i = 0; i < NRows(); i++)
        for(int j = 0; j < NCols(); j++)
        elements[i][j] *= s;
    return *this;
}

Matrix& Matrix::operator +=(const Matrix &m) {
    if((NRows() != m.NRows() ) || (NCols() != m.NCols())) throw std::domain_error("Incompatible formats");
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        elements[i][j] += m[i][j];
    return *this;
}

Matrix& Matrix::operator -=(const Matrix &m) {
    if((NRows() != m.NRows() ) || (NCols() != m.NCols())) throw std::domain_error("Incompatible formats");
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        elements[i][j] -= m[i][j];
    return *this;
}

Matrix operator *(double s, const Matrix &m) {
    Matrix m3(m.NRows(), m.NCols());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        m3[i][j] = m[i][j] * s;
    return m3;
}

Matrix operator *(const Matrix &m, double s) {
    Matrix m3(m.NRows(), m.NCols());
    for(int i = 0; i < m.NRows(); i++)
        for(int j = 0; j < m.NCols(); j++)
        m3[i][j] = m[i][j] * s;
    return m3;
}

Matrix operator +(const Matrix &m1, const Matrix &m2) {
    if((m1.NRows() != m2.NRows() ) || (m1.NCols() != m2.NCols())) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(), m1.NCols());
    for(int i = 0; i < m1.NRows(); i++)
        for(int j = 0; j < m1.NCols(); j++)
        m3[i][j] = m1[i][j] + m2[i][j];
    return m3;
}

Matrix operator -(const Matrix &m1, const Matrix &m2) {
    if((m1.NRows() != m2.NRows() ) || (m1.NCols() != m2.NCols())) throw std::domain_error("Incompatible formats");
    Matrix m3(m1.NRows(), m1.NCols());
    for(int i = 0; i < m1.NRows(); i++)
        for(int j = 0; j < m1.NCols(); j++)
        m3[i][j] = m1[i][j] - m2[i][j];
    return m3;
}

/*LUDecomposer*/
class LUDecomposer {
    Matrix mat;
    Vector w;// displays exhanges between rows
    void TestSquare(){
        if(mat.NCols() != mat.NRows()){
            throw std::domain_error("Matrix is not square");
        }
    }
    void TestFormat(Matrix mat1, Matrix  mat2) const {
        if((mat2.NRows() != mat1.NRows()) || ( mat2.NCols() != mat1.NCols()))
            throw std::domain_error("Incompatible formats");
    }
public:
    LUDecomposer(Matrix m); 
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Matrix GetCompactLU() const;
    Matrix GetL() const;
    Matrix GetU() const;
    Vector GetPermuation() const;
};

/*CONSTRUCTOR THAT DOES LU FACTORISATION*/
    LUDecomposer::LUDecomposer(Matrix m) : mat(m), w(m.NCols()) {
        TestSquare();
        int rows = mat.NRows();

        //Crout algorithm with pivotisation
        for(int j=0;j<rows;j++){
            for(int i=0;i<=j;i++){
                double s = mat[i][j];
                for(int k=0;k<i;k++)
                    s-= mat[i][k] * mat[k][j];
                mat[i][j] = s;
            }
            int p = j;
            for(int i = j + 1; i < rows; i++) {
                double s = mat[i][j];
                for(int k = 0; k <=  j - 1; k++)
                    s -= mat[i][k] * mat[k][j];
                mat[i][j] = s;
                if(std::fabs(s) > std::fabs(mat[p][j]))
                    p = i;
            }
            if(std::fabs(mat[p][j]) < mat.GetEpsilon())
                throw std::domain_error("Matrix is singular");
            if(p != j) {
                for(int r = 0; r < rows; r++) {
                    double temp = mat[p - 1][r];
                    mat[p - 1][r] = mat[j][r];
                    mat[j][r] = temp;
                }
            }
            w[j] = p;
            double mi = mat[j][j];
            for(int i = j + 1; i < rows; i++) {
                mat[i][j] /= mi;
            }
        }
    }
/*4 SOLVING METHODS FOR SYSTEMS WITH DIFFERENT RIGHT SIDE*/

    void LUDecomposer::Solve(const Vector &b, Vector &x) const{
        if(b.NElems() != mat.NRows()) 
            throw std::domain_error("Incompatible formats");
        else if(b.NElems() != x.NElems()) 
            throw std::domain_error("Incompatible formats");
        x = std::move(Solve(b));
    }
    Vector LUDecomposer::Solve(Vector b) const{
        
        int elems = b.NElems();
        if(elems != mat.NRows()) 
            throw std::domain_error("Incompatible formats");

        for(int i = 0; i < elems; i++) {
            int p = w[i];
            long double s = b[p];
            b[p] = b[i];    
            for(int j = 0; j < i; j++)
                s -= mat[i][j] * b[j];
            b[i] = s;
        }
        //supstitucija unazad
        for(int i = elems - 1; i >= 0; i--) {
            double s = b[i];
            for(int j = i + 1; j < elems; j++)
                s -= mat[i][j] * b[j];
            b[i] = s / mat[i][i];
        }      
        return b;
    }
    inline void LUDecomposer::Solve(Matrix &b, Matrix &x) const{
        TestFormat(b, x);
        int rows = b.NRows();
        int cols = b.NCols();

        for(int m = 0; m < cols; m++) {
            for(int i = 0; i < rows; i++) {
                int p = w[i];
                long double s = b[p][m];
                b[p][m] = b[i][m];
                for(int j = 0; j < i; j++)
                    s -= mat[i][j] * b[j][m];
                b[i][m] = s;
            }
            //supstitucija unazad
            for(int i = rows - 1; i >= 0; i--) {
                double s = b[i][m];
                for(int j = i + 1; j < rows; j++)
                    s -= mat[i][j] * x[j][m];
                x[i][m] = s /  mat[i][i];
            }
        }
    }
    Matrix LUDecomposer::Solve(Matrix b) const{
        TestFormat(b,mat);
        Matrix X(b.NRows(), b.NCols());
        Solve(b,X);
        return X;
    }
/*GETERS*/
    Matrix LUDecomposer::GetCompactLU() const{
        return mat;
    }
    Matrix LUDecomposer::GetL() const{
        Matrix L(mat.NRows(), mat.NCols());
        int n = L.NRows();
        for(int i = 0; i < n; i++)
            for(int j = 0; j <= i; j++) {
                if(i == j) 
                    L[i][j] = 1;
                else
                    L[i][j] = mat[i][j];
            }
        return L;
    }
    Matrix LUDecomposer::GetU() const{
        Matrix U(mat.NRows(), mat.NCols());
        int n =U.NRows();
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                if(j >= i)
                    U[i][j] = mat[i][j];
        return U;
    }
    //needs to return double vector 
    Vector LUDecomposer::GetPermuation() const{        
        int elems = w.NElems();
        Vector exchanges(elems);
        for(int i=0; i < elems; i++)
            exchanges[i] = w[i];
        return exchanges;
    }






/*QRDecomposer*/
class QRDecomposer {
    Matrix mat;
    Matrix vectors;
    Vector R_diagonal_elements;
public:
    void TestSquare() const {
        if(mat.NCols() != mat.NRows())
            throw std::domain_error("Matrix is not square");
    }
    QRDecomposer(Matrix m);
    void Solve(const Vector &b, Vector &x) const;
    Vector Solve(Vector b) const;
    void Solve(Matrix &b, Matrix &x) const;
    Matrix Solve(Matrix b) const;
    Vector MulQWith(Vector v) const;
    Matrix MulQWith(Matrix m) const;
    Vector MulQTWith(Vector v) const;
    Matrix MulQTWith(Matrix m) const;
    Matrix GetQ() const;
    Matrix GetR() const;
};


QRDecomposer::QRDecomposer(Matrix m) : mat(m), vectors(m.NRows(),m.NCols()), R_diagonal_elements(mat.NRows()) {
    int rows = m.NRows();
    int cols = m.NCols();

    if(rows < cols) 
        throw std::domain_error("Invalid matrix format");

    /*HOUSEHOLDER ALGORITHM*/
    for(int k = 0; k < cols; k++) {
        long double s = 0;
        for(int i = k; i < rows; i++)
            s += mat[i][k] * mat[i][k];
        s = std::sqrt(s);
        long double mi = std::sqrt(s * (s + std::fabs(mat[k][k])));
        if(std::fabs(mi) < mat.GetEpsilon())
            throw std::domain_error("Matrix is singular");
        if(mat[k][k] < 0) s = -s;
            mat[k][k] = (mat[k][k] + s) / mi;
        for(int i = k + 1; i < cols; i++)
            mat[i][k] = mat[i][k] / mi;
        R_diagonal_elements[k] = -s;
        for(int j = k + 1; j < cols; j++) {
            s = 0;
            for(int i = k; i < rows; i++)
                s += mat[i][k] * mat[i][j];
            for(int i = k; i < rows; i++)
                mat[i][j] -= s * mat[i][k];
        }
    }
}

    void QRDecomposer::Solve(const Vector &b, Vector &x) const{
        if(b.NElems() != mat.NRows()) 
            throw std::domain_error("Incompatible formats");
        else if(b.NElems() != x.NElems()) 
            throw std::domain_error("Incompatible formats");
        x = std::move(Solve(b));
    }
    Vector QRDecomposer::Solve(Vector b) const{
        TestSquare();
        if(b.NElems() != mat.NRows())
            throw std::domain_error("Incompatible formats");

        int rows = mat.NRows();
        for(int k = 0; k < b.NElems(); k++) {
            double s = 0;
            for(int i = k; i < rows; i++)
                s += mat[i][k] * b[i];
            for(int i = k; i < rows; i++)
                b[i] -= s * mat[i][k];
        }
        int elems = b.NElems();
        //suspstitucija unazad
        for(int i = elems - 1; i >= 0; i--) {
            double s = b[i];
            for(int j = i + 1; j < elems; j++)
                s -= mat[i][j] * b[j];
            b[i] = s / R_diagonal_elements[i];
        }
        return b;
    }
    void QRDecomposer::Solve(Matrix &b, Matrix &x) const{
        TestSquare();
        if(b.NRows() != mat.NRows() || b.NCols() != mat.NCols())
            throw std::domain_error("Incompatible formats");
        if(b.NRows() != x.NRows() || b.NCols() != x.NCols())
            throw std::domain_error("Incompatible formats");
        x = Solve(b);
    }
    Matrix QRDecomposer::Solve(Matrix b) const{
        TestSquare();
        if(b.NRows() != mat.NRows() || b.NCols() != mat.NCols())
            throw std::domain_error("Incompatible formats");
        int rows = b.NRows();
        int cols = b.NCols();
        Matrix x(rows, cols);

        for(int i = 0; i < cols; i++) {
            Vector v(rows);
            for(int j = 0; j < rows; j++) {
                v[j] = b[j][i];
            }
            v = Solve(v);
            for(int j = 0; j < rows; j++) {
                x[j][i] = v[j]; 
            }
        }
        return x;
    }
    Vector QRDecomposer::MulQWith(Vector v) const{
        TestSquare();
        if(v.NElems() != mat.NRows())
            throw std::domain_error("Incompatible formats");

        int rows = mat.NRows();
        for(int k = v.NElems() - 1; k >= 0; k--) {
            double s = 0;
            for(int i = k; i < rows; i++)
                s += mat[i][k] * v[i];
            for(int i = k; i < rows; i++)
                v[i] -= s * mat[i][k];
        }
        return v;
    }
    Matrix QRDecomposer::MulQWith(Matrix m) const{
        if(mat.NCols() != m.NRows())
            throw std::domain_error("Incompatible formats");

        int rows = mat.NRows();
        int cols = mat.NCols();
        int m_rows = m.NRows();
        for(int j = 0; j < m_rows; j++) {
            for(int k = cols - 1; k >= 0; k--) {
                double s = 0;
                for(int i = k; i < rows; i++)
                    s += mat[i][k] * m[i][j];

                for(int i = k; i < rows; i++)
                    m[i][j] -= s * mat[i][k];
            }
        }
        return m;
    }
    Vector QRDecomposer::MulQTWith(Vector v) const{
        TestSquare();
        if(v.NElems() != mat.NRows())
            throw std::domain_error("Incompatible formats");

        int rows = mat.NRows();
        for(int k = 0; k < v.NElems(); k++) {
            double s = 0;
            for(int i = k; i < rows; i++)
                s += mat[i][k] * v[i];
            for(int i = k; i < rows; i++)
                v[i] -= s * mat[i][k];
        }
        return v;
    }
    Matrix QRDecomposer::MulQTWith(Matrix m) const{
        if(mat.NCols() != m.NRows())
        throw std::domain_error("Incompatible formats");
        
        int rows = mat.NRows();
        int m_rows = m.NRows();//faster

        for(int j = 0; j < m_rows; j++) {
            for(int k = 0; k < m_rows; k++) {
                double s = 0;
                for(int i = k; i < rows; i++)
                    s += mat[i][k] * m[i][j];

                for(int i = k; i < rows; i++)
                     m[i][j] -= s * mat[i][k];
            }
        }
        return m;
    }
    Matrix QRDecomposer::GetQ() const{
        int rows =mat.NRows();//to prevent calling function in loops
        int cols = mat.NCols();//to prevent calling function in loops
        Matrix Q(rows,rows);

        for(int j = 0; j < rows; j++) {
            for(int i = 0; i < rows; i++)
                Q[i][j] = 0;
            Q[j][j] = 1;
            for(int k = cols - 1; k >= 0; k--) {
                double s = 0;
                for(int i = k; i < rows; i++)
                    s += mat[i][k] * Q[i][j];

                for(int i = k; i < rows; i++)
                    Q[i][j] -= s * mat[i][k];
            }
        }
        return Q;
    }
    Matrix QRDecomposer::GetR() const{
        int rows =mat.NRows();//to prevent calling function in loops
        int cols = mat.NCols();//to prevent calling function in loops
        Matrix R(rows,cols);

        for(int i=0;i < rows;i++){
            for(int j=0; j < cols; j++){
                if(i == j)
                    R[i][j] = R_diagonal_elements[i];
                else if(i > j)
                    R[i][j] = 0;
                else 
                    R[i][j] = mat[i][j];
            }
        }
        return R;
    }

/*UNIT TESTS*/
/*
void QRDecomposerConstructorTest(){
    try{
        QRDecomposer mat1({{1,1,1}, {2,2,2}, {3,3,3}});//singular
    }catch(exception e){
        std::cout << "matrica je singularna";
    }
    QRDecomposer mat2({{}, {}});//not square matrix
    QRDecomposer mat3({{}, {}});// square, not singular
    
if()
}*/
    
int main(){
    /*MATRIX ADDITION METHODS*/  
    /*
Matrix mat({{2,4,6}, {2 ,4,6}});
mat /= 2;
mat.Print();

try {
    mat/= 0;
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what();
}

Matrix mat1({{1,1,1}, {2,3,4}});
Matrix mat2({{1,2,3}, {3,5,7}, {5, 9, 7}});
try {
    mat1/=mat2;
}catch(std::domain_error &e ){
    std::cout << "Exception: " << e.what();
}

Matrix result = LeftDiv(mat1,mat2);
result.Print();

try{
    mat1.Det();
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what();
}
Matrix matt{{1,2,3},{4,5,6},{7,8,9}};
matt.ReduceToRREF();
Matrix m{{1,0,-1},{0,1,2},{0,0,0}};
if(matt.EqualTo(m))
    std::cout << "RREF WORKING" << std::endl;

double e =std::numeric_limits<double>::epsilon();
Matrix mat3({{1,e,4}, {4,e,3}, {1,3,4}} );
Matrix choped({{1,0,4}, {4,0,3}, {1,0,4}});
if(mat3.EqualTo(choped))
    std::cout << "chop working!" << std::endl;

Matrix mat5({{1,2,3}, {7,9,1}, {5,1,8}});
std::cout << mat5.Rank() << std::endl;

try {
    mat1.Invert();
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what() << std::endl;
}

Matrix mat6 {{3,0,2},{0,0,1},{2,-2,1}};
Matrix inverse{{0.2,-0.2, 0.2},{0.2,0.3,-0.3},{0,1,0}};
if(inverse.EqualTo(Inverse(mat6))) 
std::cout << "inverse works";*/
    /*LU testing*/
/*
try {
    LUDecomposer({{1,1,1}});
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what() << std::endl;
}

try{
    LUDecomposer({{1,1,1}, {2,2,2}, {3,3,3}});
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what() << std::endl;
}
try {
    LUDecomposer({{1,2,3}, {4,5,6}, {8,10,9}}).Solve({1,1,1,1,1});
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what();
}

*/
Matrix A{{0,7,2},{4,6,5},{2,1,7}};
Matrix x{{4,2,3},{2,3,1},{8,11,9}};
Matrix rez=A*x;
LUDecomposer lu(A);
lu.Solve(rez,rez);
if(rez.EqualTo(lu.GetL()*lu.GetU())) 
    std::cout<<"LU works!\n";
/*QR testing*/
/*
try {
    QRDecomposer qr({{1,1,1,1,1}, {1,1,1,1,1}});
}catch(std::domain_error &e){
    std::cout << "Exception: " << e.what();
}
try {
    QRDecomposer qr({{1,1,1}, {2,2,2}, {3,3,3}});
}catch(std::domain_error &e){
    std::cout << "Exeption: " << e.what();
}*/
Matrix B{{2.34,3.45,4.56},{5.67,6.78,7.89},{8.91,9.12,1.23}};
QRDecomposer qr(A);
if ((qr.GetQ()*qr.GetR()).EqualTo(B)) 
    std::cout<<"qr works"<<std::endl;

Matrix A1{{2,2,3}, {4,5,5}, {7,7,8}};
QRDecomposer qr1(A1);
Vector v{1,5,1};
Vector x2=qr1.Solve(v);
    if (x2.EqualTo({-5,4,1})) std::cout<<"solve works"<<std::endl;
return 0;
}
