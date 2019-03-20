#ifndef USEFUL_H
#define USEFUL_H

#include "Dense"

typedef Eigen::MatrixXd mat;
typedef Eigen::VectorXd vec;

typedef double db;
typedef long double dbl;

// FUNCTIONS


void Zero_vec(vec &v,int rows);
void Zero_mat(mat &m,int rows, int cols);
void tridag(vec &a, vec &b, vec &c, vec &r, vec &sol); 
void fourier(vec &data, int isign);
void Rfourier(vec &data, int isign); // length of vector must be a power of 2
void modulate(vec &v, int rows);
void modulate1(vec &v,int rows);
db heaviside(db  R, db x);
db Htheta(db x);
db Power(db x,int n);


template <class T>
inline T Nan_(T x) {
    if (abs(x) < 1e-70 || x != x) {
        return 0;
    } else
        return x;
}

template <class T>
inline T sq(T x) {
    return (T) (x * x);
}

template <class T>
inline T cu(T x) {
    return (T) (x * x * x);
}

template <class T>
inline T recip(T x) {
    return (T) 1.0 / ((T) x);
}

template <class T>
inline T frac(T x, T y) {
    return (T) (x) / (T) y;
}

template <class T>
T dirac(T R, T x, T eps) {

    if (eq(R, x, eps)) {
        return 1.0;
    } else {
        return 0;
    }

}

template <class T>
inline void SWAP(T &a, T &b) {
    T temp = a;
    a = b;
    b = temp;
}


template<class T>
inline T cube_d(T x) {
    return x * x*x;
}

template <class T>
inline T sqr_d(T x) {

    return x*x;

}

template <class T>
T dirac(T R,T x);
template <class T>
inline bool eq(T a,T b,T eps){return abs(a-b)<eps;}
db pow12(db g);

// CONSTANTS
const db pi = 3.14159265358979323846; 
const db euler = 2.71828; 
//useful CLASSES


//{} ARRAY
class dbarray
{
public:
  dbarray();
  ~dbarray();
  void setelem(int,db);
  void printall();
  void resize(int);
  void setnumelems(int);
  void zero();
  void addelem(int);
  void divideelems(db);
  db getelem(int);
private:

  db *elems;
  size_t numelems;

};

#endif
