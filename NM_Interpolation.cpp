#include <iostream>
#include <math.h>
#include <cstdlib>
#include <vector>
using namespace std;

double* Nodes(double a, double b, int n);
double f(double x);
double* F_values(double* nodes, int n);
double Lagrange(double* nodes, double* f_values, double t, int n);
double* Lagranges_Array(double* nodes, double* f_values, int n,int nyzlov);


#define M_PI           double(3.14159265358979323846)
#define EPSILON  0.000001;

int main()
{
    cout << "6 uzlov:" << endl; 
    int n = 5;
    double* nodes = Nodes(0, 1.5, 5);
    for (int i = 0; i < n+1; i++) {
        cout << "x[" << i << "]" << " = " << nodes[i] << endl;
    }
    cout << "Znacheniya funktsii v 6 uzlah" << endl;
    double* f_values = F_values(nodes, 5);
    for (int i = 0; i < n + 1; i++) {
        cout << "f(x[" << i << "])" << " = " << f_values[i] << endl;
    }
    cout << "11 uzlov(tochek):" << endl;
    int n1 = 10;
    double* nodes11 = Nodes(0, 1.5, 10);
    for (int i = 0; i < n1 + 1; i++) {
        cout << "x[" << i << "]" << " = " << nodes11[i] << endl;
    }
    cout << "Znacheniya funktsii v 11 uzlah:" << endl;
    double* f_values11 = F_values(nodes11, 10);
    for (int i = 0; i < n1 + 1; i++) {
        cout << "f(x[" << i << "])" << " = " << f_values11[i] << endl;
    }
    int coltoch = 11;
    cout << "The Lagrange interpolation polynomial at 6 points:" << endl;
    double* lagranges_11 = Lagranges_Array(nodes, f_values, coltoch,n+1);
    for (int i = 0; i < coltoch; i++) {
        cout << "L_n(x[" << i << "])" << " = " <<lagranges_11[i] << endl;
    }
    cout << "Maximalnaya pogreshost'" << endl;
    for (int i = 0; i < coltoch; i++) {
        cout << "L_n(x[" << i << "])"<< "-" << "f(x[" << i << "])"<< "=" << abs(lagranges_11[i]-f_values11[i]) << endl;
    }

}

double* Nodes(double a, double b, int n) {
    double h = (b - a) / n;
    double* nodes = new double[n+1];
    for (int i = 0; i < n+1; i++) {
        double xi = a + i * h;
        nodes[i] = xi;
    }
    return nodes;
}

double f(double x) { //Я ничего не менял здесь, значения функции снова неверно считает
    double a = x;
    double F = a;
    double q;
    int n = 0;
    while (abs(a) >= 0.000001) {
        q = ((-1) * (M_PI / 2) * (M_PI / 2) * (4 * n + 1) * x * x * x * x) / ((2 * n + 2) * (2 * n + 1) * (4 * n + 5));
        a = a * q;
        F += a;
        n++;
    }
    return F;
}

double* F_values(double* nodes, int n) {
    double* f_values = new double[n + 1];
    for (int i = 0; i < n + 1; i++) {
        f_values[i] = f(nodes[i]);
    }
    return f_values;
}

double Lagrange(double* nodes, double* f_values, double t, int n) {
    double sum, prod;
    sum = 0;
    for (int j = 0; j < n; j++) {
        prod = 1;
        for (int i = 0; i < n; i++) {
            if (i != j) {
                prod = prod * (t - nodes[i]) / (nodes[j] - nodes[i]);
            }
        }
        sum = sum + f_values[j] * prod;
    }
    return sum;
}

double* Lagranges_Array(double* nodes, double* f_values, int coltoch,int colyzlov) {
    double h = (double)(1.5 - 0) / (coltoch - 1);
    double* langranges_11 = new double[coltoch];
    for (int i = 0; i < coltoch; i++) {
        double lnx = Lagrange(nodes, f_values, 0+i*h, colyzlov);
        langranges_11[i] = lnx;
    }
    return langranges_11;
}





