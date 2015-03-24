
#include <armadillo>
#include <iostream>
#include <cmath>
#include <assert.h>



using namespace arma;
using namespace std;

typedef vector<double> Values;

ostream& operator<<(ostream& os, const vector<double>& obj) {
    for(auto elem:obj) {
        os << elem << " "; 
    }
    os << endl;
    return os;
}

class QR {
    private:
        mat A;
        mat Ainit;
        mat reflection;
        mat Q;
        double eps;
        double m;
        vector<double> s;
        vector<double> b;

    public:
        static QR from_file(string filePath) {
            ifstream input(filePath.c_str());
            int m;
            vector<double> s;
            size_t size;

            input >> size;
            input >> m;

            mat A(size, size);
            for (size_t i = 0; i < size; i++) {
                for (size_t j = 0; j < size; j++) {
                    double temp;
                    input >> temp;
                    A(i,j) = temp;
                }
            }

            for (size_t i = 0; i < size; i++) {
                double temp;
                input >> temp;
                s.push_back(temp);
            }

            QR qr(A, m, s);

            return qr;
        }

        friend ostream& operator<<(ostream& os, const QR& obj);

        QR(mat A, int m, vector<double> s) {
            this->A = A;
            this->Ainit = mat(A);
            this->eps = pow(10.0, -m);
            this->s = s;
        }

        void findTheB() {
            for (int i = 0; i < s.size(); i++) {
                double temp = 0;
                for (int j = 0; j < s.size(); j++) {
                    temp += s[j] * A(i,j);
                }
                b.push_back(temp);
            }
        }

        void household() {
            findTheB();
            int n = A.n_rows;
            mat P;
            Q = mat(A).eye();

            Values u = Values(n, 0);
            double beta;
            double k;

            for (int r = 0; r < n - 1; r++) {
                //Constructia matrici P
                double omega = 0;
                for (int i = r; i < n; i++) {
                    double val = A(i,r);
                    omega += pow(val, 2);
                }
                /* cout << omega << " " << eps << endl; */
                /* assert(omega <= eps); */

                k = sqrt(omega);
                k = A(r, r) > 0 ? -k : k;

                beta = omega - k * A(r,r);

                u[r] = A(r,r) - k;

                for (int i = r+1; i < n; i++) {
                    u[i] = A(i,r);
                }
                
                //Transformarea colnelor
                for (int j = r+1; j < n; j++) {
                    double y = 0;
                    for (int i = r; i < n; i++) {
                        y += u[i] * A(i,j);
                    }
                    y /= beta;

                    for (int i = r; i < n; i++) {
                        A(i,j) = A(i,j) - y * u[i];
                    }
                }

                //Transformar coloanei r din A
                A(r,r) = k;
                for (int i = r + 1; i < n; i++) {
                    A(i,r) = 0;
                }

                double y = 0;
                for (int i = r; i < n; i++) {
                    y += u[i] * b[i];
                }
                y /= beta;

                for (int i = r; i < n; i++) {
                    b[i] = b[i] - y * u[i];
                }

                //Q = P[r]  * Q
                for (int j = 0; j < n; j++) {
                    y = 0;
                    for (int i = r; i < n; i++) {
                        y += u[i] * Q(i,j); 
                    }
                    y /= beta;

                    for (int i = r; i < n; i++) {
                        Q(i, j) = Q(i, j) - y * u[i]; 
                    }
                }

            }
            Row<double> row = Row<double>();
            mat bMat = mat(b);
            mat x = trans(Q) * bMat / A;
        }

};



ostream& operator<<(ostream& os, const QR& obj) {
    
    os << "R: "<< obj.A << endl;
    os << "Q: " << obj.Q << endl;

    mat newQ, newR;
    qr(newQ, newR, obj.Ainit);
    os << endl << endl << "Armadillo" << endl;
    os << "R" << newR << endl;
    os << "Q" << newQ;
    return os;
}


int main() {
    QR qr = QR::from_file("data.txt");
    qr.household();
    cout << qr;
}
