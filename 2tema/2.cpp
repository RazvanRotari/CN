#include <armadillo>
#include <iostream>
#include <cmath>
#include <assert.h>


using namespace arma;
using namespace std;

class LU {
    private:
        mat A;
        mat Ainit;
        double eps;
        double m;
        vector<double> b;
    public:
        static LU from_file(string filePath) {
            ifstream input(filePath.c_str());
            double m;
            vector<double> b;
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
                b.push_back(temp);
            }

            LU lu(A, m, b);

            return lu;
        }

        LU(mat A, double m, vector<double> b) {
            this->A = mat(A);
            this->Ainit = A;
            this->m = m;
            this->b = b;
            eps = pow(10, -m);
        }

        double div1(double a) {
            assert(abs(a) > eps);
            return 1/a;
        }
        friend ostream& operator<<(ostream& os, const LU& obj);
};


ostream& operator<<(ostream& os, const LU& obj) {
    os << obj.A << endl;
    os << obj.Ainit;

    return os;
}


int main() {

    LU lu = LU::from_file("data.txt");
    cout << lu;
    return 0;
}
