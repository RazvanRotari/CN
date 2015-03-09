#include <armadillo>
#include <iostream>
#include <cmath>
#include <assert.h>


using namespace arma;
using namespace std;

ostream& operator<<(ostream& os, const vector<double>& obj) {
    for(auto elem:obj) {
        os << elem << " "; 
    }
    os << endl;
    return os;
}

class LU {
    private:
        mat A;
        mat Ainit;
        double eps;
        double m;
        vector<double> b;
        vector<double> solution;
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
            assert(check());
        }

        double div(double a, double b) {
            assert(abs(b) > eps);
            return a/b;
        }

        bool check() {
            for (size_t i = 0; i < A.n_rows; i++) {
                mat submat = A.submat(0, 0, i, i);
                double d = det(A);
                if (d == 0) return false;
            }
            return true;
        }

        void LUdecomposition() {
            for (int p = 0; p < A.n_rows; p++) {
                calculateU(p);
                calculateL(p);
            }
        
        }

        void calculateU(int p) {
            for (int i = p; i < A.n_rows; i++) {
                double u = Ainit(p, i);
                double sum = 0;
                for (int k = 0; k < p; k++) {
                    sum += getL(p, k, p) * getU(k, i, p);
                }
                u -= sum;
                A(p, i) = u;
            }
        }

        void calculateL(int p) {
            for (int i = p; i < A.n_rows; i++) {
                double l = Ainit(i, p);
                double sum = 0;
                for (int k = 0; k < p; k++) {
                    sum += getL(i, k, p) * getU(k, p, p);
                }
                l -= sum;
                double temp  =  div(l, A(p, p));
                if (i == p) continue;
                A(i, p) = temp;
            }
        }

        double getL(int i, int j, int p = -1) {
            if (i == j) return 1;
            if (j == p) return 0;
            return A(i, j);
        }

        double getU(int i, int j, int p = -1) {
            if (i == p) return 0;
            return A(i, j);

        }

        double LUdet() {
            double prod = 1;
            for (int i = 0; i < A.n_rows; i++) {
                prod *= A(i, i);
            }
            return prod * 1; // U * L;
        }

        vector<double> findTheX() {
            double size = A.n_rows;
            vector<double> ys = vector<double>(size, 0);
            ys[0] = div(b[0], getL(0, 0));

            for (int i = 1; i < size; i++) {
                double temp = 0;
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += getL(i, j) * ys[j];
                }
                temp = b[i] - sum;
                temp = div(temp, getL(i, i));
                ys[i] = temp;
            }

            vector<double> xs = vector<double>(size, 0);
            double n = size - 1;
            xs[n] = div(ys[n], getU(n ,n));

            for (int i = n - 1; i >= 0; i--) {
                double temp = 0;
                double sum = 0;
                for (int j = i + 1; j < size; j++) {
                    sum += getU(i, j) * xs[j];
                }
                temp = ys[i] - sum;
                temp = div(temp, getU(i, i));
                xs[i] = temp;
            }
            solution = xs;
            return xs;
        }

        double findTheEuclidean() {
            double errorSum = 0;
            for (int i = 0; i < A.n_rows; i ++) {
                double sum = 0;
                for (int j = 0; j < A.n_rows; j++) {
                    sum += Ainit(i, j) * solution[j];
                }
                double temp = sum - b[i];
                errorSum += pow(abs(temp), 2);
            }
            return sqrt(errorSum);


        }

        void solveSystem() {
            LUdecomposition();
            findTheX();
            cout << "Error is: " << findTheEuclidean() << endl;
        }

    public:
        friend ostream& operator<<(ostream& os, const LU& obj);
};


ostream& operator<<(ostream& os, const LU& obj) {
    os << obj.A << endl;

    return os;
}

int main() {

    LU lu = LU::from_file("data.txt");
    lu.solveSystem();
    return 0;
}
