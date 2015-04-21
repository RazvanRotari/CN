#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <limits>

using namespace std;

template <class T>
struct Element {
    size_t pos;
    T value;
};

template <class T>
class Matrix {
    typedef vector<Element<T> > Row;
    typedef vector<Row> Mat;
    typedef Element<T> element;

    size_t numberOfRows;
    size_t numberOfColumns;
    public:
    bool simetric;
    Mat mat;    

public:
    static Matrix<T> from_file(ifstream& file, size_t size) {
        Matrix result;
        result.mat = Mat(size);
        result.numberOfRows = size;
        result.numberOfColumns = size;
        result.simetric = false;


        while (!file.eof()) {
            T val;
            size_t i,j;
            char temp;
            file >> val >> temp >> i >> temp >> j;    
            element pos = {j, val};
            result.mat[i].push_back(pos);
        }
        for (auto row:result.mat) {
            sort(row.begin(), row.end(), [](Element<T>& a, Element<T>& b){return a.pos > b.pos;});
        }
        return result;
    }

    static Matrix<T> symetricRand(size_t size) {
        std::default_random_engine generator;
        std::uniform_int_distribution<size_t> distribution(0, size);
        int dice_roll = distribution(generator);
        Matrix<T> result;
        result.mat = Mat(size);
        result.simetric = true;

        //for earch row
        for (size_t i=0; i < size;i++) {
            std::uniform_int_distribution<size_t> countDistribution(0, (size - i) % 10);
            size_t count = countDistribution(generator);
            for (size_t counter=0; counter < count; counter++) {
                std::uniform_int_distribution<size_t> countDistribution(i, size);
                size_t position = countDistribution(generator);

                std::uniform_real_distribution<T> valueDistribution(0, numeric_limits<T>::max());
                T value = valueDistribution(generator);

                Element<T> elem = {position, value};
                result.mat[i].push_back(elem);
            }
        }
        return result;
    }

    Matrix<T> operator+(const Matrix& mat) const {
        Matrix<T> result;
        const Matrix<T>& a = this->numberOfRows > mat.numberOfRows ? mat : *this;
        size_t n;
        if (this->numberOfRows > mat.numberOfRows) {
            result = *this;
            n = mat.numberOfRows;
        } else {
            result = mat;
            n = this->numberOfRows;
        }
        result.numberOfRows = n;
        for (size_t idx=0; idx < n; idx++) {
            auto& rowResult = result.mat[idx], rowA = a.mat[idx];
            for (auto elemA:rowA) {
                auto elemR = find_if(rowResult.begin(), rowResult.end(), [&elemA](Element<T> elemR){return elemR.pos == elemA.pos;});
                if (elemR == rowResult.end()) {
                    auto posR = find_if(rowResult.begin(), rowResult.end(), [&elemA](Element<T> elemR){return elemA.pos > elemR.pos;});
                    rowResult.insert(posR, elemA);
                    continue;
                };
                (*elemR).value += elemA.value;
            }
        }
       return result; 
    }

    vector<T> operator*(const vector<double> b) const {
        vector<T> result;
        size_t i = 0;
        for(auto row:mat) {
            T sum = 0;
            for(auto elem:row) {
                sum += elem.value * b[elem.pos];
            }
            result.push_back(sum);
            i++;
            

        }
        return result;
    }

    Matrix<T> operator*(const Matrix<T> B) const {
        Matrix<T> result;
        result.mat = Mat(this->numberOfRows);
        size_t columnIndex = 0;
        for (size_t columnIndex=0; columnIndex < B.numberOfRows; columnIndex++) {
            auto column = B.getcolumnForIndex(columnIndex);
            
            //Compute each column from the result
            size_t rowIndex = 0;
            for(auto row:this->mat) {
                /* cout << "Mul: " << "( " << rowIndex << ", "<< columnIndex << ")" << endl; */
                T value = multiply(row, column);
                if (value != 0) result.mat[rowIndex].push_back({columnIndex, value});
                rowIndex++;
                /* cout << endl; */
            }
        }
        return result;
    }

    T multiply(const vector<Element<T> >& row, const vector<Element<T> >& column) const {
        auto columnPrevPosition = column.begin();
        T sum = 0;
        /* cout << "Row:" << row; */
        /* cout << "column: "<< column; */
        for (auto elem:row) {
            size_t posX = elem.pos;
            auto elemB = find_if(columnPrevPosition, column.end(), [&posX](Element<T> elemB){return posX == elemB.pos;});
            if (elemB != column.end()) {
                //There is a value in both matrix
                sum += elem.value * (*elemB).value;
                columnPrevPosition = elemB;
            }

        }
        return sum;
    }

    vector<Element<T> > getcolumnForIndex(const size_t index) const {
        vector<Element<T>> column;
        //Get the column
        size_t rowIndex = 0;
        for (auto row:mat) { 
            auto elem = find_if(row.begin(), row.end(), [&index](Element<T> elem){return index == elem.pos;});
            if (elem != row.end()){
            Element<T> newElem = *elem;
            newElem.pos = rowIndex;

            column.push_back(newElem);

            };
            rowIndex++;
        }
        return column;
    }

    T getValueAtPosition(size_t i, size_t j)
    {
        auto row = mat[i]; 
        auto elem = find_if(row.begin(), row.end(), [&j](Element<T> elem){return j == elem.pos;});
        if(elem != row.end()){
            return elem->value;
        }
        return 0;
    }


    template <class T>
    vector<T> Gauss_Seidel(vector<T> b, double epsilon = pow(10, -5), size_t kmax = 1000) 
    {
        for(size_t i = 0; i  < numberOfRows; i++) 
        { 
            if(abs(getValueAtPosition(i, i)) <= epsilon)
            {
                cout << "Nu se poate rezolva sistemul folosind Gauss-Seidel" << endl;
                exit(0);   
            }
        }
        vector<T> xc(numberOfRows, 0), xp(numberOfRows, 0);
      //  vector<T> xc = {1, 2, 3, 4, 5};
        //vector<T> xp = {1, 2, 3, 4, 5};

        size_t k = 0;
        double deltax;
        do 
        {
            xp = xc;
            size_t index = 0;
            for(auto row:mat)
            {
                xc[index] = b[index];
                T sum = 0;
                for(auto elem:row)
                {
                    if (elem.pos >= index) break;
                    sum += elem.value * xc[elem.pos];
                }
                xc[index] -= sum;
                sum = 0;
                auto it = find_if(row.begin(), row.end(), [&index](Element<T> elem){return elem.pos > index;});
                for(auto val=it; val != row.end(); val++)
                {
                    sum += val->value * xp[val->pos];
                }
                xc[index] -= sum;
                xc[index] /= getValueAtPosition(index, index);
                index++;
            }
            T sum = 0; 
            for(size_t i = 0; i < xc.size(); i++)
            {
                sum += pow(abs(xc[i] - xp[i]), 1);
            }
            deltax = sqrt(sum);
            k++;
        } while(deltax >= epsilon && k <= kmax && deltax <= pow(10, 8));
        cout << deltax << " " << k << endl;
        if(deltax < epsilon)
        {
            cout << xc;
        }
        else
        {
            cout << "Divergenta.";
        }
        return xc;
    }
};

template <class T>
vector<T> readRow(ifstream& file, size_t size) {
    vector<T> result;
    for (size_t i=0; i<size; i++) {
        T val;
        file >> val;
        result.push_back(val);
    }
    return result;
}

template <class T>
ostream& operator<<(ostream& os, const vector<T>& obj) {
    for(auto elem:obj) {
        os << elem << " "; 
    }
    os << endl;
    return os;
}

template <class T>
ostream& operator<<(ostream& os, const vector<Element<T> >& obj) {
    for(auto elem:obj) {
        os << "{"<< elem.pos << ": " << elem.value << "}  "; 
    }
    os << endl;
    return os;
}

template <class T>
ostream& operator<<(ostream& os, const Matrix<T>& obj) {
    size_t index = 0;
    bool symetric = obj.simetric;
    for(auto row:obj.mat) {
        for (auto elem:row) {
            os << elem.value << " ," << index << ", "<< elem.pos << endl;
            if (symetric) os << elem.value << " ," << elem.pos << ", "<<  index << endl;
        }
        index++;
    }
    os << endl;
    return os;
}



int main() {
    ifstream file("m_rar_2015_1.txt");
    size_t size;
    file >> size;
    auto b = readRow<double>(file, size);
    auto A = Matrix<double>::from_file(file, size);
    cout << "Read" << endl;
    A.Gauss_Seidel(b);
   // cout << b << endl;
  //  cout << "A:" << A << endl;

}
