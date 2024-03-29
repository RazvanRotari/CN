#include <vector>
#include <iostream>
#include <cmath>
#include <random>
#include <utility>
#include <algorithm>

#include <armadillo>

using namespace std;
using namespace arma;

using Point = pair<double, double>;
using Values = vector<double>;
using DataPoints = vector<Point>;


ostream& operator<<(ostream& os, const Point& elem) {
    os << "(" << elem.first << ", " << elem.second << ") "; 
    return os;
}

ostream& operator<<(ostream& os, const DataPoints& obj) {
    os << "[";
    for(auto elem:obj) {
        os << "(" << elem.first << ", " << elem.second << ") "; 
    }
    os << "]" << endl;
    return os;
}

double function(double input) {
    /* return pow(input, 3) + 2 * pow(input, 2) + input - 1; */
    return 2 * pow(input, 2) + input - 1;
}

double derivedFunction(double input) {
    /* return 3 * pow(input, 2) + 4 * input + 1; */
    return 4 * input + 1;
}

DataPoints generateDataPoints(double start, double end, size_t count) {
    random_device rd;
    default_random_engine random_generator(rd());
    uniform_real_distribution<double> distribution(start, end);
    DataPoints dataPoints;
    dataPoints.push_back(make_pair(start, function(start)));

    double h = (end-start) / (count + 1);
    for (size_t index=0; index < count; index++) {
        /* double input = distribution(random_generator); */
        double input = start + (index + 1) * h;

        double output = function(input);
        dataPoints.push_back(make_pair(input, output));
    }
    dataPoints.push_back(make_pair(end, function(end)));

    sort(dataPoints.begin(),dataPoints.end(), [](const Point& a, const Point& b){return (a.first < b.first);});
    return dataPoints;
 }

double lagrange(DataPoints dataPoints, double input) {
    Values values;
    values.reserve(dataPoints.size());
    for (auto elem:dataPoints) {
        values.push_back(elem.second);
    }

    for(size_t step = 1; step < dataPoints.size(); step++) {
        for(size_t index = 0; index < values.size() - 1; index++) {
            double first = dataPoints[index].first;
            double last = dataPoints[index+step].first;
            values[index] = ((last - input) * values[index] + (input - first) * values[index+1]) / (last - first);
        }
        values.pop_back();
    }
    return values[0];
}

double spline(DataPoints dataPoints, double input, double derivedFirst, double derivedLast) {
    size_t size = dataPoints.size();

    auto x = [&dataPoints](size_t index){return dataPoints.at(index).first;};
    auto y = [&dataPoints](size_t index){return dataPoints.at(index).second;};

    auto h = [&x](size_t index){return x(index + 1) - x(index);};

    sp_mat H(size, size);
    H(0,0) = 2 * h(0);
    H(0, 1) = h(0);

    vec f(size);

    double val = (y(1) - y(0)) / h(0);
    val -= derivedFirst;
    val *= 6;
    f[0] = val;

    val = derivedLast;
    val -= (y(size - 1) - y(size - 2)) / h(size - 2);
    val *= 6;
    f[size - 1] = val;

    for (size_t index = 1; index < size - 1; index++) {
        double hindex = h(index);
        double hprev = h(index-1);
        H(index, index - 1) = hprev;
        H(index, index) = 2 * (hprev + hindex);
        H(index, index + 1) = hindex;

        double val = 0;
        val = (y(index + 1) - y(index)) / hindex;
        val -= (y(index) - y(index - 1)) / hprev;
        val *= 6;
        f[index] = val;
    }
    H(size - 1, size - 2) = h(size - 2);
    H(size - 1, size - 1) = 2 * h(size - 2);

    auto A = spsolve(H, f);
    cout << A;

    size_t i0 = 0; 
    for(size_t index = 0; index < size - 1; index++) {
        if (dataPoints[index].first <= input && input <= dataPoints[index + 1].first) {
            i0 = index;
            break;
        }
    }
    double alpha = y(i0+1) - y(i0);
    alpha /= h(i0);
    alpha -= h(i0) * (A(i0 + 1) - A(i0)) / 6;

    double beta = (x(i0+1) * y(i0) - x(i0) * y(i0 + 1)) / h(i0);
    beta -= h(i0) * (x(i0+1) * A(i0) - x(i0) * A(i0+1)) / 6;

    double ret = pow(input - x(i0), 3) * A(i0 + 1);
    ret /= 6 * h(i0);
    ret += (pow(x(i0 + 1) - input, 3) * A(i0)) / (6 * h(i0));


    ret += alpha * input;
    ret += beta;

    return ret;
}

int main() {
    double start = 0, end = 8;
    double count = 10;
    double input = 2;

    /* cout << "start: "; */
    /* cin >> start; */
    /* cout << "end: "; */
    /* cin >> end; */
    /* cout << "input: "; */
    /* cin >> input; */

    DataPoints dataPoints = generateDataPoints(start, end, count);
    cout << dataPoints;
    cout << "Lagrange: "<< lagrange(dataPoints, input) << " Expected: " << function(input) << endl;
    cout << "Spline: "<< spline(dataPoints, input, derivedFunction(start), derivedFunction(end)) << " Expected: " << function(input);
    return 0;
}
