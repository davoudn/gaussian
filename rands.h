#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <armadillo>
#include <map>
#include <memory>

struct single{};
struct multi{};

template <typename T>
struct unif;

template <>
struct unif<double> {
        double m_a1,m_a2;
        unif (double _a1, double _a2):m_a1(_a1), m_a2(_a2) {
                std::srand(std::time(nullptr));
        }
        double operator ()(){
                return ((double)std::rand()/RAND_MAX)*(m_a2-m_a1);
        }
};

template <>
struct unif<int> {
        double m_a1,m_a2;
        unif (int _a1, int _a2) {
                std::srand(std::time(nullptr));
        }
        double operator () (){
                int dummy = ((double)std::rand()/RAND_MAX)*(m_a2-m_a1);
                return dummy;
        }
};

template <typename V>
struct drawGaussian;

template <>
struct drawGaussian<single> {
        drawGaussian()
        {
            std::srand(std::time(nullptr)); // use current time as seed for random generator
        }

        double box_muller(){
                double u1{0.0}, u2{0.0}, eps{0.0001};
                while (u1 < eps || u2 < eps){
                        u1 = (double)std::rand()/RAND_MAX;
                        u2 = (double)std::rand()/RAND_MAX;
                }
                double x = std::pow(-2 * std::log(u1), 0.5) * std::cos(2 * 3.1415 * u2);
                return x;
                }
};

template <>
struct drawGaussian<multi> {
       private:
         drawGaussian<single> m_grand;
         int m_dim;
         arma::Col<double> m_eta;
        public:
         drawGaussian(int _dim):m_eta(_dim), m_dim(_dim),m_grand(){}

     const arma::Col<double>& operator ()() {
               std::for_each(m_eta.begin(), m_eta.end(), [this](double& x) {x=this->m_grand.box_muller();});
               return m_eta;
       }
};


