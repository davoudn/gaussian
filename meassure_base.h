
#include <vector>
#include <armadillo>
#include <memory>

template <typename PDF>
struct meassure_base{
 using t_pdf = typename std::shared_ptr<PDF>;
 t_pdf m_pdf;
 meassure_base(t_pdf& _pdf ):m_pdf(_pdf){}


};

template <typename PDF>
struct meassure_2_point_correlation: public meassure_base <PDF>{
        using t_pdf = typename std::shared_ptr<PDF>;

        std::vector<arma::Mat<double>> m_results;
        arma::Mat<double> m_correlations;
        int m_size;
        meassure_2_point_correlation (t_pdf& _pdf):meassure_base<PDF>(_pdf),m_correlations(_pdf->size(),_pdf->size()),
                                                   m_size(_pdf->size()) {}

        void meassure () {
             for (int i=0; i<m_size; i++)
                     for (int j=0; j<m_size; j++){
                             m_correlations(i,j) = this->m_pdf->m_x(i) * this->m_pdf->m_x(j);
                     }
             m_results.push_back(m_correlations);
             return;
        }
};

