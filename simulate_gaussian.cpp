#include "simulate.h"
#include "analysis.h"


template <>
struct simulate<local,gaussian>: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulate (t_pdf _pdf):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_type<local,gaussian>(_pdf));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};
/*struct simulateGaussian_langevin: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian_langevin (t_pdf _pdf, double _D, double _dt):drawSample<gaussian>(_pdf) {
                    this->add_move(new langevin<gaussian>(_pdf,_D,_dt));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};
*/
template <>
struct simulate<langevin,gaussian>: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulate (t_pdf _pdf, double _D, double _dt):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_type<langevin,gaussian>(_pdf,_D,_dt));
//                  this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};


template <>
struct simulate<eigdecompose,gaussian>: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulate (t_pdf _pdf, int _nchanged):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_type<eigdecompose,gaussian>(_pdf,_nchanged));
//                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};

template <>
struct simulate<choldecompose,gaussian>: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulate<choldecompose,gaussian> (t_pdf _pdf, int _nchanged):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_type<choldecompose,gaussian>(_pdf,_nchanged));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};

template<>
struct simulate<stochasticrelaxation,gaussian>: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulate (t_pdf _pdf, double _omega, double _dt ):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_type<stochasticrelaxation,gaussian>(_pdf, _omega, _dt));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};



int main()
{
	int N {5};
	double delta {0.5};
        arma::mat _Q (N,N);
	for (int i=0; i <N ; i++)
	    for (int j=0; j <N ; j++)
		    if (i==j) _Q(i,j) = (1.0);
		    else _Q(i,j) = delta;

	std::shared_ptr<gaussian> g_pdf = std::make_shared<gaussian>(_Q,N);

	g_pdf->set_covariance_inv();
 	std::cout << g_pdf->get_covariance_inv() << "\n";
//	g_pdf->get_eig_decompose();
/*	
	simulateGaussian_metropolice sg_metropolice(g_pdf);
	sg_metropolice.do_it (10000);
        sg_metropolice.analyse ();

*/

        simulate<langevin,gaussian> sg_langevin(g_pdf,1.0,0.05);
        sg_langevin.do_it (10000);
        sg_langevin.analyse ();

	
/*
        simulateGaussian_eig_decompose sg_eig_decompose(g_pdf,5);
        sg_eig_decompose.do_it (1000);
        sg_eig_decompose.analyse ();
*/

/*        simulateGaussian_stochastic_relaxation sg_stochastic_relaxation(g_pdf, 1.0, 1.0);
        sg_stochastic_relaxation.do_it (5000);
        sg_stochastic_relaxation.analyse ();
*/
	return 0;
}
