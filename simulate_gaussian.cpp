#include "simulate.h"
#include "analysis.h"
struct simulateGaussian_metropolice: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian_metropolice (t_pdf _pdf):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_local<gaussian>(_pdf));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};

struct simulateGaussian_langevin: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian_langevin (t_pdf _pdf, double _D, double _dt):drawSample<gaussian>(_pdf) {
                    this->add_move(new langevin<gaussian>(_pdf,_D,_dt));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};


struct simulateGaussian_eig_decompose: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian_eig_decompose (t_pdf _pdf, int _nchanged):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_eig_decompose_gaussian(_pdf,_nchanged));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};


struct simulateGaussian_chol_decompose: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian_chol_decompose (t_pdf _pdf, int _nchanged):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_chol_decompose_gaussian(_pdf,_nchanged));
                    this->add_meassure (new meassure_2_point_correlation<gaussian>(_pdf));
       }

       void analyse () {
            std::cout << jack_knife (this->m_meassures[0]->m_results) << "\n";
       }
};

struct simulateGaussian_stochastic_relaxation: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian_stochastic_relaxation (t_pdf _pdf, double _omega, double _dt ):drawSample<gaussian>(_pdf) {
                    this->add_move(new move_stochastic_relaxation(_pdf, _omega, _dt));
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


/*
        simulateGaussian_langevin sg_langevin(g_pdf,1.0,0.05);
        sg_langevin.do_it (100000);
        sg_langevin.analyse ();
*/
	/*

        simulateGaussian_eig_decompose sg_eig_decompose(g_pdf,1);
        sg_eig_decompose.do_it (10000);
        sg_eig_decompose.analyse ();

*/
        simulateGaussian_stochastic_relaxation sg_stochastic_relaxation(g_pdf, 1.0, 1.0);
        sg_stochastic_relaxation.do_it (100000);
        sg_stochastic_relaxation.analyse ();



	return 0;
}
