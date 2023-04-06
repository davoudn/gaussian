
#include "simulate.h"






struct simulateGaussian: public drawSample<gaussian> {
       using t_pdf = typename std::shared_ptr<gaussian>;

       simulateGaussian (t_pdf _pdf):drawSample<gaussian>(_pdf) {
                    add_move(move_local<gaussian>(_pdf));
                    add_meassure (meassure_2_point_correlation<gaussian>(_pdf));
       }
};




int main()
{

	std::shared_ptr<gaussian> g_pdf = std::make_shared<gaussian>(arma::ones(3,3),3);

	std::cout << g_pdf->energy() << "\n";

	std::cout << g_pdf->delta_energy(0,3.0) << "\n";

	std::cout << g_pdf->energy(arma::ones(3)) << "\n";
	g_pdf->set_config (arma::ones(3));
     
	std::cout << g_pdf->energy() << "\n";

	std::cout << g_pdf->delta_energy(0,0.1) << "\n";

	simulateGaussian sg_pdf(g_pdf);

	return 0;
}
