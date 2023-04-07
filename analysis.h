#include <vector>
#include <armadillo>

arma::mat jack_knife (std::vector <arma::mat>& _d)
{     
        std::cout << _d.size() << "\n";
	arma::mat m_tmp (arma::size(_d[0]));
	for ( int i=0; i< _d.size(); ++i){
            m_tmp = m_tmp + _d[i];
	}
	return m_tmp/_d.size();
}
