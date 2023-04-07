#include <iostream>
#include <vector>
#include <memory>
#include <armadillo>



struct a {
	std::shared_ptr <std::vector<double>> m_v;
	a(std::shared_ptr <std::vector<double>>& _v):m_v(_v){}
	void set(int n, double val) {(*m_v)[n]=val;}
	virtual void  print () {std::cout << "this is parent print \n ";}
};

struct ca : public a {
ca(std::shared_ptr <std::vector<double>>& _v):a(_v) {}

void print () {std::cout << "this is child print  \n ";}

};


struct b {
        std::shared_ptr <std::vector<double>> m_v;
        b(std::shared_ptr <std::vector<double>>& _v):m_v(_v){}
        void set(int n, double val) {(*m_v)[n]=val;}
};


int main ()
{
	std::shared_ptr <std::vector<double>> v =  std::make_shared<std::vector<double>>(10,0.0);
        
	a v_a(v); 
	b v_b(v);
	int n_col, n_row, n_slice;
	arma::Mat<double> M(5,5);
	arma::Col<double> V(5);
        arma::Mat<double> s ( V.t() * M * V);
	std::cout << (V.t() * M * V)<< "\n";
	std::cout << arma::size(M) << std::endl;
	std::shared_ptr < a >  pa =  std::make_shared< a >(v);
        pa->set (6,30.0);
        std::cout << 6 << " " <<  (*v)[6] << "\n";	
	v_a.set(5,10);
	std::cout << (*v)[5] << "\n";
	v_b.set(5,20);
	std::cout << (*v)[5] << "\n";

	a v_ca (v);
	ca d(v);
//	v_ca.print ();
	return 0;

}
