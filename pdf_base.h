#include <cmath>
#include <vector>
#include <algorithm>
#include <armadillo>
#include "lattice.h"

template <typename T>
struct pdf {
        std::vector <arma::Col<T>> m_configs;
//	arma::Col<T> m_x;
        pdf(){}
        virtual T energy(){ return T(1);}
        virtual const arma::Col<T>& force (){ std::cout << "This is parrent pdf force (). No force!";}
        virtual T delta_energy(int l ,n ,T _xnew){}
        
        virtual int size (){ return 0;}
        virtual T ratio(int l, int n, T _xnew){}
	virtual void update(int l, int n, T _xnew) {}
        void record (){}
	void set_config () {}
};

template <typename T>
struct gaussian {
        arma::Mat<T> m_Q, m_U, m_stupid_fix,m_R, m_Q_inv;
        arma::Col<T> m_lambda, m_x, m_force;
        int m_dim, m_Ltrot, m_Nfields;
        gaussian (arma::Mat<T> _Q, int _Ltrot, int _Nfields):m_Q(_Q),m_dim(_Ltrot*_Nfields),m_U(_Ltrot*_Nfields,_Ltrot*_Nfields), 
                                                             m_R(_Ltrot*_Nfields,_Ltrot*_Nfields),m_force(_Ltrot*_Nfields), 
                                                             m_Q_inv(_Ltrot*_Nfields,_Ltrot*_Nfields), m_lambda(_Ltrot*_Nfields),
                                                             m_x(arma::Col<T>(_Ltrot*_Nfields)),m_stupid_fix(1,1) ,m_Ltrot(_Ltrot),m_Nfields(_Nfields){}

        
        gaussian                  (int _Ltrot, int _Nfields):m_Q(_Ltrot*_Nfields,_Ltrot*_Nfields),m_dim(_Ltrot*_Nfields),m_U(_Ltrot*_Nfields,_Ltrot*_Nfields), 
                                                             m_R(_Ltrot*_Nfields,_Ltrot*_Nfields),m_force(_Ltrot*_Nfields), 
                                                             m_Q_inv(_Ltrot*_Nfields,_Ltrot*_Nfields), m_lambda(_Ltrot*_Nfields),
                                                             m_x(arma::Mat<T>(_Ltrot*_Nfields)),m_stupid_fix(1,1) ,m_Ltrot(_Ltrot),m_Nfields(_Nfields){}
        T energy() {
	  m_stupid_fix = m_x.t() * ( m_Q * m_x );
          return m_stupid_fix(0,0)*0.5;
//        return 0.0;
        }

        T energy(arma::Col<T>& _x) {
          m_stupid_fix = _x.t() * ( m_Q * _x );
          return m_stupid_fix(0,0);
        }
 
        T energy(arma::Col<T> _x) {
          m_stupid_fix = _x.t() * ( m_Q * _x );
          return m_stupid_fix(0,0)*0.5;
        }
 
        arma::Col<T>& force (){
		m_force = -m_Q * m_x;
	//	std::cout << m_force << "ff\n";
          return m_force;
        }

        T delta_energy(int l, int n, T _xnew){
          T delta { _xnew - m_x(l,n)}, dummy {0.0};
          for (int i=0;i<m_dim; i++) {
		  dummy+= m_Q(n,i)*m_x(i)*delta;
          }
          dummy+= 0.5 * std::pow(delta,2.0)*m_Q(n,n);
          return dummy;
        }


         T ratio(int l, int n, T _xnew){
          T delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<m_dim; i++){
               dummy+= m_Q(n,i)*m_x(i)*delta;
	  }
          dummy+= 0.5 * std::pow(delta,2.0)*m_Q(n,n);
          return std::exp(-dummy);
        }

        void set_config (arma::Col<T> _x){
		m_x = _x; 
	}

        int size (){return  m_dim;}

	const arma::Mat<T>& get_covariance ()
	{
		return m_Q;
	}
        void set_covariance_inv ( )
	{
		m_Q_inv = m_Q.i();
		return;
	}

        arma::Mat<T>& get_covariance_inv ()
        {
                return m_Q_inv;
        }

	
	void set_eig_decompose()
	{
		eig_sym(m_lambda, m_U, m_Q);
		std :: cout << m_lambda << "\n";
	}

	void set_chol_decompose()
        {
		m_R = arma::chol(m_Q);
                std :: cout << m_R << "\n";
        }

	void update(int l, int n, T _xnew) 
	{
             m_x(l,n) = _xnew; 
	}


};
template <typename T>
struct Op {
	std::vector< std::pair< <T> > m_elements;
        std::map <std::pair< <T>, T> m_values;
	T g;
	arma::Mat<T> m_Mat;
	Op(int _N):m_N(_N), m_Mat(_N,_N)
        {}

        arma::Mat<T>& get_matrix()
	{
           std::for_each (m_elements.begin(), m_elements.end(), [this,m_Mat](double x) {
			                this->m_Mat(x.first,x.second) = this->m_values[x];
                                        });
          return m_Mat
	}
	void set_elements(int x, int y, T value)
	{
                 m_elements.push_back(std::make_pair(x, y));
                 m_value(std::make_pair(x, y)) = cone;
	}
	void set_g(T val)
	{
	   g = cone * m_dtau;
	}
};

template <typename T>
struct Bmatrices {

        private:
        std::vector< arma::Mat<T> > m_B, m_Binv;
        arma::Mat<T>  m_Bless, m_Bgreater, m_G, m_I;
        int m_L, m_N, m_l;
        Bmatrices(int _L, int _N):m_l(0), m_L(_L), m_N(_N), m_B(_L, arma::Mat<T>(_N, _N)), m_Binv(_L, arma::Mat<t_complex>(_N, _N)), m_G(_N, _N),  m_Bless(_N, _N),  m_Bgreater(_N, _N),I(_N, _N, fill::eye);
        {}

	void set_m_l(int l){
             m_l = l;
	     return;
	}

        void propagate (int l){
             if ( m_l == l ){
	        return;
	     } else if ( l > m_l ) {
                    for(int i=m_l+1; i<= l; i++ ){	
			  m_G = Binv[i] * m_G * B[i];   
			  Bless = Bless * B[i];
			  Bgreater = Binv[i] * Bgreater;
                    }
                    m_l = l;
	     } else if (l < m_l) {
                    for(int i=m_l-1; i>= l; i-- ){
                          m_G = B[i] * m_G * Binv[i];
                          Bless = Bless * Binv;
                          Bgreater = B[i] * Bgreater;
                    }
                    m_l = l;
	     }
	     return;
        }

        public:
        T ratio (int n, int m, int l, T d)
        {
             propagate (l);
             return 1.0 + d * ( m_I - m_G)(n,n);
        }

        void update (int n, int m, int l, T d)
        {
             propagate(l);
             arma::Mat<T> tmp(m_N,m_N);
	     tmp(n,m) = d;
	     m_G = m_G + (m_G * m_Bgreater) * tmp * (m_Bless * m_G) / (cone + d);
	     m_B[l] *= tmp;
	     m_Binv[l] = inv(m_B[l]);
        }
	 
};


template <typename T>
struct fermion:public pdf<T> {
       Bmatrices<T> m_Bmatrices;
       std::vector<Op<T>> m_V;
       Op m_T;
       arma::Mat<T> m_x;
       int m_Ltrot, m_Nsites;
       fermion (int _Ltrot, int _Nsites, int _Nfields):m_Ltrot(_Ltrot),m_Nsites(_Nsites),m_Nfields(_Nfields),m_Bmatrices(_Ltrot,_Nsites),m_x(_L,_Nfields),
                                                       m_T(_Nsites), m_V(_Nfields, Op<T>(_Nsites))
       {}

       arma::Mat<T> exp (arma::Mat<T>& mat)
       {
              return;
       }
       arma::Mat<T> get_Vl(int l)
       {
	      arma::Mat<T> vtmp(_N,_N)
	      for (int nfield=0; nfield<m_Nfields; nfield++){
                  vtmp = vtmp + m_x(l,nfield) * m_V[nfield].g * m_V[nfield].get_matrix();
	      }
	      return vtmp;
       }

       void init_Bmatrices ()
       {
              for (int l=0; l < m_Ltrot; l++){
                  m_Bmatrices.B[l] = exp (-m_T) * exp (-get_Vl(l));
	      }
	      return;
       }

       T ratio (int l, int nfield, T _xnew)
       {
	       T _ratio{ T(0) };
	       Bmatrices _Bmatrices {m_Bmatrices};
               arma::Mat<T> _Delta { exp( (_xnew - m_x[nfield]) * m_V[nfield].g * m_V[nfield].get_matrix() ) };
	       for (int i=0; i< m_Nsites; ++i)
		   _Delta(i,i) -= T(1);

	    //  //
	       for ( auto it = m_V[i].m_elements.begin(); it != m_V[i].m_elements.end(); ++it ){
		        _ratio *= _Bmatrices.ratio( it->first , it->second, l, _Delta(it->first , it->second));     	    
                        _Bmatrices.update (it->first , it->second, l, _Delta(it->first , it->second));
                    }
	       return _ratio;
       }

       void update (int l, int nfield, T _xnew)
       {
	       arma::Mat<T> _Delta { exp( (_xnew - m_x[nfield]) * m_V[nfield].g * m_V[nfield].get_matrix() ) };

               for ( auto it = m_V[i].m_elements.begin(); it != m_V[i].m_elements.end(); ++it ){
                   m_Bmatrices.update (it->first , it->second, l, _Delta(it->first , it->second));
               }
	       m_x(l,nfield) = _xnew;
	       return;
       }

};

template <typename BOSONTYPE, typename T>
struct boson:public pdf<T> {};

template <typename T>
  struct boson<gaussian>:public pdf <T>{
       int m_Ltrot, m_Nfields;
       boson(int _Ltrot, int _Nfields):gaussian(_Ltrot,_Nfields),m_Ltrot(_Ltrot), m_Nfields(_Nfields)
       {}
};

template <typename T>
struct SquareHubbard:public pdf<T> {
        T m_U,m_t,m_mu;
        double beta,dtau;
        fermion<T> m_fermion;
        boson<T>   m_boson;
        int m_Lattdim, m_Nsites;
        lattice<Square> m_lattice;
        SquareHubbard(T _U, T _t, T _beta, T _mu, int _Ltrot, int _Lattdim,int _Nfields):m_fermoin<T>(_Ltrot,_Lattdim*_Lattdim,_Nfields),
                                                                                  m_boson<T>(_Ltrot,_Lattdim*_Lattdim,_Nfields), m_lattice(_Lattdim,_Lattdim),
                                                                                  m_U(_U),m_t(_t),m_beta(_beta),m_mu(_mu),m_Lattdim(_Lattdim),
                                                                                  m_Nsites(_Lattdim*_Lattdim),m_dtau(_beta/_Ltrot)
        {
           make_TV();
           make_Q();           
        }
   
        void make_TV() 
        {
            int nsite{0};
            for (int nfield=0; nfield < m_fermion.m_Nfields; nfield++){
                 m_fermion.m_V[nfield].set_elements(nsite,nsite,cone);
                 m_fermion.m_V[nfield].set_g = cone * m_dtau;
                 nsite++;
            }
            
            // hubbard model in square lattice with nearest neigghbor hopping.
            for (int i=0; i < m_lattice.get_Nsites(); i++){
                for (int j=0; j< m_lattice.get_1st_nn_list(i).size(); j++){
                    fermion.m_T.set_elements(i,j,-m_t);
                }
                fermion.m_T.set_elements(i,i,-m_u);
            }
	    fermion.m_T.g = cone * m_dtau;
            m_fermions.init_Bmatrices();
            return;
        }
        
        void make_Q() 
        {
            for (int nfield=0; nfield < m_fermion.m_Nfields; nfield++){
                nsite++;
            }
                       
            return;
        }
        
        T ratio(int l, int n, T _xnew){
          return m_fermion.ratio(l,n,_xnew) * m_boson.ratio(l,n,_xnew);
        }
        void update(int l, int n, T _xnew){
          m_fermion.update(l,n,_xnew);
          m_boson.update(l,n,_xnew);
          return;
        }

};
