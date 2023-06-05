#include <cmath>
#include <vector>
#include <algorithm>
#include <armadillo>


template <typename T>
struct pdf {
        std::vector <arma::Col<T>> m_configs;
//	arma::Col<T> m_x;
        pdf(){}
        virtual T energy(){ return T(1);}
        virtual const arma::Col<T>& force (){ std::cout << "This is parrent pdf force (). No force!";}
        virtual T delta_energy(int n,T _xnew){}
        
        virtual int size (){ return 0;}
        virtual T ratio(int n, T _xnew){}
	virtual void update(int n, T _xnew) {}
        void record (){}
	void set_config () {}
};

template <typename T>
struct gaussian:public pdf<T> {
        arma::Mat<T> m_Q, m_U, m_stupid_fix,m_R, m_Q_inv;
        arma::Col<T> m_lambda, m_x, m_force;
        int m_dim;
        gaussian (arma::Mat<T> _Q, int _dim):m_Q(_Q),m_U(_dim,_dim), m_R(_dim,_dim),m_force(_dim), m_Q_inv(_dim,_dim)
                                                 ,m_lambda(_dim),m_x(arma::randu(_dim)),m_stupid_fix(1,1)
                                                 ,m_dim(_dim){}

      //  gaussian (arma::Mat<double>& _Q, int _dim):m_Q(_Q),m_U(_dim,_dim)
      //                                           ,m_lambda(_dim),m_x(_dim),m_stupid_fix(1,1)
                                                 //,m_dim(_dim){}


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

        T delta_energy(int n, T _xnew){
          T delta { _xnew - m_x[n]}, dummy {0.0};
          for (int i=0;i<m_dim; i++) {
		  dummy+= m_Q(n,i)*m_x(i)*delta;
          }
          dummy+= 0.5 * std::pow(delta,2.0)*m_Q(n,n);
          return dummy;
        }


         T ratio(int n, T _xnew){
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

	void update(int n, T _xnew) 
	{
             m_x[n] = _xnew; 
	}


};
template <typename T>
struct Op {
	std::vector< std::pair< <T> > m_elements;
	T g;
	arma::Mat<T> m_Mat;
	Op(int _N):m_N(_N), m_Mat(_N,_N)
        arma::Mat<T>& get_matrix()
	{
          
           std::for_each (m_elements.begin(), m_elements.end(), [this,m_Mat](double x) {
			                m_Mat(x.first,x.second) = 1.0;
                                      }
                                   );
          return m_Mat
	}
};




template <typename T>
struct fermion:public pdf<T> {
       Bmatrices<T> m_Bmatrices;
       std::vector<Op<T>> m_V;
       Op<T> m_T;
       arma::Mat<T> m_x;
       int m_Ltrot, m_Nsites;
       fermion (int _Ltrot, int _Nsites, int _Nfields):m_Ltrot(_Ltrot),m_Nsites(_Nsites),m_Nfields(_Nfields),m_Bmatrices(_Ltrot,_Nsites),m_x(_L,_Nfields),m_T(_Nsites), m_V(_Nfields, Op<T>(_Nsites))
       {}

       arma::Mat<T> exp (arma::Mat<T>& mat){
           return;
       }
       arma::Mat<T> get_Vl(int l)
       {
	       arma::Mat<T> vtmp(_N,_N)
	       for (int nfield=0; nfield<m_Nfields; nfield++){
                   vtmp = vtmp + m_x(l,nfield) * m_V[nfiel].g * m_V[nfield].get_matrix();
	       }
	       return vtmp;
       }

       void make_Bmatrices ()
       {
	   arma::Mat<T> expT (exp (-m_T));
           for (int l=0; l < m_Ltrot; l++){
              m_Bmatrices.B[l] = expT * exp (-get_Vl(l));
	   }
	   return;
       }


};

template <typename T>
struct Bmatrices {
        std::vector< arma::Mat<T> > m_B, m_Binv, m_Gr
        arma::Mat<T>  m_Ble, m_Blt, m_Gr;
        int m_L, m_N, m_l;
        Bmatrices(int _L, int _N):m_l(0), m_L(_L), m_N(_N), m_B(_L, arma::Mat<T>(_N, _N)), m_Binv(_L, arma::Mat<t_complex>(_N, _N)), m_Gr(_N, _N),  m_Ble(_N, _N),  m_Blr(_N, _N)
        {}

	void set_m_l(int l){
		m_l = l;
		return;
	}
        T ratio (int n, int l, T d)
        {
           propagate (l);
           return 1.0 + d * (m_Blr * m_Ble * m_Gr)(n,n);
        }

        void propagate (int l){
             m_l = l;
        }

        void update (int n, int m, int l, T d)
        {

        }
};

template <typename T>
struct VanilaHubbard:public pdf<T> {
        T m_U,m_t;
        double beta,dtau;
        fermion<T> m_fermion;
        boson<T>   m_boson;

        VanilaHubbard(T _U, T _t, T _beta, int _Ltrot, int _Nsites,int _Nfields):m_fermoin<T>(_Ltrot,_Nsites,_Nfields), m_boson<T>(_Ltrot,_Nsites,_Nfields),
                                                                                 m_U(_U),m_t(_t),m_beta(_beta)
        {
           make_V();
           make_T();
           m_fermions.make_Bmatrices();
        }
   
         void make_V() 
         {

         }
         void make_T()
         {

         }



};