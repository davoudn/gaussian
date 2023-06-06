#include <cmath>
#include <vector>
#include <algorithm>
#include <armadillo>

struct R {
   int m_Nx,m_Ny;
   int m_x,m_y;
   bool m_is_end;
   R(int _Nx,int _Ny):m_Nx(_Nx), m_Ny(_Ny), m_x(0), m_y(0), m_is_end(0)
   {}

   void increment()
   {
     m_x++;
     if ( m_x == m_Nx){
        m_x=0;
        m_y++;
        if( m_y == m_Ny){
            m_Ny=0;
            m_is_end = 1;
        }
          
     }
        
   }
   int idx() {
    return m_Nx * m_y + m_x;
   }

   bool is_end()
   {
     return m_is_end;
   }
};

template <typename LATTYPE>
struct lattice{};

struct square{};

template <>
struct lattice<square>{
   arma::Col<double> m_a1, m_a2, m_b1, m_b2;
   int m_Ldim1,m_Ldim2, m_Nsites;
   double m_d;
   std::map< std::pair<int, int>, int> m_point_to_index_map;
   std::map< int, std::pair<int, int>> m_index_to_point_map;
   std::vector < std::vector<int> > m_1st_nn_list; // nn stands for nearest neighbor
 //  lattice(arma::Col<double> _a1, arma::Col<double> _a2, int _Ldim1,int _Ldim2):m_d(1.01), m_Ldim1(_Ldim1), m_Ldim2(_Ldim2),m_Nsites(_Ldim1*_Ldim2),
 //                                                                               m_a2(_a2), m_a1(_a1)
 //  {
 //  }
   lattice(arma::Col<double> _a1, arma::Col<double> _a2, int _Ldim1,int _Ldim2):m_d(1.01), m_Ldim1(_Ldim1), m_Ldim2(_Ldim2),m_Nsites(_Ldim1*_Ldim2),
                                                                                m_a2(2), m_a1(2)
   {
    m_a1(0) = 0.0; m_a1(1) = 1.0;
    m_a2(0) = 1.0; m_a2(1) = 0.0;
   }
   void find_1st_distance()
   {
      return;
   }

   void make_1st_nn_list()
   {
    for (R r1(m_Ldim1,m_Ldim2); r1.is_end()!=true; r1.increment())
        for (R r2(m_Ldim1,m_Ldim2); r2.is_end()!=true; r2.increment()){
            if ( arma::norm((r1.m_x-r2.m_x) * m_a1 + (r1.m_y-r2.m_y) * m_a2 )  <= m_d ){
                m_1st_nn_list[r1.idx()].push_back(r2.idx());
            }
        }       
        return;  
   }

   const std::vector<int>& get_1st_nn_list(int idx)
   {
      return m_1st_nn_list[idx];
   }
   void make_maps()
   {
       for (R r(m_Ldim1,m_Ldim2); r.is_end()!=true; r.increment()){
            m_index_to_point_map[r.idx()] = std::make_pair(r.m_y,r.m_x);
            m_point_to_index_map[ std::make_pair(r.m_y,r.m_x)] = r.idx();
       }
       return;
   }

   std::pair<int,int> get_point(int idx)
   {
    return m_index_to_point_map[idx];
   }
   int get_idx(std::pair<int,int> point)
   {
    return m_point_to_index_map[point];
   }
   int get_Nsites()
   {
    return m_Nsites;
   }
};