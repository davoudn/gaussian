#include <cstdlib>
#include <fstream>
#include <ctime>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <armadillo>
#include <map>
#include <memory>
#include "rands.h"
#include "pdf_base.h"
#include "meassure_base.h"
#include "move_base.h"

template <typename PDF>
struct drawSample
{
        using t_pdf = typename std::shared_ptr<PDF>;
        t_pdf m_pdf;
        std::vector <arma::Col<double>> m_configs;
        std::vector<move_base<PDF>*> m_moves;
        std::vector<meassure_base<PDF>*> m_meassures;
        drawSample (t_pdf _pdf):m_pdf(_pdf), m_configs(1,arma::Col<double>(_pdf->size())) {}

        void meassure(){
               for (int it=0; it < m_meassures.size(); it++)
                       m_meassures[it]->meassure();
        return;
        }

        void move(){
                for (int it=0; it < m_moves.size(); it++){
                       m_moves[it]->move();
                       if (m_moves[it]->m_acc == true ){
                       //   meassure ();
                       }
		       meassure ();
                }
                return;
        }

        void record_configs () {
                m_configs.push_back(m_pdf->m_x);
                return;
        }

        void do_it(int n_iter){
             for (int it=0;it < n_iter; it++) {
                 move();
//		 std::cout << m_pdf->m_x << "\n";
             }
             return;
        }

        void add_move (move_base<PDF>* _move){
              m_moves.push_back(_move);
	      std::cout << "_move added. \n";
        }

        void add_meassure ( meassure_base<PDF>* _meassure){
              m_meassures.push_back (_meassure);
	      std::cout << "_meassure added. \n"; 
        }

        const arma::Col<double>& get_current_config () {
                return m_pdf->m_x;
        }
};

template <typename MOVETYPE, typename PDF>
struct simulate;


