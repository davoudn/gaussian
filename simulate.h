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
        std::vector<move_base<PDF>> m_moves;
        std::vector<meassure_base<PDF>> m_meassures;
        drawSample (t_pdf _pdf):m_pdf(_pdf), m_configs(1,arma::Col<double>(_pdf->size())) {}

        void meassure(){
               for (typename std::vector<meassure_base<PDF>>::iterator it = m_meassures.begin(); it!=m_meassures.end(); it++){
                       it->meassure();
                }
        return;
        }

        void move(){
                for (typename std::vector< move_base<PDF> >::iterator it = m_moves.begin(); it!=m_moves.end(); it++){
                       it->move();
                       if (it->m_acc == true ){
                          meassure ();
                          record_configs();
                       }
                }
                return;
        }

        void record_configs () {
                m_configs.push_back(m_pdf.m_x);
                return;
        }

        void do_it(int n_iter){
             for (int it=0;it < n_iter; it++) {
                 move();
             }
             return;
        }

        void add_move (move_base<PDF> _move){
              m_moves.push_back(_move);
        }

        void add_meassure ( meassure_base<PDF> _meassure){
              m_meassures.push_back (_meassure);
        }

        const arma::Col<double>& get_current_config () {
                return m_pdf->m_x;
        }
};

