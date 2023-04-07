/* ------------------------------------ */
#include <vector>
#include <fstream>
template <typename V>
struct histogram {};

template <>
struct histogram<single> {
        std::vector<double> m_v;
        std::vector<double> m_h;
        int mesh_no;
        double m_field_length;

        histogram (std::vector<double>& _v):m_v(_v), m_h(1,0)
        {}

        void make(int _mesh_no, double _field_length)
        {
                double l_delta = m_field_length/(_mesh_no+1);
                double l_norm{0};
                h.resize(_mesh_no,0);
                std::for_each (m_v.begin(), m_v.end(), [this,_field_length,_mesh_no,l_delta](double x) {
                                for(int i=0;i<_mesh_no; ++i)
                                   if ( x > ( i*l_delta - _field_length/2.0) && x < ((i+1)*l_delta - _field_length/2.0) ) {
                                      this->m_h[i]+=1.0;
                                      }
                                   } );
                for (std::vector<double>::iterator it=h.begin();it!=h.end();it++)
                        l_norm+= *it;

//              std::for_each ( h.begin(), h.end(), [this,_norm](double& x) { x=x/_norm; });
                for (std::vector<double>::iterator it=h.begin();it!=h.end();it++)
                        *it = *it/l_norm/l_delta;
                std::fstream fout;

                fout.open ("hist.dat",std::fstream::out);
                for (int i=0; i < _mesh_no; i++)
                        fout << i<< " "  << i*l_delta - _field_length/2.0 << " " << m_h[i] << "\n";
        }
};

