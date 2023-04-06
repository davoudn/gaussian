/* ------------------------------------ */
template <typename V>
struct histogram {};

template <>
struct histogram<single> {
        std::vector<double> v;
        std::vector<double> h;
        int mesh_no;
        double field_length;

        histogram (std::vector<double>& _v):v(_v), h(1,0)
        {}

        void make(int _mesh_no, double _field_length)
        {
                double _delta = _field_length/(_mesh_no+1);
                double _norm{0};
                h.resize(_mesh_no,0);
                std::for_each (v.begin(), v.end(), [this,_field_length,_mesh_no,_delta](double x) {
                                for(int i=0;i<_mesh_no; ++i)
                                   if ( x > ( i*_delta - _field_length/2.0) && x < ((i+1)*_delta - _field_length/2.0) ) {
                                      this->h[i]+=1.0;
                                      }
                                   } );
                for (std::vector<double>::iterator it=h.begin();it!=h.end();it++)
                        _norm+= *it;

//              std::for_each ( h.begin(), h.end(), [this,_norm](double& x) { x=x/_norm; });
                for (std::vector<double>::iterator it=h.begin();it!=h.end();it++)
                        *it = *it/_norm;
                std::fstream fout;

                fout.open ("hist.dat",std::fstream::out);
                for (int i=0; i < _mesh_no; i++)
                        fout << i<< " "  << i*_delta - _field_length/2.0 << " " << h[i] << "\n";
        }
};

