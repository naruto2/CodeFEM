typedef std::vector< std::map< unsigned int, double> > cpumatrix;
typedef viennacl::compressed_matrix<double>            gpumatrix;
typedef viennacl::vector<double>                       gpuvector;


void matrix2gpumatrix(sparse::matrix<double> &A, viennacl::compressed_matrix<double> &Agpu)
{
  int n = A.size();

  cpumatrix Acpu(n);

  for (unsigned int i=0; i<A.size(); i++) {
    for ( auto it : A[i] ){
      int j = it.first;
      Acpu[i][j] = A[i][j];
    }
  }
  copy(Acpu, Agpu);
}

void vector2gpuvector(vector<double>& b, gpuvector& bgpu)
{
  copy(b.begin(), b.end(), bgpu.begin());
}

void gpuvector2vector(gpuvector& xgpu, vector<double>& x)
{
  copy(xgpu.begin(), xgpu.end(), x.begin());
}

