typedef std::vector< std::map< unsigned int, double> > cpumatrix;
typedef viennacl::compressed_matrix<double>            gpumatrix;
typedef viennacl::vector<double>                       gpuvector;


void matrix2gpumatrix(sparse::matrix<double> &A, viennacl::compressed_matrix<double> &Agpu)
{
  int n = A.size();
  FILE *fp = fopen("A.matrix","w");
  cpumatrix Acpu(n-1);
  printf("n-1 = %d\n",n-1);
  for (unsigned int i=1; i<A.size(); i++) {
    for ( auto it : A[i] ){
      int j = it.first;
      if (  A.size()-1 < j ) { printf("j=%d abort()\n",j); continue;}
      if (  j <= 0 ) {printf("j=%d abort()\n",j); abort();}
      Acpu[i-1][j-1] = it.second;
      fprintf(fp,"%d %d %f\n",i-1,j-1,it.second);
    }
  }
  fclose(fp);
  copy(Acpu, Agpu);
}

void fprintvector(vector<double>& x)
{
  FILE *fp = fopen("x.vector","w");
  for ( int k=0; k<x.size(); k++) fprintf(fp,"x[%d] = %f\n",k,x[k]);
  fclose(fp);
}


void vector2gpuvector(vector<double>& b, gpuvector& bgpu)
{
  for(int k=1; k<b.size(); k++) bgpu[k-1] = b[k];
}

void gpuvector2vector(gpuvector& xgpu, vector<double>& x)
{
  for(int k=1; k<x.size(); k++) x[k] = xgpu[k-1];
}

