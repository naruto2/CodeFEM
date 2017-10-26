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
  printf("Acpu.size() = %d\n",Acpu.size());
  printf("Acpu[0][0] = %f\n",Acpu[0][0]);
  printf("Acpu[1][1] = %f\n",Acpu[1][1]);
  printf("Acpu[A.size()-3][A.size()-3] = %f\n",Acpu[A.size()-3][A.size()-3]);
  printf("Acpu[A.size()-2][A.size()-2] = %f\n",Acpu[A.size()-2][A.size()-2]);
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
  //fprintvector(b);
  //exit(0);
  printf("bgpu.size() = %d\n",bgpu.size());
  for(int k=1; k<b.size(); k++) bgpu[k-1] = b[k];
  printf("b[0] = %f\n",b[0]);
  printf("b[1] = %f\n",b[1]);
  printf("b[b.size()-3] = %f\n",b[b.size()-3]);
  printf("b[b.size()-2] = %f\n",b[b.size()-2]);
  printf("b[b.size()-1] = %f\n",b[b.size()-1]);
}

void gpuvector2vector(gpuvector& xgpu, vector<double>& x)
{
  for(int k=1; k<x.size(); k++) x[k] = xgpu[k-1];
  printf("x[0] = %f\n",x[0]);
  printf("x[1] = %f\n",x[1]);
  printf("x[x.size()-3] = %f\n",x[x.size()-3]);
  printf("x[x.size()-2] = %f\n",x[x.size()-2]);
  printf("x[x.size()-1] = %f\n",x[x.size()-1]);
  fprintvector(x);
}

