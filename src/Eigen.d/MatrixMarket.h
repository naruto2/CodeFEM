void gzip(string file) {
  string command;
  command = "gzip -f " + file;
  system(command.c_str());
}

void fprintmtx(string file, Smatrix &A) {
  FILE *fp;
  if ( file == "stdout" ) fp = stdout;
  else {
    fp = fopen(file.c_str(),"w");
    if( fp == NULL ) return;
  }
  if ( IsSymmetric(A) )
    fprintf(fp,"%%%%MatrixMarket matrix coordinate real symmetric\n");
  else
    fprintf(fp,"%%%%MatrixMarket matrix coordinate real general\n");
  
  int n = 0;
  for(int i=0;i<A.outerSize();++i)
    for(Smatrix::InnerIterator it(A,i);it;++it) n++;

  fprintf(fp,"%d %d %d\n", A.rows(), A.cols(), n);

  for(int i=0;i<A.outerSize();++i)
    for(Smatrix::InnerIterator it(A,i);it;++it)
      fprintf(fp,"%d %d %le\n",it.row()+1,it.col()+1,it.value());
  if ( file == "stdout" ) return;
  fclose(fp);
}


void fprintmtx(string file, Vector &b) {
  FILE *fp;
  if ( file == "stdout" ) fp = stdout;
  else {
    fp = fopen(file.c_str(),"w");
    if( fp == NULL ) return;
  }

  fprintf(fp,"%%%%MatrixMarket matrix array real general\n");

  fprintf(fp,"%d 1\n", b.size());
  
  int n = b.size();
  for(int i=0;i<n;++i)
    fprintf(fp,"%le\n",b[i]);

  if ( file == "stdout" ) return;
  fclose(fp);
}

void fprintmtx(string name, Smatrix &A, Vector &b) {

  string Afile, bfile;

  Afile = name + ".mtx";
  bfile = name + "_rhs1.mtx";

  fprintmtx(Afile,A);
  fprintmtx(bfile,b);

  gzip(Afile);
  gzip(bfile);
}

void fscanmtx(string file, Smatrix &A) {
  string in_str;
  ifstream ifs(file);

  while(1) {
    getline(ifs, in_str);
    if ( in_str.c_str()[0] != '%' ) break;
  }
  int rows, cols, n;
  sscanf(in_str.c_str(),"%d %d %d",&rows,&cols,&n);

  int i, row, col;
  double v;
  for ( i = 0; i<n; i++) {
    getline(ifs, in_str);
    sscanf(in_str.c_str(),"%d %d %le",&row, &col, &v);
    Tri(row-1,col-1,v);
  }
  A = MapSmatrix(rows,cols);
}  

void fscanmtx(string file, Vector &b) {
  string in_str;
  ifstream ifs(file);

  while(1) {
    getline(ifs, in_str);
    if ( in_str.c_str()[0] != '%' ) break;
  }
  int rows, n;
  sscanf(in_str.c_str(),"%d %d",&rows,&n);
  b.resize(rows);

  int i;
  double v;
  for ( i = 0; i<rows; i++) {
    getline(ifs, in_str);
    sscanf(in_str.c_str(),"%le",&v);
    b[i] = v;
  }
}  

int exists(string fname)
{
  FILE *file;
  if ((file = fopen(fname.c_str(), "r")))
    {
      fclose(file);
      return 1;
    }
    return 0;
}


void fscanmtx(string name, Smatrix &A, Vector &b) {
  string Afile = name + ".mtx";
  string bfile = name + "_rhs1.mtx";
  string command;
  command = "gunzip " + Afile + ".gz";
  system(command.c_str());
  command = "gunzip " + bfile + ".gz";
  system(command.c_str());

  fscanmtx(Afile,A);
  command = "gzip " + Afile;
  system(command.c_str());
  
  if ( exists( bfile ) ) {
    fscanmtx(bfile,b);
    command = "gzip " + bfile;
    system(command.c_str());
  }  
  else
    b = A*VectorXd::Constant(A.cols(),1);

}
