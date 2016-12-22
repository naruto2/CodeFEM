// Smatrix A; 疎行列A
typedef SparseMatrix<double> Smatrix;


// Vector b; ベクトル
typedef VectorXd Vector;


// ...内部型
typedef Triplet<double> triplet;


// ...内部型
vector<triplet> tripletList;


// Tri(i,j,v); 疎行列へ値のセット
void Tri(int i, int j, double v) {
  tripletList.push_back(triplet(i,j, v));
}


// ...内部関数
void FinalizeTri(Smatrix &A) {
  A.setFromTriplets(tripletList.begin(), tripletList.end());
  tripletList.clear();
}


// Smatrix A = MapSmatrix(m,n); m行n列の疎行列を生成する
Smatrix MapSmatrix(int m, int n) {
  Smatrix A(m,n);
  FinalizeTri(A);
  return A;
}


// Vector b = MapVector(B,n); n次元ベクトルを配列Bより生成する
Vector MapVector(double *array, int n) {
  return Map<Vector>(array,n);
}


// ...内部関数
int Error(string msg) {
  cout << msg << endl;
  return -1;
}


// catSmatrix 疎行列のキャット
Smatrix catSmatrix(Smatrix A, Smatrix B) {
  Smatrix C;
  C.resize(A.rows(), A.cols()+B.cols());
  C.middleCols(0,A.cols()) = A;
  C.middleCols(A.cols(),B.cols()) = B;
  return C;
}


// vcatSmatrix 疎行列のキャット(縦方向)
Smatrix vcatSmatrix(Smatrix C, Smatrix D) {
  Smatrix CT = C.transpose();
  Smatrix DT = D.transpose();
  Smatrix E;

  E.resize(CT.rows(),CT.cols()+DT.cols());
  E.middleCols(0,CT.cols()) = CT;
  E.middleCols(CT.cols(),DT.cols()) = DT;

  return E.transpose();
}


// catVector ベクトルのキャット
Vector catVector(Vector x, Vector y) {
  Vector z(x.rows()+y.rows());
  z << x, y;
  return z;
}


// halvVector ベクトル長を半分に
Vector halfVector(Vector xx) {
  int i, n;
  n = xx.rows()/2;
  Vector x(n);
  for(i=0; i<n; i++) x(i) = xx(i);
  return x;
}


// IsSymmetric 行列が対称か否か
int IsSymmetric(Smatrix &A) {
  int n = A.cols();
  for(int i=0; i<n; i++)
    for(int j=i+1; j<n; j++)
      if ( A.coeffRef(i,j) != A.coeffRef(j,i) )
	return 0;
  return 1;
}




// norm(b); L2ノルム
double norm(const Vector &b) {
  return sqrt(b.dot(b));
}


// dot(a,b); ベクトルの内積
double dot(const Vector &a, const Vector &b) {
  return a.dot(b);
}


// 残差チェック
double check(Vector method(Smatrix &A,Vector &b), Smatrix &A, Vector &b) {
  Vector x = method(A,b);
  return norm(A*x - b);
}


// 時間計測
void
count( string str, Vector method(Smatrix &A,Vector &b), Smatrix &A, Vector &b) {

  clock_t   start = clock();
  double residual = check(method,A,b);
  clock_t     end = clock();

  printf("%f \t", (double)(end-start)/CLOCKS_PER_SEC);
  printf("%e \t", residual);
  printf("%s\n",str.c_str());
  return;
}


// Elu LU分解による連立一次方程式の求解
Vector Elu(Smatrix &A, Vector &b) {
  SparseLU<Smatrix> solver;
  solver.compute(A);
  return solver.solve(b);
}


// Ebicgstab BiCGSTAB法による連立一次方程式の求解
Vector Ebicgstab(Smatrix &A, Vector &b) {
  BiCGSTAB<Smatrix> solver;
  solver.compute(A);
  return solver.solve(b);
}


// Ecg  CG法による連立一次方程式の求解
Vector Ecg(Smatrix &A, Vector &b) {
  ConjugateGradient<Smatrix> solver;
  solver.compute(A);
  return solver.solve(b);
}


// Eogita 荻田氏の方法(対称行列化)
Vector Eogita(Smatrix &A, Vector &b) {

  Smatrix AT = A.transpose();

  Smatrix AA = vcatSmatrix( catSmatrix(A+AT, A-AT),
			    catSmatrix(AT-A,-A-AT));

  Vector  bb = catVector(2*b, -2*b);

  return halfVector(Ecg(AA,bb));
}
