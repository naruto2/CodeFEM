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
