int halfbw(Matrix& a)
{
  int n, i, j, w ;

  n = a.size();
  for( i=0; i<n; i++ )
    for( auto it : a[i] ) {
      j = it.first;
      if ( a[i][j] != 0.0 ) if( w < abs(i-j) ) w = abs(i-j);
    }
  return w;
}


/*!
 * LU分解(ピボット交換なし)
 *  - 行列A(n×n)を下三角行列(L:Lower triangular matrix)と上三角行列(U:Upper triangular matrix)に分\\
解する
*  - L: i >= j,  U: i < j の要素が非ゼロでUの対角成分は1
*  - LとUを一つの行列にまとめた形で結果を返す
* @param[inout] A n×nの係数行列．LU分解した結果を格納する．
* @param[in] n 行列の大きさ
* @return 1:成功,0:失敗
*/
int le(int m, int n) {
  return ( m < n )? m : n;
}

int ge(int m, int n) {
  return ( m > n )? m : n;
}


int LUDecomp(matrix &A)
{
  int n = A.size();
  if(n <= 0) return 0;
  int w = halfbw(A);
  int k;

  
  for(int i = 0; i < n; ++i){
    printf("i = %d\n",i);
    // l_ijの計算(i >= j)
    for(int j = ge(0,i-w-1); j <= i; ++j){
      double lu = A[i][j];
      for(int k = ge(0,i-w-1); k < j; ++k){
	  lu -= A[i][k]*A[k][j];    // l_ik * u_kj
      }
      A[i][j] = lu;
      //if (A[i][j] != 0.0 ) printf("%f\n",A[i][j]);
    }
    // u_ijの計算(i < j)
    for(int j = i+1; j < le(i+w+1,n); ++j){
      double lu = A[i][j];

      //for(int k = 0; k < i; ++k){
      for ( auto it: A[i] ) {
	k = it.first;
	if ( k < i )
	  lu -= A[i][k]*A[k][j];    // l_ik * u_kj
      }
      A[i][j] = lu/A[i][i];
    }
  }

  return 1;
}


int LUSolver(matrix &A, vector<double> &b, vector<double> &x)
{
  int n = A.size();
  if(n <= 0) return 0;
  int j;
  // 前進代入(forward substitution)
  //  LY=bからYを計算
  for(int i = 0; i < n; ++i){
    double bly = b[i];
    //for(int j = 0; j < i; ++j){
    for ( auto it : A[i] ) {
      j = it.first;
      if ( j < i )
	bly -= A[i][j]*x[j];
    }
    x[i] = bly/A[i][i];
  }

  // 後退代入(back substitution)
  //  UX=YからXを計算
  for(int i = n-1; i >= 0; --i){
    double yux = x[i];
    //for(int j = i+1; j < n; ++j){
    for ( auto it : A[i] ) {
      j = it.first;
      if ( i < j && j < n )
	yux -= A[i][j]*x[j];
    }
    x[i] = yux;
  }

  return 1;
}
