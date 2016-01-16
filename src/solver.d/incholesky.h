/*!
 * 不完全コレスキー分解(incomplete Cholesky decomposition)
 *  - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
 *  - l_ii = 1とした場合
 *  - L: i > jの要素が非ゼロで対角成分は1
 *  - 行列Aの値が0である要素に対応する部分を飛ばす
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[out] d 対角行列(対角成分のみ)
 * @param[in] n 行列の大きさ
 * @return 1:成功,0:失敗
 */
int IncompleteCholeskyDecomp(matrix &A, matrix  &L, vector<double> &d, int n)
{
  if(n <= 0) return 0;

  d[0] = A[0][0];
  L[0][0] = 1.0;

  for(int i = 1; i < n; ++i){
    // i < k の場合
    internal_map Ai = A[i];
    for ( internal_map::iterator j = Ai.begin(); j != Ai.end() && j->first < i; j++ ) {
      if(fabs(A[i][j->first]) < 1.0e-10) continue;

      double lld = A[i][j->first];

      L[i][0]; // L.sync();
      internal_map Li = L[i];
      for( internal_map::iterator k = Li.begin(); k != Li.end() && k->first < j->first; k++)
	lld -= L[i][k->first]*L[j->first][k->first]*d[k->first];

      L[i][j->first] = (1.0/d[j->first])*lld;
    }

    // i == k の場合
    double ld = A[i][i];

    L[i][0]; // L.sync()
    internal_map Li = L[i];
    for ( internal_map::iterator k = Li.begin(); k != Li.end(); k++) 
      ld -= L[i][k->first]*L[i][k->first]*d[k->first];
  
    d[i] = ld;
    L[i][i] = 1.0;
  }

  return 1;
}



/*!
 * 不完全コレスキー分解(incomplete Cholesky decomposition)
 *  - 対称行列A(n×n)を下三角行列(L:Lower triangular matrix)と対角行列の積(LDL^T)に分解する
 *  - l_ii * d_i = 1とした場合
 *  - L: i > jの要素が非ゼロで対角成分は1
 *  - 行列Aの値が0である要素に対応する部分を飛ばす
 * @param[in] A n×nの対称行列
 * @param[out] L 対角成分が1の下三角行列
 * @param[out] d 対角行列(対角成分のみ)
 * @param[in] n 行列の大きさ
 * @return 1:成功,0:失敗
 */
int IncompleteCholeskyDecomp2(matrix &A, matrix &L, vector<double> &d, int n)
{
  if(n <= 0) return 0;

  L[0][0] = A[0][0];
  d[0] = 1.0/L[0][0];

  for(int i = 1; i < n; ++i){


    for ( auto jt : A[i] ) { int j = jt.first;

      if(fabs(A[i][j]) < 1.0e-10) continue;

      double lld = A[i][j];

      for ( auto it : A[i] ) { int k = it.first;
	if ( k < j )
	  lld -= L[i][k]*L[j][k]*d[k];
      }
      L[i][j] = lld;
      L[j][i] = lld;
    }


    d[i] = 1.0/L[i][i];
  }

  return 1;
}
