// https://detail.chiebukuro.yahoo.co.jp/qa/question_detail/q1423481709
/*
ソースの検討が完了しました。どうぞ御覧下さい。

解説：
完全ピボット選択はピボットを左上頂点とする正方行列
の中から絶対値最大の要素を選び、それがピボット位置
に来るように行と列を交換して消去を実行します。当然
完全ピボット選択の方が、ピボット選択の目的実現に適
しています。
しかし部分ピボット選択ではじめて精度が保たれるよう
な連立方程式の出現は、ピボット選択を必要としない連
立方程式にくらべて稀であり、完全ピボット選択でなけ
ればならないという場合は、更に頻度が低いと考えるの
が常識的です。
部分選択と完全選択によって救われる精度がいかほどか
というような定量的な情報を知らないので、ごく定性的
な話しか出来ません。私は通常、選択をする場合も部分
選択で満足しています。
以下に、完全ピボット選択付きのガウス消去法を示しま
す。テストは y=1+2x+3x^2+4x^3+5x^4 の係数を未知数
として、x=-2,-1,0,1,2 の５点での x,y から逆算する
問題を使いました。このソースはその時点のままです。
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// 静的配列サイズ仮設定
#define   N   100

// 係数配列参照マクロ（記述上の便宜のため）
#define   A(i,j)    a[N*i+j]
#define   B(i)      b[i]

// 連立方程式行列を表示
void printEq(double *a,double *b,int n){
    int   i,j;

    for(i=0;i<n;i++) {
        for(j=0;j<n;j++) {
            printf(" %9f",A(i,j));
          }
        printf("   %10f\n",B(i));
      }
}
// 完全ピボット選択
double perfectsel(
		    double  *a, // 連立方程式の係数行列：[N*N]
		    double  *b, // 連立方程式の右辺：[N]
		    int     *ib,// 変数整列インデクス[N]
		    int     n,  // 連立方程式の元数
		    int     k)  // ピボットインデクス
{
    double  abskk;  // ピボット比較絶対値
		      int     ik,jk;  // ピボット候補位置
				        int     ibw;    // ib[]交換退避ワード
						        double  aw,bw;  // a,b 交換退避ワード
									  int     i,j;
    double  Akk;

    // 変数整列インデクスを初期化
      if(k==0) {
        for(i=0;i<n;i++) ib[i] = i;
      }
    // No.kピボットの完全ピボット選択
      ik = jk = k;
    abskk = fabs(A(k,k));
    for(j=k;j<n;j++) {
        for(i=k;i<n;i++) {
            Akk = A(k,k);
            if(fabs(Akk)>abskk) {
	        abskk = fabs(Akk);
	        ik = i; jk = j;
	      }
          }
      }
    // 行の交換
      if(ik!=k) {
        for(j=k;j<n;j++) {
            aw=A(k,j);A(k,j)=A(ik,j);A(ik,j)=aw;
          } bw=B(k);B(k)=B(ik);B(ik)=bw;
      }
    // 列の交換
      if(jk!=k) {
        for(i=0;i<n;i++) {
            aw=A(i,k);A(i,k)=A(i,jk);A(i,jk)=aw;
          } ibw=ib[k];ib[k]=ib[jk];ib[jk]=ibw;
      }
    // 確定ピボットを戻す
      return A(k,k);
}
// ガウス消去：完全ピボット選択付
// 戻り値：0=normal, 1=0-pivot
int gausselim(
	        double *a,  // 連立方程式の係数行列：[N*N]
	        double *b,  // 連立方程式の右辺：[N]
	        double *x,  // 連立方程式の解ベクトル：[N]
	        int n)      // 連立方程式の元数
{
    int     ib[N];      // 未知数整列インデクス
		      double  eps=1.e-6;  // ０ピボット判定限界
					    int     i,j,k;
    double  Akk;
    double  Aik;

    // 前進消去
      for(k=0;k<n;k++) {
        // ピボット選択
          Akk = perfectsel(a,b,ib,n,k);
        if(fabs(Akk)<eps) {
            printf("Pivot = %.1e\n",Akk);
            return 1;
          }
        // ピボット行正規化
          for(j=k;j<n;j++) A(k,j) /= Akk;
        B(k) /= Akk;
        // ピボット下の消去
          for(i=k+1;i<n;i++) {
            Aik = A(i,k);
            for(j=k;j<n;j++) A(i,j) -= Aik*A(k,j);
            B(i) -= Aik*B(k);
          }
      }
    // 後退代入
      for(i=n-2;i>=0;i--) {
        for(j=i+1;j<n;j++) {
            B(i) -= A(i,j)*B(j);
          }
      }
    // 解ベクトル整列
      for(k=0;k<n;k++) {
        x[ib[k]] = b[k];
      }
    return 0;
}
// テスト用４次関数
double f(double x){
    return 1+x*(2+x*(3+x*(4+x*5)));
}
int main(){
    double  a[N*N]; // 連立方程式の係数行列：[N*N]
		       double  b[N];   // 連立方程式の右辺：[N]
				        double  x[N];   // 連立方程式の解ベクトル：[N]
						         int     n=5;    // 連立方程式の元数
								         int     err;
    int     i,j,k;
    double  xv;

    // 連立方程式設定
      // ４次代数式の係数(1,2,3,4,5)を求める
      printf("Equations to fit y=1+2x+3x^2+4x^3+5x^4\n");
    for(i=0;i<n;i++) {
        xv = (i-2);
        for(j=0;j<n;j++) {
            A(i,j) = pow(xv,j);
          } B(i) = f(xv);
      }
    printf("Matrix of equations:\n");
    printEq(a,b,n);
    // 解法
      err = gausselim(a,b,x,n);
    if(err) return 1;
    // 解ベクトル
      printf("Vector of solution:\n");
    for(k=0;k<n;k++) {
        printf("%10f\n",x[k]);
      }
    return 0;
}
#undef  A
#undef  N
/* 出力：
文字数制限のため省略します
*/

