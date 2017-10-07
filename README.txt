CodeFEM - Codes for FEM

### 1st.
$ cd bin/; ./mkall; cd ..

### 2nd.
$ cd src/xmesh.d; make

### xmesh is a generator of FEM mesh.

### 3rd.
# mkdir /usr/include/est
# cp include/est/sparse.hpp   /usr/include/est
# cp include/est/bicgstab.hpp /usr/include/est
# cp src/dataparallel.d/cl_bicgstab_kernel.cl /usr/include/est
# cp lib/libbicgstab.a /usr/lib
# c++ sample.c++ -lbicgstab -lOpenCL
# Usage: ./a.out [OPTION]
# argv[1] = Number of the PE. s.t. ./a.out 16

### cl_bicgstab() is a linier solver using OpenCL.
/* sample.c++ --- A sample source for cl_bicgstab(). */
#include <cstdio>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>

int main(int argc, char **argv)
{
  cl_bicgstab_init(argc,argv);

  sparse::matrix<double> A(3);
  vector<double> b(3), x;

  A[1][1] = 2.0; A[1][2] = 1.0; b[1] = 3.0;
                 A[2][2] = 1.0; b[2] = 1.0;

  x = cl_bicgstab(A,b);
  printf("x[1]=%f\n",x[1]);
  printf("x[2]=%f\n",x[2]);
  return 0;
}


### 4th.
# cp include/est/navierstokes.hpp   /usr/include/est
# cp lib/libnavierstokes.a /usr/lib
# cp lib/libxmesh.a /usr/lib
# cp lib/libforeach.a /usr/lib
# c++ main.c++ -lbicgstab -lOpenCL -lnavierstokes -lxmesh -lforeach
# Usage: ./a.out [OPTION]
# argv[1] = Number of the PE. s.t. ./a.out 16

### This is a simulation of Navier-Stokes equations.
/* main.c++ --- A sample source for Navier-Stokes equations. */
#include <cstdio>
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/navierstokes.hpp>

int main(int argc, char **argv)
{
  double Re, dt;

  cl_bicgstab_init(argc,argv);
  if(0!=navierstokes_init("cavity32.mesh",Re=5000,dt=0.001))
    return 0;

  sparse::matrix<double> A; vector<double> U, b;

  for ( int T=0; T<= 36000000; T++) {
    fprintf(stderr,"T");
    navierstokes(A,U,b);

    fprintf(stderr,"=");
    U = cl_bicgstab(A,b);

    fprintf(stderr,"%05d\n",T);
    plotuv(U);
    if ( 0 == (T%1000)) fprintuv(U);
  }
  return 0;
}


### 5th.
# cp include/est/psc98.hpp   /usr/include/est
# cp lib/libpsc98.a /usr/lib
# c++ cl_main.c++ -lpsc98 -lbicgstab -lOpenCL -o cl_a.out
# Usage: ./cl_a.out [OPTION]
# argv[1] = Number of the PE. s.t. ./a.out 16

### Those are Parallel Software Contest '98 (PSC98)'s problems.
/* cl_main.c++ --- The problem solver of PSC98 by OpenCL. */
#include <vector>
#include <est/sparse.hpp>
#include <est/bicgstab.hpp>
#include <est/psc98.hpp>

int main(int argc, char **argv){
  cl_bicgstab_init(argc,argv);

  sparse::matrix<double> A; vector<double> x, b;
  psc98_init(A,b);
  x = cl_bicgstab(A,b);
  psc98_check(x);
  return 0;
}

# cd src/psc98.d
# make time
export PSC98=0; time ./cl_a.out 120
Problem NO : 0
|b - Ax|_inf = 1.16993e-07  (OK)

real 0m1.035s
user 0m0.702s
sys  0m0.334s
export PSC98=1; time ./cl_a.out 120
Problem NO : 1
|b - Ax|_inf = 2.58968e-07  (OK)

real 0m53.176s
user 0m39.063s
sys  0m14.147s
export PSC98=2; time ./cl_a.out 120
Problem NO : 2
|b - Ax|_inf = 9.93914e-08  (OK)

real 0m56.498s
user 0m41.382s
sys  0m15.151s
export PSC98=3; time ./cl_a.out 120
Problem NO : 3
|b - Ax|_inf = 3.2271e-07  (OK)

real 1m7.244s
user 0m50.220s
sys  0m17.066s
export PSC98=4; time ./cl_a.out 120
|b - Ax|_inf = 1.76181e-07  (OK)

real 9m18.472s
user 6m53.404s
sys  2m25.411s
export PSC98=5; time ./cl_a.out 120
Problem NO : 5
|b - Ax|_inf = 8.60549e+146  (NG)

real 15m40.922s
user 11m29.524s
sys  4m11.966s

# The problem NO.5 is ill condition matrix. So cl_bicgstab() can't solve it.


参考文献

[ 1] 株式会社フィックスターズ: OpenCL入門, インプレスジャパン, 2010.

[ 2] 小国 力 編著: 行列計算ソフトウェア, 丸善株式会社, 1992.

[ 3] 中山 司: 流れ解析のための有限要素法入門, 東京大学出版会, 2008.

[ 4] 粂井康孝: 猫でもわかるC++プログラミング, ソフトバンク クリエイティブ株式会社, 2009.

[ 5] 佃 良生, 海津 聰: デローニー法によるFEM格子の自動形成, 数値流体力学会講演論文集, 1996.

[ 6] Richard Barrant/Michael Berry/Tony F. Chan/James Demmel/June Donato/Jack Dongarra/Victor Eijhout/Roldan Pozo/Charles Romine/Henk van der Vorst 著, 長谷川 里美/長谷川 秀彦/藤野 清次 訳: 反復法Templates, 朝倉書店, 1996.

[ 7] 儀我 美一, 儀我 美保: 非線形偏微分方程式, 共立出版株式会社, 1999.

[ 8] 佃 良生: 複雑形状の領域における流れの数値計算について, 電気通信大学平成9年度情報工学科卒業論文, 1998.

[ 9] 矢川 元基, 奥田 洋司, 中林 靖: 有限要素法流れ解析, 朝倉書店, 1998.

[10] 日本計算工学会 流れの有限要素法研究委員会 編: 有限要素法による流れのシミュレーション, シュプリンガー・フェアラーク東京, 1998.

[11] 日本計算工学会 続・流れの有限要素法研究委員会 編: 有限要素法による流れのシミュレーション, シュプリンガー・フェアラーク東京, 2008.

[12] 川原 睦人: 有限要素法流体解析, 日科技連, 1985.

--
ダウンロード

[ 8] については, http://kud.dip.jp/thesis.d/UEC_thesis.pdf

