#include <cstdlib>
#include "est/sparse.hpp"
#include <f2c.h>
extern "C" {
  int acutualmain(int argc, char **argv, integer *ib, integer l, integer *mb,
		  integer m);
};


int stwart(sparse::matrix<double>&A);
integer *generate_ma(sparse::matrix<double>&A);
integer generate_ia(sparse::matrix<double>&A, integer **ia);

int main(int argc, char **argv){
    static integer ib[100] =
      {
 	   2,                      10, 11,    13,    15, 16,   //  6
	         4, 5,             10,                   16,   // 10
	                        9,     11, 12,                 // 13
	1,          5,    7, 8,                                // 17
	1, 2,    4, 5, 6,       9, 10,                         // 24
	         4,             9,         12,                 // 27
	   2, 3, 4,                               14,          // 31
	1, 2,                   9,                        16,  // 35
	         4,                                   15,      // 37
	               6, 7, 8,                   14, 15,      // 42
	   2,                                             16,  // 44
	                        9,     11, 12,                 // 47
	               6,                         14, 15,      // 50
	                                              15,      // 51
	         4,    6,                  12,        15,      // 55
	   2,    4,                10,        13,              // 59
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0,0,0,0,
	0 };

    
    static integer mb[21] = { 0,6,10,13,17,24,27,31,35,37,42,44,47,50,51,55,
			      59,0,0,0,0 };

    sparse::matrix<double> A(17);
    A[ 1][1]=0; A[ 1][2]=2; A[ 1][3]=0; A[ 1][4]=0; A[ 1][5]=0; A[ 1][6]=0; A[ 1][7]=0; A[ 1][8]=0; A[ 1][9]=0; A[ 1][10]=10;A[ 1][11]=11;A[ 1][12]=0; A[ 1][13]=13;A[ 1][14]=0; A[ 1][15]=15;A[ 1][16]=16;
    A[ 2][1]=0; A[ 2][2]=0; A[ 2][3]=0; A[ 2][4]=4; A[ 2][5]=5; A[ 2][6]=0; A[ 2][7]=0; A[ 2][8]=0; A[ 2][9]=0; A[ 2][10]=10;A[ 2][11]=0; A[ 2][12]=0; A[ 2][13]=0; A[ 2][14]=0; A[ 2][15]=0; A[ 2][16]=16; 
    A[ 3][1]=0; A[ 3][2]=0; A[ 3][3]=0; A[ 3][4]=0; A[ 3][5]=0; A[ 3][6]=0; A[ 3][7]=0; A[ 3][8]=0; A[ 3][9]=9; A[ 3][10]=0; A[ 3][11]=11;A[ 3][12]=12;A[ 3][13]=0; A[ 3][14]=0; A[ 3][15]=0; A[ 3][16]=0;
    A[ 4][1]=1; A[ 4][2]=0; A[ 4][3]=0; A[ 4][4]=0; A[ 4][5]=5; A[ 4][6]=0; A[ 4][7]=7; A[ 4][8]=8; A[ 4][9]=0; A[ 4][10]=0; A[ 4][11]=0; A[ 4][12]=0; A[ 4][13]=0; A[ 4][14]=0; A[ 4][15]=0; A[ 4][16]=0; 
    A[ 5][1]=1; A[ 5][2]=2; A[ 5][3]=0; A[ 5][4]=4; A[ 5][5]=5; A[ 5][6]=6; A[ 5][7]=0; A[ 5][8]=0; A[ 5][9]=9; A[ 5][10]=10;A[ 5][11]=0; A[ 5][12]=0; A[ 5][13]=0; A[ 5][14]=0; A[ 5][15]=0; A[ 5][16]=0;
    A[ 6][1]=0; A[ 6][2]=0; A[ 6][3]=0; A[ 6][4]=4; A[ 6][5]=0; A[ 6][6]=0; A[ 6][7]=0; A[ 6][8]=0; A[ 6][9]=9; A[ 6][10]=0; A[ 6][11]=0; A[ 6][12]=12;A[ 6][13]=0; A[ 6][14]=0; A[ 6][15]=0; A[ 6][16]=0; 
    A[ 7][1]=0; A[ 7][2]=2; A[ 7][3]=3; A[ 7][4]=4; A[ 7][5]=0; A[ 7][6]=0; A[ 7][7]=0; A[ 7][8]=0; A[ 7][9]=0; A[ 7][10]=0; A[ 7][11]=0; A[ 7][12]=0; A[ 7][13]=0; A[ 7][14]=14;A[ 7][15]=0; A[ 7][16]=0;
    A[ 8][1]=1; A[ 8][2]=2; A[ 8][3]=0; A[ 8][4]=0; A[ 8][5]=0; A[ 8][6]=0; A[ 8][7]=0; A[ 8][8]=0; A[ 8][9]=9; A[ 8][10]=0; A[ 8][11]=0; A[ 8][12]=0; A[ 8][13]=0; A[ 8][14]=0; A[ 8][15]=0; A[ 8][16]=16; 
    A[ 9][1]=0; A[ 9][2]=0; A[ 9][3]=0; A[ 9][4]=4; A[ 9][5]=0; A[ 9][6]=0; A[ 9][7]=0; A[ 9][8]=0; A[ 9][9]=0; A[ 9][10]=0; A[ 9][11]=0; A[ 9][12]=0; A[ 9][13]=0; A[ 9][14]=0; A[ 9][15]=15;A[ 9][16]=0;
    A[10][1]=0; A[10][2]=0; A[10][3]=0; A[10][4]=0; A[10][5]=0; A[10][6]=6; A[10][7]=7; A[10][8]=8; A[10][9]=0; A[10][10]=0; A[10][11]=0; A[10][12]=0; A[10][13]=0; A[10][14]=14;A[10][15]=15;A[10][16]=0; 
    A[11][1]=0; A[11][2]=2; A[11][3]=0; A[11][4]=0; A[11][5]=0; A[11][6]=0; A[11][7]=0; A[11][8]=0; A[11][9]=0; A[11][10]=0; A[11][11]=0; A[11][12]=0; A[11][13]=0; A[11][14]=0; A[11][15]=0; A[11][16]=16;
    A[12][1]=0; A[12][2]=0; A[12][3]=0; A[12][4]=0; A[12][5]=0; A[12][6]=0; A[12][7]=0; A[12][8]=0; A[12][9]=9; A[12][10]=0; A[12][11]=11;A[12][12]=12;A[12][13]=0; A[12][14]=0; A[12][15]=0; A[12][16]=0; 
    A[13][1]=0; A[13][2]=0; A[13][3]=0; A[13][4]=0; A[13][5]=0; A[13][6]=6; A[13][7]=0; A[13][8]=0; A[13][9]=0; A[13][10]=0; A[13][11]=0; A[13][12]=0; A[13][13]=0; A[13][14]=14;A[13][15]=15;A[13][16]=0;
    A[14][1]=0; A[14][2]=0; A[14][3]=0; A[14][4]=0; A[14][5]=0; A[14][6]=0; A[14][7]=0; A[14][8]=0; A[14][9]=0; A[14][10]=0; A[14][11]=0; A[14][12]=0; A[14][13]=0; A[14][14]=0; A[14][15]=15;A[14][16]=0; 
    A[15][1]=0; A[15][2]=0; A[15][3]=0; A[15][4]=4; A[15][5]=0; A[15][6]=6; A[15][7]=0; A[15][8]=0; A[15][9]=0; A[15][10]=0; A[15][11]=0; A[15][12]=12;A[15][13]=0; A[15][14]=0; A[15][15]=15;A[15][16]=0;
    A[16][1]=0; A[16][2]=2; A[16][3]=0; A[16][4]=4; A[16][5]=0; A[16][6]=0; A[16][7]=0; A[16][8]=0; A[16][9]=0; A[16][10]=10;A[16][11]=0; A[16][12]=0; A[16][13]=13;A[16][14]=0; A[16][15]=0; A[16][16]=0; 


    stwart(A);

#if 1
    integer *ma, m;
    ma = generate_ma(A);
    m  = A.size()-1; 
    
    integer *ia, l;
    l = generate_ia(A,&ia);
    
    acutualmain(argc, argv, ia, l, ma, m);
#endif

    return 0;
}
