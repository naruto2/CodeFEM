// ogita 荻田氏の方法(対称行列化)
Vector ogita(Smatrix &A, Vector &b) {

  Smatrix AT = A.transpose();

  Smatrix AA = vcatSmatrix( catSmatrix(A+AT, A-AT),
			    catSmatrix(AT-A,-A-AT));

  Vector  bb = catVector(2*b, -2*b);
  
  return halfVector(cg(AA,bb));
}

