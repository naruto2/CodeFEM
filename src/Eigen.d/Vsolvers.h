gpuvector btogpu(Vector &b) {
  int n = b.size();
  gpuvector bgpu(n);
  for(int i=0;i<n;i++) bgpu[i]=b[i];
  return bgpu;
}


Vector fromgpu(gpuvector xgpu) {
  Vector x(xgpu.size());
  copy(xgpu.begin(), xgpu.end(), &x[0]);
  return x;
}


gpumatrix Atogpu(Smatrix &A) {
  int n = A.cols();
  gpumatrix Agpu(n,0);
  cpumatrix Acpu(n);

  Acpu.clear();
  Acpu.resize(n);

  for(int i=0;i<A.outerSize();++i){
    for(Smatrix::InnerIterator it(A,i);it;++it){
      Acpu[it.row()][it.col()] = it.value();
    }
  }
  copy(Acpu, Agpu);
  return Agpu;
}


// BiCGSTAB法による求解
Vector Vbicgstab(Smatrix &A, Vector &b) {
  bicgstab_tag tag(1e-8,1000000);  
  return fromgpu(solve(Atogpu(A), btogpu(b), tag));
}


// CG法による求解
Vector Vcg(Smatrix &A, Vector &b) {
  cg_tag tag(1e-8,1000000);
  if ( IsSymmetric(A) )
    return fromgpu(solve(Atogpu(A), btogpu(b), tag));
  return b;
}


// GMRES法による求解
Vector Vgmres(Smatrix &A, Vector &b) {
  gmres_tag tag(1e-8, 100000, 40);
  return fromgpu(solve(Atogpu(A), btogpu(b), tag));
}


// Vogita 荻田氏の方法(対称行列化)
Vector Vogita(Smatrix &A, Vector &b) {

  Smatrix AT = A.transpose();

  Smatrix AA = vcatSmatrix( catSmatrix(A+AT, A-AT),
			    catSmatrix(AT-A,-A-AT));

  Vector  bb = catVector(2*b, -2*b);

  return halfVector(Vcg(AA,bb));
}
