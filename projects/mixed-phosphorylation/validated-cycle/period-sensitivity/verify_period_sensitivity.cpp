// Copyright 2026 Zijian Zeng. Original code under the MIT license.
// CAPD is GPL-3.0. Mac runs are development evidence pending Linux reproduction.
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <stdexcept>
#include "capd/capdlib.h"
#include "seed.h"
using namespace capd;

interval matrix_norm(const IMatrix& A) {
 interval out(0.);
 for(int i=0;i<4;i++) {
  interval row(0.);
  for(int j=0;j<4;j++) row+=abs(A[i][j]);
  if(row.rightBound()>out.rightBound())out=interval(row.rightBound());
 }
 return out;
}
interval vector_norm(const IVector& x) {
 interval out(0.);
 for(int i=0;i<4;i++) if(abs(x[i]).rightBound()>out.rightBound())out=interval(abs(x[i]).rightBound());
 return out;
}
interval determinant4(const IMatrix& A) {
 int p[4]={0,1,2,3}; interval out(0.);
 do {
  int inversions=0; interval term(1.);
  for(int i=0;i<4;i++) {
   term*=A[i][p[i]];
   for(int j=i+1;j<4;j++)inversions+=p[i]>p[j];
  }
  out+=(inversions%2 ? -term : term);
 }while(std::next_permutation(p,p+4));
 return out;
}
int main(int argc,char** argv) {
 std::cout<<std::setprecision(17);
 try {
  const double radius=argc>1 ? std::stod(argv[1]) : 1e-9;
  const int order=argc>2 ? std::stoi(argv[2]) : 30;
  if(!(radius>0.) || order<5)throw std::runtime_error("Invalid radius/order");
  // e is a constant state; its variational column is the parameter derivative.
  // The fixed certified orbit is at e=0. This polynomial family is the first
  // order jet of the actual slow-manifold vector field, not its full formula.
  IMap vf("par:r;var:u,v,w,a,b,e;fun:100*b-u+e*(u+10000*b-300*(8-a-b)*(40-u-v-w-a-b)),u-v,v-100*w*(8-a-b)+a,100*w*(8-a-b)-r*a,3*(40-u-v-w-a-b)*(8-a-b)-101*b-3*e*(8-a-b)*(u+100*b),0;");
  vf.setParameter("r",interval(19)/interval(10));
  IOdeSolver solver(vf,order);
  solver.setAbsoluteTolerance(1e-14);solver.setRelativeTolerance(1e-14);
  ICoordinateSection section(6,0,interval(13)/interval(4));
  IPoincareMap pm(solver,section,poincare::MinusPlus);
  IVector qX(6),c(4),X(4);
  qX[0]=interval(13)/interval(4);qX[5]=interval(0.);
  for(int i=0;i<4;i++) {
   c[i]=interval(center_values[i]);
   X[i]=c[i]+interval(-radius,radius);qX[i+1]=X[i];
  }
  std::cout<<"INITIAL_BOX "<<qX<<"\n";
  IVector q0(qX);
  for(int i=0;i<4;i++)q0[i+1]=c[i];
  C0HOTripletonSet pointSet(q0);interval point_time;
  const IVector ret0=pm(pointSet,point_time);
  IVector fixed_residual(4);
  for(int i=0;i<4;i++)fixed_residual[i]=ret0[i+1]-c[i];
  C1HORect2Set set(qX);
  interval return_time;IMatrix flowDerivative(6,6);
  const IVector ret=pm(set,flowDerivative,return_time);
  const IMatrix DP=pm.computeDP(ret,flowDerivative);
  IMatrix A(4,4),R(4,4),eye(4,4);eye.setToIdentity();
  IVector rhs(4);
  for(int i=0;i<4;i++) {
   rhs[i]=-DP[i+1][5];
   for(int j=0;j<4;j++) {
    A[i][j]=DP[i+1][j+1]-eye[i][j];R[i][j]=interval(inverse_values[i][j]);
   }
  }
  const IMatrix E=eye-R*A;
  const IVector krawczyk=c-R*fixed_residual+E*(X-c);
  const bool fixed_point_inclusion=subsetInterior(krawczyk,X);
  const interval eta=matrix_norm(E),detR=determinant4(R);
  const bool invertibleR=detR.leftBound()>0. || detR.rightBound()<0.;
  if(!fixed_point_inclusion || !invertibleR || !(eta.rightBound()<1.)) {
   std::cout<<"KRAWCZYK "<<krawczyk<<"\nFIXED_POINT_INCLUSION "<<fixed_point_inclusion<<"\nKRAWCZYK "<<krawczyk<<"\nFIXED_POINT_INCLUSION "<<fixed_point_inclusion
           <<"\nDETERMINANT_R "<<detR<<"\nPRECONDITIONED_REMAINDER_NORM "<<eta<<"\n";
   throw std::runtime_error("Preconditioned sensitivity solve not validated");
  }
  const IVector r=R*rhs;
  // For the actual A and rhs at the certified fixed point, z_e = r+E*z_e.
  // Neumann convergence proves existence/invertibility, and this norm bound
  // encloses the unique solution before interval fixed-point refinement.
  const interval bound=vector_norm(r)/(interval(1.)-eta);
  IVector ze(4);
  for(int i=0;i<4;i++)ze[i]=interval(-bound.rightBound(),bound.rightBound());
  for(int k=0;k<12;k++) {
   const IVector next=r+E*ze;
   for(int i=0;i<4;i++) {
    const double low=std::max(ze[i].leftBound(),next[i].leftBound());
    const double high=std::min(ze[i].rightBound(),next[i].rightBound());
    if(low>high)throw std::runtime_error("Empty sensitivity enclosure");
    ze[i]=interval(low,high);
   }
  }
  const interval normal=100*ret[4]-ret[0]; // e=0 exactly
  if(!(normal.leftBound()>0.))throw std::runtime_error("Return is not certified transverse");
  IVector tau_z(4);
  for(int j=0;j<4;j++)tau_z[j]=-flowDerivative[0][j+1]/normal;
  const interval tau_e=-flowDerivative[0][5]/normal;
  interval T1=tau_e;
  for(int j=0;j<4;j++)T1+=tau_z[j]*ze[j];
  const bool nonzero=T1.leftBound()>0. || T1.rightBound()<0.;
  std::cout<<"RETURN_TIME "<<return_time<<"\nRETURN_BOX "<<ret
           <<"\nFLOW_DERIVATIVE "<<flowDerivative<<"\nRETURN_MAP_DERIVATIVE "<<DP
           <<"\nKRAWCZYK "<<krawczyk<<"\nFIXED_POINT_INCLUSION "<<fixed_point_inclusion
           <<"\nDETERMINANT_R "<<detR<<"\nPRECONDITIONED_REMAINDER_NORM "<<eta
           <<"\nRAW_SENSITIVITY_NORM_BOUND "<<bound<<"\nINITIAL_POINT_SENSITIVITY "<<ze
           <<"\nNORMAL_VELOCITY "<<normal<<"\nTAU_Z "<<tau_z
           <<"\nTAU_E "<<tau_e<<"\nT1 "<<T1
           <<"\nT1_EXCLUDES_ZERO "<<nonzero<<"\n";
  std::cout<<"SENSITIVITY_CHECKS "<<(nonzero ? "PASSED" : "FAILED")<<"\n";
  return nonzero ? 0 : 2;
 } catch(const std::exception& ex) {
  std::cerr<<"VALIDATION_EXCEPTION "<<ex.what()<<"\n";return 3;
 }
}
