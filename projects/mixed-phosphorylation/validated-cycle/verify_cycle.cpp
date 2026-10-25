// Copyright 2026 Zijian Zeng. Original checker code under the MIT license.
// Uses CAPD (GPL-3.0); see its retained upstream COPYING before redistribution.
// This is a research-stage certificate attempt, not a claim of priority.
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <stdexcept>
#include "capd/capdlib.h"
#include "seed.h"
using namespace capd;

interval max_row_sum(const IMatrix& A) {
  interval result(0.);
  for(int i=0;i<4;i++) {
    interval row(0.);
    for(int j=0;j<4;j++) row += abs(A[i][j]);
    if(row.rightBound()>result.rightBound()) result=interval(row.rightBound());
  }
  return result;
}

interval determinant4(const IMatrix& A) {
 int permutation[4]={0,1,2,3}; interval result(0.);
 do {
  int inversions=0; interval term(1.);
  for(int i=0;i<4;i++) { term*=A[i][permutation[i]];
   for(int j=i+1;j<4;j++) inversions += permutation[i]>permutation[j];
  }
  result += (inversions%2 ? -term : term);
 } while(std::next_permutation(permutation,permutation+4));
 return result;
}

int main(int argc,char** argv) {
 std::cout << std::setprecision(17);
 try {
  const double radius=argc>1 ? std::stod(argv[1]) : 1e-8;
  const int order=argc>2 ? std::stoi(argv[2]) : 30;
  if(!(radius>0.) || order<5) throw std::runtime_error("Invalid radius/order");
  IMap vf("par:r;var:u,v,w,a,b;fun:100*b-u,u-v,v-100*w*(8-a-b)+a,100*w*(8-a-b)-r*a,3*(40-u-v-w-a-b)*(8-a-b)-101*b;");
  vf.setParameter("r",interval(19)/interval(10));
  IOdeSolver solver(vf,order);
  solver.setAbsoluteTolerance(1e-14);
  solver.setRelativeTolerance(1e-14);
  ICoordinateSection section(5,0,interval(13)/interval(4));
  IPoincareMap pm(solver,section,poincare::MinusPlus);
  IVector c(4),X(4),q0(5),qX(5);
  q0[0]=qX[0]=interval(13)/interval(4);
  for(int i=0;i<4;i++) {
    c[i]=interval(center_values[i]);
    X[i]=c[i]+interval(-radius,radius);
    q0[i+1]=c[i]; qX[i+1]=X[i];
  }
  // Exact algebraic forward invariance of the physical domain is documented
  // separately. Here certify the whole initial box is strictly inside it.
  bool initial_positive=true;
  interval p=8-qX[3]-qX[4],s=40-qX[0]-qX[1]-qX[2]-qX[3]-qX[4];
  for(int i=0;i<5;i++) initial_positive &= qX[i].leftBound()>0.;
  initial_positive &= p.leftBound()>0. && s.leftBound()>0.;
  std::cout << "RADIUS " << radius << " ORDER " << order << "\nINITIAL_BOX " << qX << "\n";
  C0HOTripletonSet pointSet(q0);
  interval centerTime;
  const IVector ret0=pm(pointSet,centerTime);
  std::cout << "CENTER_RETURN_TIME " << centerTime << "\nCENTER_RETURN " << ret0 << "\n";
  IVector residual(4);
  for(int i=0;i<4;i++) residual[i]=ret0[i+1]-c[i];
  C1HORect2Set fullSet(qX);
  interval boxTime;
  IMatrix flowDerivative(5,5);
  const IVector retX=pm(fullSet,flowDerivative,boxTime);
  const IMatrix allDP=pm.computeDP(retX,flowDerivative);
  IMatrix D(4,4),R(4,4),eye(4,4);
  eye.setToIdentity();
  for(int i=0;i<4;i++)for(int j=0;j<4;j++) {
    D[i][j]=allDP[i+1][j+1]; R[i][j]=interval(inverse_values[i][j]);
  }
  IMatrix correction=eye-R*(D-eye);
  IVector krawczyk=c-R*residual+correction*(X-c);
  const bool inclusion=subsetInterior(krawczyk,X);
  const interval determinantR=determinant4(R);
  const bool preconditioner_invertible=determinantR.leftBound()>0. || determinantR.rightBound()<0.;
  const interval correctionNorm=max_row_sum(correction);
  const bool newton_contraction=correctionNorm.rightBound()<1.;
  // Any fixed-point Jacobian belongs to D. An enclosure of D^4 with
  // infinity norm <1 implies its spectral radius is strictly below 1.
  IMatrix power=D;
  for(int k=1;k<4;k++)power=power*D;
  const interval norm4=max_row_sum(power);
  const bool attracting=norm4.rightBound()<1.;
  const interval transverse=100*retX[4]-retX[0];
  const bool crossing=transverse.leftBound()>0.;
  std::cout << "BOX_RETURN_TIME " << boxTime << "\nBOX_RETURN " << retX
            << "\nRESIDUAL " << residual << "\nDP " << D
            << "\nDETERMINANT_R " << determinantR << "\nNEWTON_CONTRACTION_NORM " << correctionNorm
            << "\nKRAWCZYK " << krawczyk << "\nKRAWCZYK_INCLUSION " << inclusion
            << "\nDP_FOURTH_POWER_NORM_INF " << norm4
            << "\nTRANSVERSE " << transverse
            << "\nINITIAL_PHYSICAL_INTERIOR " << initial_positive
            << "\nATTRACTING_LINEARIZATION " << attracting
            << "\nPOSITIVE_CROSSING " << crossing << "\n";
  const bool ok=inclusion && preconditioner_invertible && newton_contraction && attracting && initial_positive && crossing && boxTime.leftBound()>0.;
  std::cout << "CERTIFICATE_CHECKS " << (ok ? "PASSED" : "FAILED") << "\n";
  return ok ? 0 : 2;
 } catch(const std::exception& e) {
  std::cerr << "VALIDATION_EXCEPTION " << e.what() << "\n"; return 3;
 }
}
