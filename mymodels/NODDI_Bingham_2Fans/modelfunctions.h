#include "diffusivities.h"
#include "hyp_Sapprox.h"
#include "hyp_Sapprox_doublePrecision.h"
#include "eigenValues.h"

// Parameters in NODDI-Bingham
// P[0]: fiso
// P[1]: fintra
// P[2]: ffibre2
// P[3]: kappa1
// P[4]: beta1
// P[5]: th1
// P[6]: ph1
// P[7]: psi1
// P[8]: kappa2
// P[9]: beta2
// P[10]: th2
// P[11]: ph2
// P[12]: psi2

// CFP[0:2] are bvecs 
// CFP[3] are bvals

// FixP[0] is S0

MACRO T Calculate_ExtraCellular(int nFan,
				T* P,
				T* CFP,
				T &xv,
				T &gradXLocal,
				T &gradYLocal)
{
  // Calculate Signal from ExtraCellular compartment //
  T kappa, beta;
  if(nFan==1){
    kappa=P[3];
    beta=P[4];
  }else{
    kappa=P[8];
    beta=P[9];
  }

  //dPerp = dPar*(1-f);
  T Dperpendicular = Dparallel*(1-P[1]);

  //BinghamHinderedDiffusionCoeff
  T dPar; // diffusivity along the dominant orientation
  T dPerp1; // diffusivity along the fanning direction
  T dPerp2; // diffusivity along the final direction
  //[dPar dPerp1 dPerp2] = BinghamHinderedDiffusionCoeff(x(1), x(2), x(3), x(4));
  T dParMdPerp = Dparallel - Dperpendicular;
  T delta = (T)0.00005;

  T nc =  hyp_Sapprox_doublePrecision(kappa,beta,(T)0.0);
  //computed using the derivative relation with finite difference
  T diff =  hyp_Sapprox_doublePrecision(kappa+delta,beta,(T)0.0) - hyp_Sapprox_doublePrecision(kappa-delta,beta,(T)0.0);
  diff = diff/((T)2.0*delta);
  dPar = Dperpendicular + dParMdPerp*diff/nc;

  diff = hyp_Sapprox_doublePrecision(kappa,beta+delta,(T)0.0) - hyp_Sapprox_doublePrecision(kappa,beta-delta,(T)0.0);
  diff = diff/((T)2.0*delta);
  dPerp1 = Dperpendicular + dParMdPerp*diff/nc;
 
  //computed using the trace constancy
  dPerp2 = Dparallel + 2*Dperpendicular - dPar - dPerp1;
  ////////

  T cosThetaSq = xv*xv;
  T sinThetaSq = 1-(cosThetaSq);
  T phi= atan2_gpu(gradYLocal, gradXLocal);
  T cosPhi = cos_gpu(phi);
  T cosPhiSq = cosPhi*cosPhi;
  T sinPhiSq = 1-cosPhiSq;
  
  T ExtraCell = exp_gpu(-CFP[3]*(dPar*cosThetaSq + sinThetaSq*(dPerp1*sinPhiSq + dPerp2*cosPhiSq)));
      
  return ExtraCell;
}

MACRO T Calculate_IntraCellular(int nFan,
				T* P,
				T* CFP,
				T& xv,
				T &gradXLocal,
				T &gradYLocal)
{
  // Calculate Signal from IntraCellular compartment //
  T kappa, beta;
  if(nFan==1){
    kappa=P[3];
    beta=P[4];
  }else{
    kappa=P[8];
    beta=P[9];
  }
  
  T parComp =-CFP[3]*Dparallel; // Parallel component: -bval * dintra
  // Radius is 0, so no Perpendicular Component

  T BinghamNC = 1/hyp_Sapprox_doublePrecision(kappa,beta,(T)0.0);

  T Matrix[9];
  // matrix = [0 0 0; 0 beta 0; 0 0 kappa];
  // projG = [grad_dirs(i,:)*normaldir grad_dirs(i,:)*fanningdir grad_dirs(i,:)*fibredir]';
  // matrix = matrix + Lpmp(i)*projG*projG'; 
  Matrix[0] = gradXLocal*gradXLocal*parComp;
  Matrix[1] = gradXLocal*gradYLocal*parComp;
  Matrix[2] = gradXLocal*xv*parComp;
  Matrix[3] = Matrix[1];
  Matrix[4] = beta + gradYLocal*gradYLocal*parComp;
  Matrix[5] = gradYLocal*xv*parComp;
  Matrix[6] = Matrix[2];
  Matrix[7] = Matrix[5];
  Matrix[8] = kappa + xv*xv*parComp;

  T eigenvalues[3];
  getEigenValues_symm_3x3(Matrix,eigenvalues);
  
  T IntraCell = (hyp_Sapprox_doublePrecision(eigenvalues[0],eigenvalues[1],eigenvalues[2])*BinghamNC);
  
  return IntraCell;
}

MACRO void get_GradLocals_1(
			    T* P,
			    T* CFP,
			    T* fibredir,
			    T* gradsLocal)
{  
  // gradYLocal = grad_dirs*fanningdir;
  // fanningdir = [ -cosPsi.*sinPhi - cosTheta.*cosPhi.*sinPsi cosPhi.*cosPsi - cosTheta.*sinPhi.*sinPsi sinTheta.*sinPsi]';
  T fanningdir[3];
  fanningdir[0] = -cos_gpu(P[7])*sin_gpu(P[6]) - cos_gpu(P[5])*cos(P[6])*sin_gpu(P[7]);
  fanningdir[1] = cos_gpu(P[6])*cos_gpu(P[7]) - cos_gpu(P[5])*sin_gpu(P[6])*sin_gpu(P[7]);
  fanningdir[2] = sin_gpu(P[5])*sin_gpu(P[7]);

  gradsLocal[1] =  CFP[0]*fanningdir[0] +  CFP[1]*fanningdir[1] + CFP[2]*fanningdir[2];
  
  T crossProduct_fanfib[3];
  crossProduct_fanfib[0] = fanningdir[1]*fibredir[2] - fanningdir[2]*fibredir[1];
  crossProduct_fanfib[1] = fanningdir[2]*fibredir[0] - fanningdir[0]*fibredir[2];
  crossProduct_fanfib[2] = fanningdir[0]*fibredir[1] - fanningdir[1]*fibredir[0];
  
  // gradXLocal = grad_dirs*cross(fanningdir,fibredir);
  gradsLocal[0] =  CFP[0]*crossProduct_fanfib[0] +  CFP[1]*crossProduct_fanfib[1] + CFP[2]*crossProduct_fanfib[2];
}

MACRO void get_GradLocals_2(
			    T* P,
			    T* CFP,
			    T* fibredir,
			    T* gradsLocal)
{  
  // gradYLocal = grad_dirs*fanningdir;
  // fanningdir = [ -cosPsi.*sinPhi - cosTheta.*cosPhi.*sinPsi cosPhi.*cosPsi - cosTheta.*sinPhi.*sinPsi sinTheta.*sinPsi]';
  T fanningdir[3];
  fanningdir[0] = -cos_gpu(P[12])*sin_gpu(P[11]) - cos_gpu(P[10])*cos(P[11])*sin_gpu(P[12]);
  fanningdir[1] = cos_gpu(P[11])*cos_gpu(P[12]) - cos_gpu(P[10])*sin_gpu(P[11])*sin_gpu(P[12]);
  fanningdir[2] = sin_gpu(P[10])*sin_gpu(P[12]);

  gradsLocal[1] =  CFP[0]*fanningdir[0] +  CFP[1]*fanningdir[1] + CFP[2]*fanningdir[2];
  
  T crossProduct_fanfib[3];
  crossProduct_fanfib[0] = fanningdir[1]*fibredir[2] - fanningdir[2]*fibredir[1];
  crossProduct_fanfib[1] = fanningdir[2]*fibredir[0] - fanningdir[0]*fibredir[2];
  crossProduct_fanfib[2] = fanningdir[0]*fibredir[1] - fanningdir[1]*fibredir[0];
  
  // gradXLocal = grad_dirs*cross(fanningdir,fibredir);
  gradsLocal[0] =  CFP[0]*crossProduct_fanfib[0] +  CFP[1]*crossProduct_fanfib[1] + CFP[2]*crossProduct_fanfib[2];
}

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, 	// Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{
  T isoterm = exp_gpu(-CFP[3]*Diso); 	// exp(-bval*d)

  // cosTheta
  T fibredir[3];
  fibredir[0] = sin_gpu(P[5])*cos_gpu(P[6]);
  fibredir[1] = sin_gpu(P[5])*sin_gpu(P[6]);
  fibredir[2] = cos_gpu(P[5]);
  T xv_1 = CFP[0]*fibredir[0] + CFP[1]*fibredir[1] + CFP[2]*fibredir[2];
  T gradsLocal_1[2];
  get_GradLocals_1(P,CFP,fibredir,gradsLocal_1);

  T ExtraCell = Calculate_ExtraCellular(1,P,CFP,xv_1,gradsLocal_1[0],gradsLocal_1[1]);
  T IntraCell = Calculate_IntraCellular(1,P,CFP,xv_1,gradsLocal_1[0],gradsLocal_1[1]);
  T Fib_1 = (1-P[1])*ExtraCell + P[1]*IntraCell;

  fibredir[0] = sin_gpu(P[10])*cos_gpu(P[11]);
  fibredir[1] = sin_gpu(P[10])*sin_gpu(P[11]);
  fibredir[2] = cos_gpu(P[10]);
  T xv_2 = CFP[0]*fibredir[0] + CFP[1]*fibredir[1] + CFP[2]*fibredir[2];
  T gradsLocal_2[2];
  get_GradLocals_2(P,CFP,fibredir,gradsLocal_2);

  ExtraCell = Calculate_ExtraCellular(2,P,CFP,xv_2,gradsLocal_2[0],gradsLocal_2[1]);
  IntraCell = Calculate_IntraCellular(2,P,CFP,xv_2,gradsLocal_2[0],gradsLocal_2[1]);
  T Fib_2 =  (1-P[1])*ExtraCell + P[1]*IntraCell;

  T Anisoterm = (1-P[2])*Fib_1 + P[2]*Fib_2;

  T pred_signal=FixP[0]*(P[0]*isoterm + (1-P[0])*Anisoterm);
    
  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) // Estimated parameters
{
  if(P[3]<P[4]){
    // kappa1 < beta1
    return false;
  }

  if(P[8]<P[9]){
    // kappa2 < beta2
    return false;
  }

  return true;
}

// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimate
			       T* P, // Estimated parameters, use P*
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
  derivatives[0]=NUMERICAL(0);
  derivatives[1]=NUMERICAL(1);
  derivatives[2]=NUMERICAL(2);
  derivatives[3]=NUMERICAL(3);
  derivatives[4]=NUMERICAL(4);
  derivatives[5]=NUMERICAL(5);
  derivatives[6]=NUMERICAL(6);
  derivatives[7]=NUMERICAL(7);
  derivatives[8]=NUMERICAL(8);
  derivatives[9]=NUMERICAL(9);
  derivatives[10]=NUMERICAL(10);
  derivatives[11]=NUMERICAL(11);
  derivatives[12]=NUMERICAL(12);

  
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{
  if(P[3]<P[4]){
    // kappa1 < beta1
    P[4]=P[3];
  }

  if(P[8]<P[9]){
    // kappa2 < beta2
    P[9]=P[8];
  }
}


