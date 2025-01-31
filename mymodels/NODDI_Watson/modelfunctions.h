
#include "diffusivities.h"
#include "WatsonFunctions.h"
#include "LegendreFunctions.h"

// Parameters in NODDI-Watson 
// P[0]: fiso
// P[1]: fintra
// P[2]: kappa
// P[3]: th
// P[4]: ph
// P[5]:d_par

// CFP[0:2] are bvecs 
// CFP[3] are bvals
// FixP[0] is S0

MACRO T step_size(T param, T scale){
  T step;
  step = (T)TINY*nonzerosign(param)*fabs_gpu(param)*scale;
  step = min_gpu(max_gpu(fabs_gpu(step),(T)MINSTEP),(T)MAXSTEP);
  
  return step;
}

MACRO T Calculate_ExtraCellular(T* P,
				T* CFP,
				T& xv)
{
  // Calculate Signal from ExtraCellular compartment //
  
  //dPerp = dPar*(1-f);
  T Dperpendicular = P[5]*(1-P[1]);
  
  T Dpar_equivalent;
  T Dperp_equivalent;
  
  WatsonHinderedDiffusionCoeff(P[2],P[5],Dperpendicular,Dpar_equivalent,Dperp_equivalent);
  
  //exp(-bval.*((dPar - dPerp)*cosThetaSq + dPerp));
  T ExtraCell = exp_gpu(-CFP[3]*((Dpar_equivalent - Dperp_equivalent)*xv*xv + Dperp_equivalent));
  
  return ExtraCell;
}

MACRO T Calculate_IntraCellular(T* P,
				T* CFP,
				T& xv)
{
  // Calculate Signal from IntraCellular compartment //

  T parComp =-CFP[3]*P[5]; // Parallel component: -bval * dintra
  // Radius is 0, so no Perpendicular Component

  // Compute Legendre weighted signal
  T lgi[7]; // legendre  gaussian integral
  legendreGaussianIntegral(-parComp,lgi);
  
  // Compute the spherical harmonic coefficients of the Watson's distribution
  T coeff[7];
  WatsonSHCoeff(P[2],coeff);
  
  for(int i=0;i<7;i++){
    lgi[i]*=coeff[i];
  }
  
  if(xv>1) xv=1;
  if(xv<-1) xv=-1;
  // Compute the SH values at cosTheta
  T SH[7];
  computeSH_values(xv,SH);
  
  for(int i=0;i<7;i++){
    lgi[i]*=SH[i];
  }
  
  T IntraCell = (T)0.0;
  for(int i=0;i<7;i++){
    IntraCell+=lgi[i];
  }
  
  // Remove negative values
  if(IntraCell<0) IntraCell=1e-3;
  IntraCell *= (T)0.5;

  return IntraCell;
}

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimate
			 T* P, 	// Estimated parameters
			 T* CFP, // Fixed Parameters common to all the voxels
			 T* FixP) // Fixed Parameters for each voxel
{
  T isoterm = exp_gpu(-CFP[3]*Diso); 	// exp(-bval*d)

  T xv = CFP[0]*sin_gpu(P[3])*cos_gpu(P[4])	// (bvec(1)*sinth*cosph
	+ CFP[1]*sin_gpu(P[3])*sin_gpu(P[4]) 	// + bvec(2)*sinth*sinph
	+ CFP[2]*cos_gpu(P[3]);			// + bvec(3)*costh)

  T ExtraCell = Calculate_ExtraCellular(P,CFP,xv);
  T IntraCell = Calculate_IntraCellular(P,CFP,xv);
		
  T anisoterm = (1-P[1])*ExtraCell+P[1]*IntraCell;
  
  T pred_signal=FixP[0]*((1-P[0])*anisoterm + P[0]*isoterm);
    
  return pred_signal;
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
		       int npar, // Number of Parameters to estimate
		       T* P) // Estimated parameters
{
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
  T isoterm = exp_gpu(-CFP[3]*Diso); 	// exp(-bval*d)

  T xv = CFP[0]*sin_gpu(P[3])*cos_gpu(P[4])	// (bvec(1)*sinth*cosph
	+ CFP[1]*sin_gpu(P[3])*sin_gpu(P[4]) 	// + bvec(2)*sinth*sinph
	+ CFP[2]*cos_gpu(P[3]);			// + bvec(3)*costh)

  T ExtraCell = Calculate_ExtraCellular(P,CFP,xv);
  T IntraCell = Calculate_IntraCellular(P,CFP,xv);
		
  T anisoterm = (1-P[1])*ExtraCell+P[1]*IntraCell;
  
  //T pred_signal=FixP[0]*((1-P[0])*anisoterm + P[0]*isoterm);
 
  // df/fiso
  derivatives[0]= FixP[0]*(isoterm-anisoterm);
 
  // df/fintra
  derivatives[1]= FixP[0]*(1-P[0])*(IntraCell - ExtraCell);

  // df/dk   
  // numerical differentiation
  derivatives[2]=NUMERICAL(2);

  // df/dth   
  // numerical differentiation
  derivatives[3]=NUMERICAL(3);

  // df/dph   
  // numerical differentiation
  derivatives[4]=NUMERICAL(4);
}

// Constraints run after LevenbergMarquardt (if LevenbergMarquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{
}

MACRO T custom_priors(
       int id_p,   // the number of parameter in the model (starts at 0)
			 T* P, 		// Estimated parameters
       int nmeas, // Number of measurements per voxel
			 T* CFP, 	// Fixed Parameters common to all the voxels for all measurements !!
			 T* FixP) 	// Fixed Parameters for each voxel
{
	return (T)0.0;
}

