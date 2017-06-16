// 4 functions to be implemented:
// - The model predicted signal
// - Constraints during MCMC (Optional)
// - Partial derivatives for Levenberg-Marquardt (Optional)
// - Constraints after Levenberg-Marquardt (Optional)

// P[0]: a
// P[1]: b
// CFP[0]: Xs

MACRO T Predicted_Signal(
			 int npar, 	// Number of Parameters to estimate
			 T* P, 		// Estimated parameters
			 T* CFP, 	// Fixed Parameters common to all the voxels
			 T* FixP) 	// Fixed Parameters for each voxel
{
  return exp_gpu(-P[0]*CFP[0])*P[1];
}

// Constraints checked during MCMC (if MCMC is used)
MACRO bool ConstraintsMCMC(
			   int npar, // Number of Parameters to estimate
			   T* P) // Estimated parameters
{
  // You can add any constraints (but remember that you can also specify bounds in a different file)
  // These are useful for specify relation between parameters. For instance:
  // if (P[3]>(P[4]+3.14)) return false;

  // Do not modify this.
  return true;
}

// Partial derivatives respect each model parameter
// If Levenberg–Marquardt algorithm is used, the value of the partial derivative for each parameter has to be stored in the outpu array derivatives
// You can use Numerical differentiation using the keyword "NUMERICAL"
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimate
			       T* P, // Estimated parameters, use P*
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* FixP, // Fixed Parameters for each voxel
			       T* derivatives) // Derivative respect each model estimated parameter
{
	derivatives[0]=NUMERICAL(0);
	derivatives[1]=NUMERICAL(1);

}

// Constraints run after LevenbergMarquardt (if Levenberg-Marquardt is used)
MACRO void FixConstraintsLM(	
			    int npar, // Number of Parameters to estimate
			    T* P) // Estimated parameters
{
  // You can specify what to do with some parameters after Levenberg-Marquardt if a constraint is not satisfied.
  // For instance:
  // if(P[2]>1.0]) P[2]=1.0; 
}


