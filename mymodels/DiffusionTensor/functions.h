/* This file must contain:

 - Predicted_Signal function. It should be written as a C function using the marker __device__ at the beginning. It can call to subroutines included in this file, Subroutines must been defined before called and use the marker __device__ as well. 
   Do not use the Cost Function !!

 - Partial_Derivatives() function of the predicted_signal function (NO the Cost Function !!) respect each estimated parameter in the model. Use the marker __device__
   It must write an array of values (results), one for each derivate. Write as many results as model estimated parameters.

 Do not change the name of the functions.
 Do not change the function parameters.
 Do nit change the MACRO and do not remove it.

 Use T datatype (can be defined at copililation_parameters/types.h, float or double)
 
 The toolboox will use automatically the sum of squares as cost function SUM( (y-f(x))^2 )

 http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__SINGLE.html#group__CUDA__MATH__SINGLE
 http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__DOUBLE.html#group__CUDA__MATH__DOUBLE
*/

#define MACRO template <typename T> __device__


// In DTI model the parameters P are:
// S0: 0
// Dxx: 1
// Dxy: 2
// Dxz: 3 
// Dyy: 4 
// Dyz: 5
// Dzz: 6
// CFP[0:2] are bvecs 
// CFP[3] are bvals

#define EXAMPLE_CONSTANT 1.0f

MACRO T A_subroutine(T input){
  return input+0.0;
}

MACRO T Predicted_Signal(
			 int npar, // Number of Parameters to estimated
			 T* P, // Model estimated parameters, use P*
			 T* FP, // Fixed Parameters for each voxel
			 T* CFP) // Fixed Parameters common to all the voxels
{
  T pred_signal= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]);
  
  return A_subroutine<T>(pred_signal)*EXAMPLE_CONSTANT;
}

// Partial derivatives respect each model parameter
MACRO void Partial_Derivatives(
			       int npar, // Number of Parameters to estimated
			       T* P, // Model estimated parameters, use P*
			       T* FP, // Fixed Parameters for each voxel
			       T* CFP, // Fixed Parameters common to all the voxels
			       T* derivatives) // Derivative respect each model estimated parameter
{
  // df/dS0
  derivatives[0]= exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]);

  // df/dDxx
  derivatives[1]= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-CFP[3]*powf(CFP[0],2));


  // df/dDxy
  derivatives[2]= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-2*CFP[3]*CFP[0]*CFP[1]);

  // df/dDxz
  derivatives[3]= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-2*CFP[3]*CFP[0]*CFP[2]);

  // df/dDyy
  derivatives[4]= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-CFP[3]*powf(CFP[1],2));

  // df/dDyz
  derivatives[5]= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-2*CFP[3]*CFP[1]*CFP[2]);

  // df/dDzz
  derivatives[6]= P[0]*exp(-CFP[3]*powf(CFP[0],2)*P[1] 
	- CFP[3]*powf(CFP[1],2)*P[4] 
	- CFP[3]*powf(CFP[2],2)*P[6] 
	- 2*CFP[3]*CFP[0]*CFP[1]*P[2] 
	- 2*CFP[3]*CFP[0]*CFP[2]*P[3] 
	- 2*CFP[3]*CFP[1]*CFP[2]*P[5]) 
	* (-CFP[3]*powf(CFP[2],2));


}

