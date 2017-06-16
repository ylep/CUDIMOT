#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
#define NPARAMS 20
// fiso, fintra, ffibre2 
// kappa1, beta1, th1, ph1, psi1
// kappa2, beta2, th2, ph2, psi2
// LR_fiso, LR_fintra, LR_ffibre2 
// LR_kappa1, LR_beta1
// LR_kappa2, LR_beta2
#define NCFP 3 // bvecs, bvals, LR_map
#define NFIXP 1 // S0
/////////////////////

///// Do not edit this /////
struct MODEL
{
  static int CFP_size[NCFP];
  static int FixP_size[NFIXP];
};
////////////////////////////

#endif

