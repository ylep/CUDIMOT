#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
#define NPARAMS 11
// fiso, fintra, kappa, beta, th, ph, psi
// LR_fiso, LR_fintra, LR_kappa, LR_beta
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

