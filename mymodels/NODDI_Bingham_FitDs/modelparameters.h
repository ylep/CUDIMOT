#ifndef MODEL_H_INCLUDED
#define MODEL_H_INCLUDED

///// Edit this /////
typedef double MyType;
#define NPARAMS 9 // fiso, fintra, kappa, beta, th, ph, psi, Diso, Dparallel
#define NCFP 2 // bvecs, bvals
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

