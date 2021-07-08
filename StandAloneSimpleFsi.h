/*
 * StandAloneSimpleFsi.h
 *
 *  Created on: 18 cze 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef CORRFIT_CODE_STANDALONESIMPLEFSI_H_
#define CORRFIT_CODE_STANDALONESIMPLEFSI_H_

#include "Pair.h"


#define PIPID 211
#define KPID 321
#define PPID 2212

#define COULOMBSTEPS 170

#define KAONMASS 0.493677
#define PIONMASS 0.13957


struct _dcomplex {
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;


class StandAloneSimpleFsi {

  double mKStarLong, mKStarOut, mKStarSide, mKStarSigned;
  double mROS, mRSS, mRLS, mRSt;

  int pairtype;
  int partpid;
  int partpid2;

  // end external variables


  dcomplex d0s;
  dcomplex f0s;
  dcomplex d0t;
  dcomplex f0t;

  double gamov[2000];
  long double pionac;
  long double oneoveracsq;
  long double twopioverac;
  long double coulqscpart;
  int twospin;
  int writegrps;
  long double euler;
  long double f0;
  long double d0;

  double GetCoulomb();
  void GetFFsingle(dcomplex* ffp, int sign = 1);
  dcomplex mult(dcomplex arga, long double argb);
  dcomplex mult(dcomplex arga, dcomplex argb);
  double Gamow(double arg);
  long double modl2(dcomplex arg);
  void Kinematics(Pair& pair);

public:
  StandAloneSimpleFsi();
  void InitializeGamow();
  double getKStarSigned() { return mKStarSigned; }
  double getWeight(Pair& aPair);
  void setPairType(int aPairType);
  int getPairType();
  inline double GetRout() const { return mROS * 0.197326968; };
  inline double GetRside() const { return mRSS * 0.197326968; };
  inline double GetRlong() const { return mRLS * 0.197326968; };
  ~StandAloneSimpleFsi();
};


#endif /* CORRFIT_CODE_STANDALONESIMPLEFSI_H_ */
