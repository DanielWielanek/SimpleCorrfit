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
#define SPPID 3222
#define SZPID 3212
#define LAPID 3122
#define XIZERPID 3322
#define XIMINPID 3312
#define NEPID 2112
#define ABSRAP 0.7
#define ETAABS 0.35
#define PTMIN 0.12
#define PTMAX 1.5
#define TMAX 500

#define COULOMBSTEPS 170

#define KAONMASS 0.493677
#define PIONMASS 0.13957

#define FOURPI 1.25663706143591725e+01


struct _dcomplex {
  long double re;
  long double im;
};

typedef struct _dcomplex dcomplex;


class StandAloneSimpleFsi {

  double mKStarLong, mKStarOut, mKStarSide, mKStarSigned, mBetat;
  double pairphi;

  double deta, dphi;

  double mROS, mRSS, mRLS, mRSt;
  double mKT;
  double ktmin, ktmax;

  int pairtype;
  int partpid;
  int partpid2;

  int CounterCosPhi;
  int allCounterCosPhi;


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
  void SelectPairType() {};
  void InitializeGamow();
  // int IsResidualPair(int,int,int,int);
  int IsIdentical(int) { return 0; };
  double getKStarSigned() { return mKStarSigned; }
  virtual void setQuantumOn() {};
  virtual void setStrongOff() {};
  virtual void setCoulOn() {};
  virtual double getWeight(Pair& aPair);
  virtual void setDefaultCalcPar() {};
  virtual void setStrongOn() {};
  virtual void setPairType(int aPairType);
  virtual void setQuantumOff() {};
  virtual void set3BodyOn() {};
  virtual void setCoulOff() {};
  virtual void set3BodyOff() {};
  virtual int getPairType();
  inline double GetRout() const { return mROS * 0.197326968; };
  inline double GetRside() const { return mRSS * 0.197326968; };
  inline double GetRlong() const { return mRLS * 0.197326968; };
  virtual ~StandAloneSimpleFsi();
};


#endif /* CORRFIT_CODE_STANDALONESIMPLEFSI_H_ */
