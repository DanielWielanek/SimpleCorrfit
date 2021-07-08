/*
 * StandAloneSimpleFsi.cxx
 *
 *  Created on: 18 cze 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "StandAloneSimpleFsi.h"
#include <TMath.h>

using namespace std;

StandAloneSimpleFsi::StandAloneSimpleFsi() {
  // TODO Auto-generated constructor stub
}

double StandAloneSimpleFsi::GetCoulomb() {
  double kstrst = mKStarOut * mROS + mKStarSide * mRSS + mKStarLong * mRLS;

  // Classical limit - if distance is larger than Coulomb radius,
  // the interaction does not matter
  if (fabs(mRSt) > fabs(pionac)) return (1.0);
  if (fabs(mRSt) == 0.0) return (Gamow(fabs(mKStarSigned)));

  // Classical limit - in the case of large k*r* product we go to
  // classical coulomb interaction
  if (fabs(mKStarSigned) * mRSt * (1.0 + kstrst / (mRSt * fabs(mKStarSigned))) > 15.0)
    return (1.0 - 1.0 / (mRSt * pionac * mKStarSigned * mKStarSigned));

  // Calculate the F function
  dcomplex ffplus;
  GetFFsingle(&ffplus);

  if (!finite(ffplus.re)) {
    cout << "FFPlus Re not a number !"
         << " " << endl;
    cout << mRSt << " " << mKStarSigned << " " << mROS << " " << mRSS << " " << mRLS << endl;
  }

  if (!finite(ffplus.im))
    cout << "FFPlus Im not a number !"
         << " " << endl;


  return (Gamow(fabs(mKStarSigned)) * modl2(ffplus));
}

long double StandAloneSimpleFsi::modl2(dcomplex arg) { return arg.re * arg.re + arg.im * arg.im; }

void StandAloneSimpleFsi::GetFFsingle(dcomplex* ffp, int sign) {
  double comprep[COULOMBSTEPS];
  double compimp[COULOMBSTEPS];
  double eta, ksip;
  dcomplex alfa, zetp;

  int nsteps;

  double kstar  = fabs(mKStarSigned);
  double kstrst = mKStarOut * mROS + mKStarSide * mRSS + mKStarLong * mRLS;
  double coskr  = sign * kstrst / (fabs(mKStarSigned) * mRSt);

  if (kstar * mRSt * (1 + coskr) > 15.0)
    nsteps = 170;
  else if (kstar * mRSt * (1 + coskr) > 10.0)
    nsteps = 45;
  else if (kstar * mRSt * (1 + coskr) > 5.0)
    nsteps = 35;
  else
    nsteps = 15;

  eta     = 1.0 / (kstar * pionac);
  alfa.re = 0.0;
  alfa.im = -eta;

  dcomplex fcomp, scompp;
  double tcomp;
  dcomplex sump;
  dcomplex fcmult;

  double rad = mRSt;

  ksip = kstar * rad * (1 + coskr);

  zetp.re = 0.0;
  zetp.im = ksip;

  fcomp.re  = 1.0;
  fcomp.im  = 0.0;
  scompp.re = 1.0;
  scompp.im = 0.0;
  tcomp     = 1.0;

  for (int istep = 0; istep < nsteps; istep++) {
    sump = mult(fcomp, scompp);

    sump = mult(sump, 1.0 / (tcomp * tcomp));

    if (istep == 0) {
      comprep[istep] = sump.re;
      compimp[istep] = sump.im;
    } else {
      comprep[istep] = comprep[istep - 1] + sump.re;
      compimp[istep] = compimp[istep - 1] + sump.im;
    }

    fcmult.re = alfa.re + istep;
    fcmult.im = alfa.im;

    fcomp  = mult(fcomp, fcmult);
    scompp = mult(scompp, zetp);
    tcomp *= (istep + 1);

    if ((sump.re * sump.re + sump.im * sump.im) < 1.0e-14) {
      nsteps = istep;
      break;
    }
  }

  ffp->re = comprep[nsteps - 1];
  ffp->im = compimp[nsteps - 1];
}

dcomplex StandAloneSimpleFsi::mult(dcomplex arga, long double argb) {
  dcomplex res;

  res.re = arga.re * argb;
  res.im = arga.im * argb;

  return res;
}

dcomplex StandAloneSimpleFsi::mult(dcomplex arga, dcomplex argb) {
  dcomplex res;

  res.re = arga.re * argb.re - arga.im * argb.im;
  res.im = arga.re * argb.im + argb.re * arga.im;

  return res;
}

double StandAloneSimpleFsi::Gamow(double arg) {
  //   if (arg<0.0001)
  //     return 0.0;
  //   if (arg > 0.4)
  long double eta = twopioverac / arg;
  return (eta) *1.0 / (exp(eta) - 1.0);
  //   int bin = arg / 0.0002;
  //   double b1 = bin*0.0002 + 0.0001;
  //   double b2 = bin*0.0002 + 0.0003;
  //   return ((gamov[bin+1] * (arg - b1) + gamov[bin] * (b2 - arg)) / 0.0002);
}

void StandAloneSimpleFsi::InitializeGamow() {
  twospin = 0;

  euler = 0.577215665;
  f0    = 7.77 / 0.197327;
  d0    = 2.77 / 0.197327;

  d0s.re = 0.0;
  d0s.im = 0.0;
  d0t.re = 0.0;
  d0t.im = 0.0;

  f0s.re = 0.0;
  f0s.im = 0.0;
  f0t.re = 0.0;
  f0t.im = 0.0;
  switch (pairtype) {
    case 2: {
      pionac    = 248.519 / 0.197327;
      writegrps = 0;
    } break;
    case 3: {
      pionac    = -248.519 / 0.197327;
      writegrps = 0;
    } break;
  }
  oneoveracsq = 1.0 / (pionac * pionac);
  twopioverac = 2.0 * TMath::Pi() / pionac;
  double tpaoverk;
  for (int iter = 0; iter < 2000; iter++) {
    tpaoverk    = twopioverac / (iter * 0.0002 + 0.0001);
    gamov[iter] = tpaoverk * 1.0 / (exp(tpaoverk) - 1);
  }
}

double StandAloneSimpleFsi::getWeight(Pair& aPair) {
  aPair.CalculateKinematics();
  Kinematics(aPair);
  InitializeGamow();
  return GetCoulomb();
}

void StandAloneSimpleFsi::setPairType(int aPairType) {
  if (aPairType == 11) {  // pi+K+
    pionac    = 248.519 / 0.197327;
    partpid   = PIPID;
    partpid2  = KPID;
    writegrps = 0;
    pairtype  = 2;
  } else if (aPairType == 10) {  // pi+K-
    pionac    = -248.519 / 0.197327;
    partpid   = -PIPID;
    partpid2  = KPID;
    writegrps = 0;
    pairtype  = 3;
  }
}

int StandAloneSimpleFsi::getPairType() { return pairtype; }

void StandAloneSimpleFsi::Kinematics(Pair& aPair) {
  mKStarOut    = aPair.GetKStarOut();
  mKStarSide   = aPair.GetKStarSide();
  mKStarLong   = aPair.GetKStarLong();
  mKStarSigned = aPair.GetKStarSigned();
  mRSS         = aPair.GetRSide();
  mROS         = aPair.GetROut();
  mRLS         = aPair.GetRLong();
  mRSt         = aPair.GetRStar();
}

StandAloneSimpleFsi::~StandAloneSimpleFsi() {
  // TODO Auto-generated destructor stub
}
