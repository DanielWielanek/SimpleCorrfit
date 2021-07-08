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
  double p1x = aPair.p1().x;
  double p1y = aPair.p1().y;
  double p1z = aPair.p1().z;
  double p1e = aPair.p1().t;


  double p2x = aPair.p2().x;
  double p2y = aPair.p2().y;
  double p2z = aPair.p2().z;
  double p2e = aPair.p2().t;

  double tPx = p1x + p2x;
  double tPy = p1y + p2y;
  double tPz = p1z + p2z;
  double tE  = p1e + p2e;
  double tPt = tPx * tPx + tPy * tPy;
  double tMt = tE * tE - tPz * tPz;  // mCVK;
  double tM  = sqrt(tMt - tPt);
  tMt        = sqrt(tMt);
  tPt        = sqrt(tPt);

  mKT     = tPt / 2.0;
  pairphi = TMath::ATan2(tPy, tPx);

  //  if ((mKT < ktmin) || (mKT > ktmax)) return 0;

  mBetat = tPt / tMt;

  // Boost to LCMS
  double tBeta  = tPz / tE;
  double tGamma = tE / tMt;
  mKStarLong    = tGamma * (p1z - tBeta * p1e);
  double tE1Lm  = tGamma * (p1e - tBeta * p1z);


  // Rotate in transverse plane
  mKStarOut  = (p1x * tPx + p1y * tPy) / tPt;
  mKStarSide = (-p1x * tPy + p1y * tPx) / tPt;

  // Boost to pair cms
  mKStarOut = tMt / tM * (mKStarOut - tPt / tMt * tE1Lm);


  Double_t tGammat = 1.0 / TMath::Sqrt(1.0 - mBetat * mBetat);


  // mKO = particle1lcmsm.fpx - particle2lcmsm.fpx;

  mKStarSigned = mKStarOut > 0. ? 1. : -1.;
  mKStarSigned *= sqrt(mKStarSide * mKStarSide + mKStarOut * mKStarOut + mKStarLong * mKStarLong);

  const double mFmToGeV = 0.197326968;


  double tDX    = aPair.x1().x - aPair.x2().x;
  double tDY    = aPair.x1().y - aPair.x2().y;
  double tRLong = aPair.x1().z - aPair.x2().z;
  double tDTime = aPair.x1().t - aPair.x2().t;

  /* double tDX    = aPair.pos1[0] - aPair.pos2[0];
   double tDY    = aPair.pos1[1] - aPair.pos2[1];
   double tRLong = aPair.pos1[2] - aPair.pos2[2];
   double tDTime = aPair.pos1[3] - aPair.pos2[3];
 */
  double tROut         = (tDX * tPx + tDY * tPy) / tPt;
  double tRSide        = (-tDX * tPy + tDY * tPx) / tPt;
  double tRSidePairCMS = tRSide;
  mRSS                 = tRSidePairCMS / mFmToGeV;

  double tRLongPairCMS  = tGamma * (tRLong - tBeta * tDTime);
  mRLS                  = tRLongPairCMS / mFmToGeV;
  double tDTimePairLCMS = tGamma * (tDTime - tBeta * tRLong);

  tBeta                = tPt / tMt;
  tGamma               = tMt / tM;
  double tROutPairCMS  = tGamma * (tROut - tBeta * tDTimePairLCMS);
  mROS                 = tROutPairCMS / mFmToGeV;
  double tDTimePairCMS = tGamma * (tDTimePairLCMS - tBeta * tROut);
  // mRTS                 = tDTimePairCMS / mFmToGeV;
  double tRStar = ::sqrt(tROutPairCMS * tROutPairCMS + tRSidePairCMS * tRSidePairCMS + tRLongPairCMS * tRLongPairCMS);
  mRSt          = tRStar / mFmToGeV;
  // std::cout << "KSTARS" << mKStarOut << " " << mKStarSide << " " << mKStarLong << std::endl;  // t jest ok
  //  std::cout << "RSTARS" << mRSS << " " << mROS << " " << mRLS << " " << mRSt << std::endl;
}

StandAloneSimpleFsi::~StandAloneSimpleFsi() {
  // TODO Auto-generated destructor stub
}
