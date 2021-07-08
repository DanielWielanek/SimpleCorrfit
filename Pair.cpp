/***************************************************************************
 *
 * $Id:
 *
 * Author: Laurent Conin, Fabrice Retiere, Subatech, France
 ***************************************************************************
 *
 * Description : Create pair from HiddenInfo
 *
 ***************************************************************************
 *
 * $Log:
 *
 ***************************************************************************/
#include "Pair.h"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <strstream>


Pair::Pair() {};

Pair& Pair::operator=(Pair& aPair) {
  // Only momentum is copied. Not position.
  for (int ti = 0; ti < 4; ti++) {
    mMom = aPair.mMom;
    // mX1[ti] = aPair.mX1[ti];
    // mX2[ti] = aPair.mX2[ti];
  }

  mKStarSigned = aPair.mKStarSigned;

  mKStarOut  = aPair.mKStarOut;
  mKStarSide = aPair.mKStarSide;
  mKStarLong = aPair.mKStarLong;
  mRLS       = aPair.mRLS;
  mROS       = aPair.mROS;
  mRSS       = aPair.mRSS;
  mRSt       = aPair.mRSt;

  return *this;
}

Pair::~Pair() {}

void Pair::CalculateKinematics() {
  double p1x = p1().x;
  double p1y = p1().y;
  double p1z = p1().z;
  double p1e = p1().t;


  double p2x = p2().x;
  double p2y = p2().y;
  double p2z = p2().z;
  double p2e = p2().t;

  double tPx = p1x + p2x;
  double tPy = p1y + p2y;
  double tPz = p1z + p2z;
  double tE  = p1e + p2e;
  double tPt = tPx * tPx + tPy * tPy;
  double tMt = tE * tE - tPz * tPz;  // mCVK;
  double tM  = sqrt(tMt - tPt);
  tMt        = sqrt(tMt);
  tPt        = sqrt(tPt);

  mKT = tPt / 2.0;

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


  double tGammat = 1.0 / sqrt(1.0 - mBetat * mBetat);


  // mKO = particle1lcmsm.fpx - particle2lcmsm.fpx;

  mKStarSigned = mKStarOut > 0. ? 1. : -1.;
  mKStarSigned *= sqrt(mKStarSide * mKStarSide + mKStarOut * mKStarOut + mKStarLong * mKStarLong);

  const double mFmToGeV = 0.197326968;


  double tDX    = x1().x - x2().x;
  double tDY    = x1().y - x2().y;
  double tRLong = x1().z - x2().z;
  double tDTime = x1().t - x2().t;

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
}
