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
#include <TRandom2.h>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <strstream>

#include <TFile.h>
#include <TH2.h>
#include <TMath.h>

MomResParameters* Pair::mMomRes1 = 0;
MomResParameters* Pair::mMomRes2 = 0;

double Pair::mMomResScale = 0.0;


TRandom2* Pair::mRand = new TRandom2;

double* Pair::mom1 = (double*) malloc(sizeof(double) * 4);
double* Pair::mom2 = (double*) malloc(sizeof(double) * 4);
double* Pair::pos1 = (double*) malloc(sizeof(double) * 4);
double* Pair::pos2 = (double*) malloc(sizeof(double) * 4);

int Pair::mRelocated        = 0;
int Pair::mRelocateTable[8] = {0, 1, 2, 3, 4, 5, 6, 7};

// int Pair::mRelocateTable[1] = 1;
// int Pair::mRelocateTable[2] = 2;
// int Pair::mRelocateTable[3] = 3;
// int Pair::mRelocateTable[4] = 4;
// int Pair::mRelocateTable[5] = 5;
// int Pair::mRelocateTable[6] = 6;
// int Pair::mRelocateTable[7] = 7;


Pair::Pair() : mBad(0) {
  mRandVar = 0;
  //  mRelocated=0;
  //   for (int ti=0; ti<8; ti++)
  //     mRelocateTable[ti] = ti;
};

Pair& Pair::operator=(Pair& aPair) {
  // Only momentum is copied. Not position.
  for (int ti = 0; ti < 4; ti++) {
    mMom = aPair.mMom;
    // mX1[ti] = aPair.mX1[ti];
    // mX2[ti] = aPair.mX2[ti];
    mBad = aPair.mBad;
  }
  // mRandVar=0;
  //  mRelocated = aPair.mRelocated;
  //   for (int ti=0; ti<8; ti++)
  //     mRelocateTable[ti] = aPair.mRelocateTable[ti];

  mKStarSigned = aPair.mKStarSigned;

  mKStarOut  = aPair.mKStarOut;
  mKStarSide = aPair.mKStarSide;
  mKStarLong = aPair.mKStarLong;

  mBin = aPair.mBin;
  return *this;
}

Pair::~Pair() {
  // if(mMomRes1) delete mMomRes1;
  // if(mMomRes2)delete mMomRes2;
  if (mRandVar) delete[] mRandVar;
}

void Pair::SetMomentum(const float* aVect) {


  // ___________________________________________________
  // --- Calculate signed kStar and bin it belongs to
  VectorPair tMom;
  tMom.part1.x = aVect[mRelocateTable[0]];
  tMom.part1.y = aVect[mRelocateTable[1]];
  tMom.part1.z = aVect[mRelocateTable[2]];
  tMom.part1.t = aVect[mRelocateTable[3]];
  mMom.part1.x = aVect[mRelocateTable[0]];
  mMom.part1.y = aVect[mRelocateTable[1]];
  mMom.part1.z = aVect[mRelocateTable[2]];
  mMom.part1.t = aVect[mRelocateTable[3]];
  tMom.part2.x = aVect[mRelocateTable[4]];
  tMom.part2.y = aVect[mRelocateTable[5]];
  tMom.part2.z = aVect[mRelocateTable[6]];
  tMom.part2.t = aVect[mRelocateTable[7]];
  mMom.part2.x = aVect[mRelocateTable[4]];
  mMom.part2.y = aVect[mRelocateTable[5]];
  mMom.part2.z = aVect[mRelocateTable[6]];
  mMom.part2.t = aVect[mRelocateTable[7]];


  if (mMomRes1 && mMomRes2) {


    if (1) {

      double tP = mMom.part1.x * mMom.part1.x + mMom.part1.y * mMom.part1.y + mMom.part1.z * mMom.part1.z;
      tP        = sqrt(tP);

      double per     = mMomRes1->PResA + mMomRes1->PResB * pow(tP, mMomRes1->PResAlpha) + mMomRes1->PResC * tP;
      double thetaan = TMath::ATan2(hypot(mMom.part1.x, mMom.part1.y), mMom.part1.z);
      double phier   = 0.0;
      if (fabs(mMomRes1->PhiMu) < 0.00000001)
        phier = mMomRes1->PhiA + mMomRes1->PhiB * pow(tP, mMomRes1->PhiAlpha);
      else
        phier = mMomRes1->PhiA + mMomRes1->PhiB * exp(mMomRes1->PhiAlpha * (tP - mMomRes1->PhiMu));
      double thetaer = mMomRes1->ThetaA + mMomRes1->ThetaB * pow(tP, mMomRes1->ThetaAlpha);

      double tmass = TMath::Sqrt(mMom.part1.t * mMom.part1.t - tP * tP);

      double pshift  = (mMomRes1->PMeanA + mMomRes1->PMeanB * pow(1.0 + tmass * tmass / (tP * tP), mMomRes1->PMeanAlpha));
      double rescale = (tP - pshift * mMomResScale) / tP;

      double Deltapx = pow(TMath::Abs(mMom.part1.x) / tP * per, 2) + pow(TMath::Abs(mMom.part1.y) * phier, 2)
                       + pow(TMath::Abs(mMom.part1.x * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      double Deltapy = pow(TMath::Abs(mMom.part1.y) / tP * per, 2) + pow(TMath::Abs(mMom.part1.x) * phier, 2)
                       + pow(TMath::Abs(mMom.part1.y * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      double Deltapz =
        pow(TMath::Abs(mMom.part1.z) / tP * per, 2) + pow(TMath::Abs(mMom.part1.z * TMath::Tan(thetaan)) * thetaer, 2);

      Deltapx = sqrt(Deltapx);
      Deltapy = sqrt(Deltapy);
      Deltapz = sqrt(Deltapz);


      tMom.part1.x = ((mMom.part1.x + mRand->Gaus(0, fabs(Deltapx) * mMomResScale)) * rescale);
      tMom.part1.y = ((mMom.part1.y + mRand->Gaus(0, fabs(Deltapy) * mMomResScale)) * rescale);
      tMom.part1.z = ((mMom.part1.z + mRand->Gaus(0, fabs(Deltapz) * mMomResScale)) * rescale);
      tMom.part1.t =
        TMath::Sqrt(tMom.part1.x * tMom.part1.x + tMom.part1.y * tMom.part1.y + tMom.part1.z * tMom.part1.z + tmass * tmass);


      //      phian = atan2(tMom.part2.y,tMom.part2.x);
      tP = mMom.part2.x * mMom.part2.x + mMom.part2.y * mMom.part2.y + mMom.part2.z * mMom.part2.z;
      tP = sqrt(tP);

      per = mMomRes2->PResA + mMomRes2->PResB * pow(tP, mMomRes2->PResAlpha) + mMomRes2->PResC * tP;

      thetaan = TMath::ATan2(hypot(mMom.part2.x, mMom.part2.y), mMom.part2.z);
      if (fabs(mMomRes2->PhiMu) < 0.00000001)
        phier = mMomRes2->PhiA + mMomRes2->PhiB * pow(tP, mMomRes2->PhiAlpha);
      else
        phier = mMomRes2->PhiA + mMomRes2->PhiB * exp(mMomRes2->PhiAlpha * (tP - mMomRes2->PhiMu));
      thetaer = mMomRes2->ThetaA + mMomRes2->ThetaB * pow(tP, mMomRes2->ThetaAlpha);
      tmass   = TMath::Sqrt(mMom.part2.t * mMom.part2.t - tP * tP);

      //      pshift = (mMomRes2->PMeanA+mMomRes2->PMeanB*pow(tP,mMomRes2->PMeanAlpha));
      // New Energy loss from 0808.2041
      pshift  = (mMomRes2->PMeanA + mMomRes2->PMeanB * pow(1.0 + tmass * tmass / (tP * tP), mMomRes2->PMeanAlpha));
      rescale = (tP - pshift * mMomResScale) / tP;

      Deltapx = pow(TMath::Abs(mMom.part2.x) / tP * per, 2) + pow(TMath::Abs(mMom.part2.y) * phier, 2)
                + pow(TMath::Abs(tMom.part2.x * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      Deltapy = pow(TMath::Abs(mMom.part2.y) / tP * per, 2) + pow(TMath::Abs(mMom.part2.x) * phier, 2)
                + pow(TMath::Abs(mMom.part2.y * (1 / TMath::Tan(thetaan))) * thetaer, 2);
      Deltapz = pow(TMath::Abs(mMom.part2.z) / tP * per, 2) + pow(TMath::Abs(mMom.part2.z * TMath::Tan(thetaan)) * thetaer, 2);

      tMom.part2.x = ((mMom.part2.x + mRand->Gaus(0, fabs(Deltapx) * mMomResScale)) * rescale);
      tMom.part2.y = ((mMom.part2.y + mRand->Gaus(0, fabs(Deltapy) * mMomResScale)) * rescale);
      tMom.part2.z = ((mMom.part2.z + mRand->Gaus(0, fabs(Deltapz) * mMomResScale)) * rescale);
      tMom.part2.t =
        TMath::Sqrt(tMom.part2.x * tMom.part2.x + tMom.part2.y * tMom.part2.y + tMom.part2.z * tMom.part2.z + tmass * tmass);
    }
  }

  SmearMom.part1.x = tMom.part1.x;
  SmearMom.part1.y = tMom.part1.y;
  SmearMom.part1.z = tMom.part1.z;
  SmearMom.part1.t = tMom.part1.t;

  SmearMom.part2.x = tMom.part2.x;
  SmearMom.part2.y = tMom.part2.y;
  SmearMom.part2.z = tMom.part2.z;
  SmearMom.part2.t = tMom.part2.t;

  RawMom.part1.x = mMom.part1.x;
  RawMom.part1.y = mMom.part1.y;
  RawMom.part1.z = mMom.part1.z;
  RawMom.part1.t = mMom.part1.t;

  RawMom.part2.x = mMom.part2.x;
  RawMom.part2.y = mMom.part2.y;
  RawMom.part2.z = mMom.part2.z;
  RawMom.part2.t = mMom.part2.t;

  double tPx = tMom.part1.x + tMom.part2.x;
  double tPy = tMom.part1.y + tMom.part2.y;
  double tPz = tMom.part1.z + tMom.part2.z;
  double tE  = tMom.part1.t + tMom.part2.t;
  double tPt = tPx * tPx + tPy * tPy;
  // mCVK = tPz*tPz;
  double tMt = tE * tE - tPz * tPz;  // mCVK;
  // mCVK += tPt;
  // mCVK = sqrt(mCVK);
  double tM = sqrt(tMt - tPt);
  tMt       = sqrt(tMt);
  tPt       = sqrt(tPt);

  // Boost to LCMS
  double tBeta  = tPz / tE;
  double tGamma = tE / tMt;
  mKStarLong    = tGamma * (tMom.part1.z - tBeta * tMom.part1.t);
  double tE1L   = tGamma * (tMom.part1.t - tBeta * tMom.part1.z);

  // Rotate in transverse plane
  mKStarOut  = (tMom.part1.x * tPx + tMom.part1.y * tPy) / tPt;
  mKStarSide = (-tMom.part1.x * tPy + tMom.part1.y * tPx) / tPt;

  // Boost to pair cms
  mKStarOut = tMt / tM * (mKStarOut - tPt / tMt * tE1L);

  mKStarSigned = mKStarOut > 0. ? 1. : -1.;
  mKStarSigned *= sqrt(mKStarSide * mKStarSide + mKStarOut * mKStarOut + mKStarLong * mKStarLong);
}

void Pair::Set(const float* aVect) { SetMomentum(aVect); }

void Pair::InitRandVar() {
  //  PRINT_DEBUG_3("Initialized RandVar to: " << mRandVar[0] << " " << mRandVar[1] << " " << mRandVar[2] << " " << mRandVar[3]);
}


double Pair::GetBetat() {
  double tPx = mMom.part1.x + mMom.part2.x;
  double tPy = mMom.part1.y + mMom.part2.y;
  double tPz = mMom.part1.z + mMom.part2.z;
  double tE  = mMom.part1.t + mMom.part2.t;
  return sqrt((tPx * tPx + tPy * tPy) / (tE * tE - tPz * tPz));
}
double Pair::GetKt() {
  double tPx = mMom.part1.x + mMom.part2.x;
  double tPy = mMom.part1.y + mMom.part2.y;
  double tPz = mMom.part1.z + mMom.part2.z;
  return sqrt(tPx * tPx + tPy * tPy);
}

void Pair::SetMomRes(double aMomResScale, int aPairSystem) {
  mMomResScale = aMomResScale;

  switch (aPairSystem) {
    case 1:
    case 2:
    case 3:
    case 4:
      mMomRes1 = GetPionMomResParameters();
      mMomRes2 = GetKaonMomResParameters();
      break;
    case 11:
      mMomRes1 = GetPionMomResParameters();
      mMomRes2 = GetProtonMomResParameters();
      break;
    case 12:
      mMomRes1 = GetPionMomResParameters();
      mMomRes2 = GetAntiProtonMomResParameters();
      break;
    case 13:
      mMomRes1 = GetPionMomResParameters();
      mMomRes2 = GetProtonMomResParameters();
      break;
    case 14:
      mMomRes1 = GetPionMomResParameters();
      mMomRes2 = GetAntiProtonMomResParameters();
      break;
    case 21:
      mMomRes1 = GetKaonMomResParameters();
      mMomRes2 = GetProtonMomResParameters();
      break;
    case 22:
      mMomRes1 = GetKaonMomResParameters();
      mMomRes2 = GetAntiProtonMomResParameters();
      break;
    case 23:
      mMomRes1 = GetKaonMomResParameters();
      mMomRes2 = GetProtonMomResParameters();
      break;
    case 24:
      mMomRes1 = GetKaonMomResParameters();
      mMomRes2 = GetAntiProtonMomResParameters();
      break;
    case 31:
    case 32:
    case 33:
      mMomRes1 = GetPionMomResParameters();
      mMomRes2 = GetPionMomResParameters();
      break;
    case 41:
    case 42:
    case 43:
      mMomRes1 = GetKaonMomResParameters();
      mMomRes2 = GetKaonMomResParameters();
      break;
    case 51:
      mMomRes1 = GetProtonMomResParameters();
      mMomRes2 = GetProtonMomResParameters();
      break;
  }
}

double Pair::GetMomRes() { return mMomResScale; }


MomResParameters* Pair::GetPionMomResParameters() {
  MomResParameters* t = new MomResParameters;
  //   t->PhiA = 0.00204489;
  //   t->PhiB = 0.00100642;
  //   t->PhiAlpha = -1.28315;
  //   t->ThetaA = 0.000930104;
  //   t->ThetaB = 0.00123727;
  //   t->ThetaAlpha = -1.14998;
  //   t->PResA = 0.0113174;
  //   t->PResB = 0.00168459;
  //   t->PResAlpha = -1.04594;
  //   t->PResC = 0.00948437;
  //   t->PMeanA = 0.; // no energy loss (Kalman filter)
  //   t->PMeanB = 0.;
  //   t->PMeanAlpha = 0.;
  t->PhiA       = 0.00201;
  t->PhiB       = 0.001018;
  t->PhiAlpha   = -1.274;
  t->PhiMu      = 0.0;
  t->ThetaA     = 0.000908;
  t->ThetaB     = 0.001255;
  t->ThetaAlpha = -1.141;
  t->PResA      = 0.01074;
  t->PResB      = 0.001918;
  t->PResAlpha  = -0.9895;
  t->PResC      = 0.009706;
  t->PMeanA     = 0.;  // no energy loss (Kalman filter)
  t->PMeanB     = 0.;
  t->PMeanAlpha = 0.;
  return t;
}

MomResParameters* Pair::GetKaonMomResParameters() {
  MomResParameters* t = new MomResParameters;
  //   t->PResA = 0.8015;
  //   t->PResB = -2.362;
  //   t->PResC = 1.846;
  //   t->PResAlpha = 0.3445;
  //   t->PhiA = 0.0006575;
  //   t->PhiB = 0.002813;
  //   t->PhiAlpha = -1.583;
  //   t->ThetaA = 0.0002846;
  //   t->ThetaB = 0.002458;
  //   t->ThetaAlpha = -1.475;
  //   t->PMeanA = -0.006509;
  //   t->PMeanB = 0.008757;
  //   t->PMeanAlpha = -1.373;
  t->PResA     = 0.01981;
  t->PResB     = 0.001371;
  t->PResC     = 0.0;
  t->PResAlpha = -2.112;

  t->PhiA     = 0.001791;
  t->PhiB     = 0.001319;
  t->PhiAlpha = -1.686;
  t->PhiMu    = 0.0;

  t->ThetaA     = 0.0005202;
  t->ThetaB     = 0.001752;
  t->ThetaAlpha = -1.352;

  //   t->PMeanA = -0.004136;
  //   t->PMeanB = 0.003511;
  //   t->PMeanAlpha = -1.192;

  // New parametrization of energy loss from STAR 0808.2041
  t->PMeanA     = 0.006;
  t->PMeanB     = -0.0038;
  t->PMeanAlpha = 1.10;

  return t;
}

MomResParameters* Pair::GetProtonMomResParameters() {
  MomResParameters* t = new MomResParameters;
  //   t->PhiA = 0.00179782;
  //   t->PhiB = 0.00133481;
  //   t->PhiAlpha = -1.68111;
  //   t->ThetaA = 0.000495737;
  //   t->ThetaB = 0.0017755;
  //   t->ThetaAlpha = -1.34328;
  //   t->PResA = 0.0198917;
  //   t->PResB = 0.00145757;
  //   t->PResAlpha = -2.07794;
  //   t->PResC = 0;
  //   t->PMeanA = -0.00419502;
  //   t->PMeanB = 0.00353296;
  //   t->PMeanAlpha = -1.19465;

  //   t->PResA = 0.7989;
  //   t->PResB = -2.128;
  //   t->PResAlpha = 0.5282;
  //   t->PResC = 1.41;

  // Year 1 parameterization

  //   t->PResA = 0.01708;
  //   t->PResB = 0.006794;
  //   t->PResC = 0.0;
  //   t->PResAlpha = -1.78;

  //   t->PhiA = 0.0006575;
  //   t->PhiB = 0.002813;
  //   t->PhiAlpha = -1.583;

  //   t->ThetaA = 0.0002846;
  //   t->ThetaB = 0.002458;
  //   t->ThetaAlpha = -1.475;

  //   t->PMeanA = -0.006509;
  //   t->PMeanB = 0.008757;
  //   t->PMeanAlpha = -1.373;

  // Year 2 parameterization
  t->PResA     = -0.03516;
  t->PResB     = 0.005075;
  t->PResC     = 0.07254;
  t->PResAlpha = -2.015;

  t->PhiA     = 0.008905;
  t->PhiB     = 0.01229;
  t->PhiAlpha = -6.221;
  t->PhiMu    = 0.09998;

  t->ThetaA     = 0.006189;
  t->ThetaB     = 0.001732;
  t->ThetaAlpha = -1.981;

  //   t->PMeanA = -0.00905;
  //   t->PMeanB = 0.007839;
  //   t->PMeanAlpha = -1.69;

  // New parametrization of energy loss from STAR 0808.2041
  t->PMeanA     = 0.013;
  t->PMeanB     = -0.0081;
  t->PMeanAlpha = 1.03;

  return t;
}

MomResParameters* Pair::GetAntiProtonMomResParameters() {
  MomResParameters* t = new MomResParameters;
  // Year 2 parameterization
  t->PResA     = -0.01857;
  t->PResB     = 0.002276;
  t->PResC     = 0.05261;
  t->PResAlpha = -2.545;

  t->PhiA     = 0.008614;
  t->PhiB     = 0.0174;
  t->PhiAlpha = -5.717;
  t->PhiMu    = 0.3319;

  t->ThetaA     = 0.002326;
  t->ThetaB     = 0.00596;
  t->ThetaAlpha = -1.052;

  t->PMeanA     = -0.006365;
  t->PMeanB     = 0.006887;
  t->PMeanAlpha = -1.747;

  return t;
}


void Pair::SetPosMom() {
  mom1[0] = mMom.part1.x;
  mom1[1] = mMom.part1.y;
  mom1[2] = mMom.part1.z;
  mom1[3] = mMom.part1.t;

  mom2[0] = mMom.part2.x;
  mom2[1] = mMom.part2.y;
  mom2[2] = mMom.part2.z;
  mom2[3] = mMom.part2.t;

  pos1[0] = mPos.part1.x;
  pos1[1] = mPos.part1.y;
  pos1[2] = mPos.part1.z;
  pos1[3] = mPos.part1.t;

  pos2[0] = mPos.part2.x;
  pos2[1] = mPos.part2.y;
  pos2[2] = mPos.part2.z;
  pos2[3] = mPos.part2.t;
}

void Pair::SetPosMomEFirstNotation() {
  mom1[1] = mMom.part1.x * 1000.0;
  mom1[2] = mMom.part1.y * 1000.0;
  mom1[3] = mMom.part1.z * 1000.0;
  mom1[0] = mMom.part1.t * 1000.0;

  mom2[1] = mMom.part2.x * 1000.0;
  mom2[2] = mMom.part2.y * 1000.0;
  mom2[3] = mMom.part2.z * 1000.0;
  mom2[0] = mMom.part2.t * 1000.0;

  pos1[1] = mPos.part1.x;
  pos1[2] = mPos.part1.y;
  pos1[3] = mPos.part1.z;
  pos1[0] = mPos.part1.t;

  pos2[1] = mPos.part2.x;
  pos2[2] = mPos.part2.y;
  pos2[3] = mPos.part2.z;
  pos2[0] = mPos.part2.t;
}

void Pair::SetPairNum(int aNum) { mPairNum = aNum; }

int Pair::GetPairNum() { return mPairNum; }


void Pair::GetSeed(UInt_t& s1, UInt_t& s2) { s1 = mRand->GetSeed(); }

void Pair::SetSeed(UInt_t s1, UInt_t s2) { mRand->SetSeed(s1); }
