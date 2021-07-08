/***************************************************************************
 *
 * $Id:
 *
 * Author: Laurent Conin, Fabrice Retiere, Subatech, France
 ***************************************************************************
 *
 * Description : Create pair from StHbtFsiHiddenInfo
 *
 ***************************************************************************
 *
 * $Log:
 *
 ***************************************************************************/

#ifndef ST_HBT_DATA_PAIR_HH
#define ST_HBT_DATA_PAIR_HH

#include <iostream>
#include <math.h>


struct _FourVector {
  double x, y, z, t;
};

typedef struct _FourVector FourVector;

struct _VectorPair {
  FourVector part1, part2;
};

typedef struct _VectorPair VectorPair;

struct MomResParameters {
  double PhiA;
  double PhiB;
  double PhiAlpha;
  double PhiMu;
  double ThetaA;
  double ThetaB;
  double ThetaAlpha;
  double PMeanA;
  double PMeanB;
  double PMeanAlpha;
  double PResA;
  double PResB;
  double PResAlpha;
  double PResC;
};

class Pair {
public:
  Pair();
  ~Pair();
  Pair& operator=(Pair& aPair);

  void SetMomentum(const float* aVect);
  void Set(const float* aVect);

  static void SetMomRes(double aMomResScale, int aPairSystem);
  static double GetMomRes();

  void SetPairNum(int aNum);
  int GetPairNum();

  static int nPar();


  double GetKStarTimeOutSign();
  double GetKStarTimeSideSign();
  double GetKStarTimeLongSign();
  double GetKStar();
  double GetBetat();
  double GetKt();

  double GetKStarTimeOutSignSmeared();
  double GetKStarTimeSideSignSmeared();
  double GetKStarTimeLongSignSmeared();
  double GetKStarSmeared();

  double GetKStarOut();
  double GetKStarSide();
  double GetKStarLong();

  FourVector p1();
  FourVector p2();
  FourVector x1();
  FourVector x2();

  VectorPair PairMomentum();
  VectorPair PairPosition();

  int bad() { return mBad; }
  int GetBin();
  void SetBin(int aBin);

  void SetPosMom();
  void SetPosMomEFirstNotation();
  void InitRandVar();

  static double* mom1;
  static double* mom2;
  static double* pos1;
  static double* pos2;

  VectorPair mMom;
  VectorPair mPos;


  VectorPair RawMom;
  VectorPair SmearMom;

private:
  int mBad;

  static int mRelocateTable[8];
  static int mRelocated;
  static double mMomResScale;

  double* mRandVar;

  double mKStarSigned;

  double mKStarOut;
  double mKStarSide;
  double mKStarLong;

  int mBin;

  int mPairNum;

  static MomResParameters* GetPionMomResParameters();
  static MomResParameters* GetKaonMomResParameters();
  static MomResParameters* GetProtonMomResParameters();
  static MomResParameters* GetAntiProtonMomResParameters();
  static MomResParameters* mMomRes1;
  static MomResParameters* mMomRes2;
};

inline int Pair::GetBin() { return mBin; }

inline void Pair::SetBin(int aBin) { mBin = aBin; }

inline double Pair::GetKStarTimeOutSign() { return mKStarOut; }
inline double Pair::GetKStarTimeSideSign() { return mKStarSide; }
inline double Pair::GetKStarTimeLongSign() { return mKStarLong; }
inline double Pair::GetKStar() { return mKStarSigned > 0.0 ? mKStarSigned : -mKStarSigned; }

inline double Pair::GetKStarOut() { return mKStarOut; };
inline double Pair::GetKStarSide() { return mKStarSide; };
inline double Pair::GetKStarLong() { return mKStarLong; };

inline FourVector Pair::p1() { return mMom.part1; }
inline FourVector Pair::p2() { return mMom.part2; }
inline FourVector Pair::x1() { return mPos.part1; }
inline FourVector Pair::x2() { return mPos.part2; }


inline VectorPair Pair::PairMomentum() { return mMom; }
inline VectorPair Pair::PairPosition() { return mPos; }

#endif
