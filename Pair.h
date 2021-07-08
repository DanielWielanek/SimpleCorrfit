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

class Pair {
public:
  Pair();
  ~Pair();
  Pair& operator=(Pair& aPair);
  void SetMomentum(const FourVector& p1, const FourVector& p2);
  void SetFreezout(const FourVector& x1, const FourVector& x2);

  inline double GetKStarOut() const { return mKStarOut; };
  inline double GetKStarSide() const { return mKStarSide; };
  inline double GetKStarLong() const { return mKStarLong; };
  inline double GetKStarSigned() const { return mKStarSigned; }
  inline double GetROut() const { return mROS; }
  inline double GetRSide() const { return mRSS; }
  inline double GetRLong() const { return mRLS; }
  inline double GetRStar() const { return mRSt; }
  void CalculateKinematics();

  FourVector p1();
  FourVector p2();
  FourVector x1();
  FourVector x2();

  VectorPair PairMomentum();
  VectorPair PairPosition();

  VectorPair mMom;
  VectorPair mPos;

private:
  double mKStarSigned;

  double mKStarOut;
  double mKStarSide;
  double mKStarLong;
  double mROS;
  double mRSS;
  double mRLS;
  double mRSt;

  double mKT;
  double mBetat;
};


inline FourVector Pair::p1() { return mMom.part1; }
inline FourVector Pair::p2() { return mMom.part2; }
inline FourVector Pair::x1() { return mPos.part1; }
inline FourVector Pair::x2() { return mPos.part2; }


inline VectorPair Pair::PairMomentum() { return mMom; }

inline void Pair::SetMomentum(const FourVector& p1, const FourVector& p2) {
  mMom.part1 = p1;
  mMom.part2 = p2;
}

inline void Pair::SetFreezout(const FourVector& x1, const FourVector& x2) {
  mPos.part1 = x1;
  mPos.part2 = x2;
}


inline VectorPair Pair::PairPosition() { return mPos; }

#endif
