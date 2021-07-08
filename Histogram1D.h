/*
 * Histogram1D.h
 *
 *  Created on: 8 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#ifndef SIMPLECORRFIT_HISTOGRAM1D_H_
#define SIMPLECORRFIT_HISTOGRAM1D_H_

#include <string>

class Histogram1D {
  int mBinsNo;
  double mBins[100];
  double mMin;
  double mMax;
  double mWidth;

public:
  Histogram1D(double min, double max);
  double GetMin() const { return mMin; };
  double GetMax() const { return mMax; };
  int GetNBins() const { return mBinsNo; };
  void Fill(double val, double weight);
  double GetBinContent(int i) const { return mBins[i]; }
  ~Histogram1D();
};


#endif /* SIMPLECORRFIT_HISTOGRAM1D_H_ */
