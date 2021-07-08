/*
 * Histogram1D.cpp
 *
 *  Created on: 8 lip 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "Histogram1D.h"

#include <iostream>

Histogram1D::Histogram1D(double min, double max) {
  mBinsNo = 100;
  mMin    = min;
  mMax    = max;
  mWidth  = (mMax - mMin) / ((double) (mBinsNo));
  mWidth  = 1.0 / mWidth;
  for (int i = 0; i < 100; i++)
    mBins[i] = 0;
}

void Histogram1D::Fill(double val, double weight) {
  if (val < mMin) return;
  if (val >= mMax) return;
  int bin = (val - mMin) * mWidth;
  mBins[bin] += weight;
}

Histogram1D::~Histogram1D() {}
