/*
 * CorrFitSimple.cxx
 *
 *  Created on: 19 cze 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "StandAloneFsiLednicky.h"
#include "StandAloneSimpleFsi.h"
#include "TRandom.h"
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <iostream>
using namespace std;

const double gFmToGeVToFm = 5.067731036;
const double gFmToGeV     = 0.197326968;

struct Lorentz {
  double_t x, y, z, t;
};

Lorentz GetMom(double mass) {
  Lorentz l;
  l.x      = gRandom->Gaus(0, 0.1);
  l.y      = gRandom->Gaus(0, 0.1);
  l.z      = gRandom->Gaus(0, 0.1);
  double p = sqrt(l.x * l.x + l.y * l.y + l.z * l.z + mass * mass);
  l.t      = p;
  return l;
}

void GetPos(double* p) {
  p[0] = gRandom->Gaus(0, 1);
  p[1] = gRandom->Gaus(0, 1);
  p[2] = gRandom->Gaus(0, 1);
  p[3] = 0;  // gRandom->Gaus(0, 1);
}


void Boost(Pair* p) {

  double tPx = p->p1().x + p->p2().x;
  double tPy = p->p1().y + p->p2().y;
  double tPz = p->p1().z + p->p2().z;
  double tE  = p->p1().t + p->p2().t;
  double tPt = tPx * tPx + tPy * tPy;
  double tMt = tE * tE - tPz * tPz;
  double tM  = sqrt(tMt - tPt);

  double gammat = 1.0 / TMath::Sqrt(1.0 - tPt / tMt);

  tMt = sqrt(tMt);
  tPt = sqrt(tPt);

  double betat = tPt / tMt;

  double betaz  = tPz / tE;
  double gammaz = 1.0 / TMath::Sqrt(1.0 - betaz * betaz);

  double tROutS, tRSideS, tRLongS, tRTimeS;
  double tROut, tRSide, tRLong;

  tROutS  = p->x1().x * gFmToGeVToFm;
  tRSideS = p->x1().y * gFmToGeVToFm;
  tRLongS = p->x1().z * gFmToGeVToFm;

  tRTimeS = 0.0;

  tROut       = gammat * (tROutS + betat * tRTimeS);
  double tDtL = gammat * (tRTimeS + betat * tROutS);

  tRLong     = gammaz * (tRLongS + betaz * tDtL);
  double tDt = gammaz * (tDtL + betaz * tRLongS);

  tPx /= tPt;
  tPy /= tPt;

  p->mPos.part1.x = (tROut * tPx - tRSideS * tPy) * gFmToGeV;
  p->mPos.part1.y = (tROut * tPy + tRSideS * tPx) * gFmToGeV;
  p->mPos.part1.z = tRLong * gFmToGeV;
  p->mPos.part1.t = tDt * gFmToGeV;
}

TFile* sOutFile;


int main(int argc, char** argv) {


  TString systemName = "";

  StandAloneFsiLednicky* tWeightLed  = new StandAloneFsiLednicky();
  StandAloneSimpleFsi* tWeightKisiel = new StandAloneSimpleFsi();  // changed ! DW
  const Int_t kPiOpos                = 10;
  const Int_t kPiSame                = 11;
  const Int_t theChosenOne           = kPiOpos;

  if (theChosenOne == kPiOpos) {
    systemName = "opos.root";
  } else {
    systemName = "same.root";
  }


  tWeightLed->setCoulOn();
  tWeightLed->setQuantumOff();
  tWeightLed->setStrongOff();
  tWeightLed->set3BodyOff();
  tWeightLed->setPairType(theChosenOne);  // oposite

  tWeightKisiel->setCoulOn();
  tWeightKisiel->setQuantumOff();
  tWeightKisiel->setStrongOff();
  tWeightKisiel->setPairType(theChosenOne);

  Pair* pair = new Pair();

  TH1D* numL = new TH1D("numL", "numL", 100, 0, 1);
  TH1D* denL = new TH1D("denL", "denL", 100, 0, 1);
  TH1D* numK = new TH1D("numK", "numK", 100, 0, 1);
  TH1D* denK = new TH1D("denK", "denK", 100, 0, 1);

  double poz1[4], poz2[4];
  Int_t pairs  = 100000;
  Int_t kPairs = pairs / 1000;
  for (int i = 0; i < pairs; i++) {
    if (i % 1000 == 0) {
      cout << i / 1000 << " of " << kPairs << endl;
      ;
    }
    Lorentz p1         = GetMom(0.139);
    Lorentz p2         = GetMom(0.493);
    pair->mMom.part1.x = p1.x;
    pair->mMom.part1.y = p1.y;
    pair->mMom.part1.z = p1.z;
    pair->mMom.part1.t = p1.t;

    pair->mMom.part2.x = p2.x;
    pair->mMom.part2.y = p2.y;
    pair->mMom.part2.z = p2.z;
    pair->mMom.part2.t = p2.t;


    GetPos(poz1);

    pair->mPos.part1.x = poz1[0];
    pair->mPos.part1.y = poz1[1];
    pair->mPos.part1.z = poz1[2];
    pair->mPos.part1.t = poz1[3];

    Boost(pair);


    double wK = tWeightKisiel->getWeight(*pair);
    double wL = tWeightLed->getWeight(*pair);


    double kstar = tWeightKisiel->getKStarSigned();

    // cout << "REQ" << endl;
    // cout << poz1[0] << " " << poz1[1] << " " << poz1[2] << " " << poz1[3] << endl;
    // cout << tWeightKisiel->GetRout() << " " << tWeightKisiel->GetRside() << " " << tWeightKisiel->GetRlong() << endl;
    // cout << wK << " " << wL << "----------- " << kstar << endl;

    numL->Fill(kstar, wL);
    numK->Fill(kstar, wK);
    denL->Fill(kstar);
    denK->Fill(kstar);
  }

  TFile* fx = new TFile(systemName, "recreate");
  fx->cd();
  numL->Divide(denL);
  numL->Write();
  numK->Divide(denK);
  numK->Write();
  fx->Close();

  return 0;
};
