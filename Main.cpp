/*
 * CorrFitSimple.cxx
 *
 *  Created on: 19 cze 2021
 *      Author: Daniel Wielanek
 *		E-mail: daniel.wielanek@gmail.com
 *		Warsaw University of Technology, Faculty of Physics
 */
#include "Histogram1D.h"
#include "StandAloneFsiLednicky.h"
#include "StandAloneSimpleFsi.h"
#include "TRandom.h"
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <iostream>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

using namespace std;


FourVector GetMom(double mass) {
  FourVector l;
  l.x      = gRandom->Gaus(0, 0.1);
  l.y      = gRandom->Gaus(0, 0.1);
  l.z      = gRandom->Gaus(0, 0.1);
  double p = sqrt(l.x * l.x + l.y * l.y + l.z * l.z + mass * mass);
  l.t      = p;
  return l;
}

void GetPos(double* p) {
  const Double_t R = 4.0;
  p[0]             = gRandom->Gaus(0, R);
  p[1]             = gRandom->Gaus(0, R);
  p[2]             = gRandom->Gaus(0, R);
  p[3]             = 0;  // gRandom->Gaus(0, 1);
}


void Boost(Pair* p) {
  const double gFmToGeVToFm = 5.067731036;
  const double gFmToGeV     = 0.197326968;

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

TH1D* GetRootHist(TString name, const Histogram1D& h) {
  TH1D* res = new TH1D(name, name, h.GetNBins(), h.GetMin(), h.GetMax());
  for (int i = 0; i < h.GetNBins(); i++) {
    res->SetBinContent(i + 1, h.GetBinContent(i));
  }
  return res;
}


int main(int argc, char** argv) {
  /** for random **/
  std::random_device rd {};
  std::mt19937 gen {rd()};
  std::normal_distribution<double> d {0, 0.5};
  /**end to random **/


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

  tWeightKisiel->setPairType(theChosenOne);

  Pair* pair = new Pair();

  TH1D* numL         = new TH1D("numL", "numL", 100, 0, 1);
  TH1D* denL         = new TH1D("denL", "denL", 100, 0, 1);
  TH1D* numK         = new TH1D("numK", "numK", 100, 0, 1);
  TH1D* denK         = new TH1D("denK", "denK", 100, 0, 1);
  TH1D* dif          = new TH1D("dif", "dif", 100, -0.1, 0.1);
  Histogram1D* numLH = new Histogram1D(0, 1);
  Histogram1D* denLH = new Histogram1D(0, 1);

  double poz1[4], poz2[4];
  Int_t pairs  = 1000000;
  Int_t kPairs = pairs / 1000;
  for (int i = 0; i < pairs; i++) {
    if (i % 1000 == 0) {
      cout << i / 1000 << " of " << kPairs << endl;
      ;
    }
    FourVector p1      = GetMom(0.139);
    FourVector p2      = GetMom(0.493);
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
    double wL = 0;  // tWeightLed->getWeight(*pair);


    double kstar = tWeightKisiel->getKStarSigned();


    // cout << "REQ" << endl;
    // cout << poz1[0] << " " << poz1[1] << " " << poz1[2] << " " << poz1[3] << endl;
    // cout << tWeightKisiel->GetRout() << " " << tWeightKisiel->GetRside() << " " << tWeightKisiel->GetRlong() << endl;
    // if (fabs((wK - wL) * 100.) > 0.01) cout << wK << " " << wL << "----------- " << kstar << endl;

    dif->Fill(100.0 * (wK - wL) / wK);

    numL->Fill(kstar, wL);
    numK->Fill(kstar, wK);
    numLH->Fill(kstar, wL);
    denLH->Fill(kstar, 1);
    denL->Fill(kstar);
    denK->Fill(kstar);
  }

  TFile* fx = new TFile(systemName, "recreate");
  fx->cd();
  numL->Divide(denL);
  numL->Write();
  numK->Divide(denK);
  numK->Write();

  TH1D* nc = new TH1D("checkNum", "checkNum", 100, 0, 1);
  TH1D* dc = new TH1D("checkDen", "checkDen", 100, 0, 1);
  for (int i = 0; i < 100; i++) {
    nc->SetBinContent(i + 1, numLH->GetBinContent(i));
    dc->SetBinContent(i + 1, denLH->GetBinContent(i));
  }
  nc->Write();
  dc->Write();

  dif->Write();
  fx->Close();

  return 0;
};
