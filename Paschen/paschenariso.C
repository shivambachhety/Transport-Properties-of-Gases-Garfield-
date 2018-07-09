#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TGraph.h>

#include "MediumMagboltz.hh"

using namespace Garfield;
using namespace std;

void frame(TCanvas* c1, double xmin, double ymin, double xmax, double ymax) {
   c1->Range(340.8701,-0.1167765,1233.468,4.061288);
   c1->SetFillColor(0);
   c1->SetBorderMode(0); c1->SetBorderSize(2);
   c1->SetGridx(); c1->SetGridy();
   c1->SetTickx(); c1->SetTicky();
   c1->SetRightMargin(0.018); c1->SetLeftMargin(0.098);
   c1->SetTopMargin(0.015); c1->SetBottomMargin(0.09);
   c1->SetFrameLineWidth(2); c1->SetFrameBorderMode(0); c1->SetFrameBorderSize(12);
   // draw a frame to define the range
   TH1 *hframe = new TH1F("hframe","",10, xmin, xmax);
   hframe->SetMinimum(ymin); hframe->SetMaximum(ymax);
   hframe->SetDirectory(0);
   hframe->SetStats(0);
   hframe->GetXaxis()->SetTitle("d.p [cm.bar]");
   hframe->GetXaxis()->SetLabelFont(22);
   hframe->GetXaxis()->SetLabelOffset(0.007);
   hframe->GetXaxis()->SetLabelSize(0.038);
   hframe->GetXaxis()->SetTitleSize(0.044);
   hframe->GetXaxis()->SetTickLength(0.035);
   hframe->GetXaxis()->SetTitleOffset(1.04);
   hframe->GetXaxis()->SetTitleFont(22);
   hframe->GetYaxis()->SetTitle("V_{B} [V]");
   hframe->GetYaxis()->SetLabelFont(22);
   hframe->GetYaxis()->SetLabelOffset(0.003);
   hframe->GetYaxis()->SetLabelSize(0.038);
   hframe->GetYaxis()->SetTitleSize(0.044);
   hframe->GetYaxis()->SetTickLength(0.035);
   hframe->GetYaxis()->SetTitleOffset(1.12);
   hframe->GetYaxis()->SetTitleFont(22);
   hframe->Draw(" ");
}

int main(int argc, char *argv[]) {
  TApplication app("app", &argc, argv);
  // Open a canvas with a frame
  TCanvas* c1 = new TCanvas("c1", "Paschen curves for Ar/iC_{4}H_{10}",21,28,500,527);
  c1->SetLogx();
  c1->SetLogy();
  frame(c1, 3e-4, 1.0e2, 1.0e0, 1.0e4);

  // Create a legend box
  TLegend* leg = new TLegend(0.11, 0.6, 0.7, 0.97);
  leg->SetHeader("Paschen curves in Ar/iC_{4}H_{10}");

  // Paschen curve Townsend parametrisation for Ar
  for (int ig = -3; ig < 0; ig++) {
    double gamma = pow(10.0,ig);
    double A_ar = 8.63*1000.0;                // 1/Pa.m
    double B_ar = 132.0*1000.0;               // V/Pa.m
    double pdmin = 1.01*log(1+1/gamma)/A_ar;
    double pdmax = 1e0;

    double pd1[200], vb1[200];
    for (int i=0; i<200; i++) {
      pd1[i] = pdmin*pow(pdmax/pdmin, i/199.0);
      vb1[i] = B_ar*pd1[i]/(log(A_ar/log(1+1/gamma))+log(pd1[i]));
    }
    TGraph* g1 = new TGraph(200, pd1, vb1);
    g1->SetLineColor(kAzure + (ig+1));
    g1->SetLineWidth(1);
    g1->Draw();
    char buf[100];
    sprintf(buf, "Townsend, Pure Ar, #gamma = 10^{%d}", ig);
    leg->AddEntry(g1, buf, "l");
  }

  // Curves from M.J. Druyvesteyn and F.M. Penning, Rev Mod Phys 12 (1940) 87-176.
  // The Mechanism of Electrical Discharges in Gases of Low Pressure.
  double pdFe[] = {0.71755, 0.9843,  1.42324, 1.95233, 2.93121, 4.88989, 9.77222, 19.5293, 29.3212, 48.9141, 97.7525, 195.354, 237.575};
  double vbFe[] = {277.775, 261.486, 259.518, 267.479, 286.297, 330.484, 427.263, 600.244, 747.254, 1010.87, 1602.59, 2618.59, 3000};
  double pdNi[] = {0.97692, 1.42324, 1.96708, 2.93121, 4.92683, 9.77222, 19.5293, 29.3212, 48.9141, 97.7525, 195.354, 237.575};
  double vbNi[] = {194.76,  202.257, 216.486, 244.299, 299.572, 414.545, 595.726, 747.254, 1010.87, 1602.59, 2618.59, 3000};
  double pdPt[] = {0.77948, 0.9843,  1.40198, 1.95233, 2.93121, 4.88989, 9.77222, 19.3829, 22.873};
  double vbPt[] = {196.237, 191.84,  194.76,  206.893, 229.972, 279.882, 384.384, 552.382, 614};
  double pdBa[] = {0.97692, 1.36041, 1.95233, 2.93121, 4.92683, 9.77222, 19.5293, 29.3212, 48.9141, 98.4909, 195.354, 210.623};
  double vbBa[] = {115.646, 127.579, 146.161, 173.896, 226.524, 330.484, 493.206, 642.472, 909.426, 1485.98, 2465.03, 2579.33};
  double pdNa[] = {0.50378, 0.70683, 0.99173, 1.42324, 1.98194, 2.93121, 4.92683, 9.77222, 19.5293, 29.3212, 48.9141, 98.4909, 195.354, 210.623};
  double vbNa[] = {90.1294, 91.5014, 98.6813, 113.912, 133.495, 162.466, 216.486, 325.528, 493.206, 642.472, 909.426, 1485.98, 2465.03, 2579.33};

  for (int i=0; i<sizeof(pdFe)/8; i++) {pdFe[i] /= 760.0;}
  for (int i=0; i<sizeof(pdNi)/8; i++) {pdNi[i] /= 760.0;}
  for (int i=0; i<sizeof(pdPt)/8; i++) {pdPt[i] /= 760.0;}
  for (int i=0; i<sizeof(pdBa)/8; i++) {pdBa[i] /= 760.0;}
  for (int i=0; i<sizeof(pdNa)/8; i++) {pdNa[i] /= 760.0;}

  TGraph* gFe = new TGraph(sizeof(pdFe)/8, pdFe, vbFe);
  gFe->SetLineColor(kGreen+0);
  gFe->SetLineWidth(3);
  leg->AddEntry(gFe,"Druyvesteyn + Penning, Ar on Fe","l");
  gFe->Draw();

  TGraph* gNi = new TGraph(sizeof(pdNi)/8, pdNi, vbNi);
  gNi->SetLineColor(kGreen+1);
  gNi->SetLineWidth(3);
  leg->AddEntry(gNi,"Druyvesteyn + Penning, Ar on Ni","l");
  gNi->Draw();

  TGraph* gPt = new TGraph(sizeof(pdPt)/8, pdPt, vbPt);
  gPt->SetLineColor(kGreen+2);
  gPt->SetLineWidth(3);
  leg->AddEntry(gPt,"Druyvesteyn + Penning, Ar on Pt","l");
  gPt->Draw();

  TGraph* gBa = new TGraph(sizeof(pdBa)/8, pdBa, vbBa);
  gBa->SetLineColor(kGreen+3);
  gBa->SetLineWidth(3);
  leg->AddEntry(gBa,"Druyvesteyn + Penning, Ar on Ba","l");
  gBa->Draw();

  TGraph* gNa = new TGraph(sizeof(pdNa)/8, pdNa, vbNa);
  gNa->SetLineColor(kGreen+4);
  gNa->SetLineWidth(3);
  leg->AddEntry(gNa,"Druyvesteyn, Ar on Na","l");
  gNa->Draw();

  // Read the Magboltz transport tables
  MediumMagboltz* gas = new MediumMagboltz();
  double gamma = 1e-2, p = 1;

  for (int igas=0; igas<3; igas++) {
    if (igas==0) {
      gas->LoadGasFile("ar-99-ic4h10-1-p1.gas");
    } else if (igas==1) {
      gas->LoadGasFile("ar-95-ic4h10-5-p1.gas");
    } else if (igas==2) {
      gas->LoadGasFile("ar-90-ic4h10-10-p1.gas");
    }

    double dvec[100];
    double vbvec[100];

    std::vector<double> efields,  bfields, angles;
    gas->GetFieldGrid(efields, bfields, angles);

    for (int i=0; i < efields.size(); i++) {
      double logalpha;
      gas->GetElectronTownsend(i, 0, 0, logalpha);
      double a = exp(logalpha);
      double e = efields[i];
      double d = log(1.0 + 1.0/gamma)/a;
      double vb = e*d;
      dvec[i] = d*p;
      vbvec[i] = vb;
    }

    TGraph* g = new TGraph(efields.size(), dvec, vbvec);
    g->SetLineColor(90+4*igas);
    g->SetLineWidth(2);
    g->Draw();
    if (igas==0) {
      char buf[100]; sprintf(buf, "Magboltz, Ar-iC_{4}H_{10} 99-1, #gamma = 10^{%g}",log(gamma)/log(10.0));
      leg->AddEntry(g, buf, "l");
    } else if (igas==1) {
      char buf[100]; sprintf(buf, "Magboltz, Ar-iC_{4}H_{10} 95-5, #gamma = 10^{%g}",log(gamma)/log(10.0));
      leg->AddEntry(g, buf,"l");
    } else if (igas==2) {
      char buf[100]; sprintf(buf, "Magboltz, Ar-iC_{4}H_{10} 90-10, #gamma = 10^{%g}",log(gamma)/log(10.0));
      leg->AddEntry(g, buf,"l");
    }
  }

  leg->Draw();
  c1->Print("paschenariso.ps");

  app.Run(kTRUE);
  return 0;
}

