#include <iostream>

#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>

#include "ViewSignal.hh"
#include "ComponentAnalyticField.hh"
#include "MediumMagboltz.hh"
#include "SolidTube.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "Sensor.hh"
#include "AvalancheMC.hh"
#include "AvalancheMicroscopic.hh"
#include "TrackHeed.hh"
#include "Random.hh"
#include "Plotting.hh"

using namespace Garfield;

// Finds ion creation time based on distance from the electron cluster 
// to the wire (found by fitting the ion creation time of many avalanches)
double IonTiming(double dist) {

  double p0 =  1.49880e-13;   
  double p1 =  2.09250e+02;
  double p2 =  2.61998e+02;
  double p3 = -1.24766e+02;
  
  return p0 + p1 * dist + p2 * dist * dist + p3 * dist * dist * dist;

}

int main() {

  plottingEngine.SetDefaultStyle();
 
  // Variables used to describe the geometry
  const double dWire = 50.e-4;
  const double rAnode = 1.46;

  // Variables used to describe the field
  const double vWire = 3270.;
  const double vAnode = 0.;

  // Variables describing the signal histogram
  const double tmin = 0.;
  const double tstep = 1.;
  const int nTimeBins = 1000; 

  // Average ion creation point
  const double xIon = 0.00253;

  // Location of track
  double x0 = 1.2;
  double y0 = -(sqrt(rAnode * rAnode - x0 * x0) - 0.002);
  double z0 = 0.;
  double t0 = 0.;

  // Cluster 
  double xc, yc, zc, tc, ec, extra;
  int nc;

  // Variables used to describe the avalanche
  double gain = 1190;
  double theta = 0.2654;
  
  // Make a gas medium
  MediumMagboltz* gas = new MediumMagboltz();
  const std::string path = getenv("GARFIELD_HOME");
  gas->LoadIonMobility(path + "/Data/IonMobility_Ar+_Ar.txt");
 
  // Build the geometry.
  GeometrySimple* geo = new GeometrySimple();
  SolidBox* box = new SolidBox(0., 0., 0., rAnode, rAnode, rAnode);
  geo->AddSolid(box, gas);

  // Make a component with analytic electric field.
  ComponentAnalyticField* comp = new ComponentAnalyticField();
  comp->SetGeometry(geo);
  comp->AddWire(0., 0., dWire, vWire, "s", 100., 50., 19.3);
  comp->AddTube(rAnode, vAnode, 0, "t");
  comp->AddReadout("s");

  // Make a sensor.
  Sensor* sensor = new Sensor();
  sensor->AddComponent(comp);
  sensor->AddElectrode(comp, "s");
  sensor->SetArea(-rAnode, -rAnode, -rAnode, rAnode, rAnode, rAnode);
  sensor->SetTimeWindow(tmin, tstep, nTimeBins);
  sensor->ClearSignal();
 
  // MC integration
  AvalancheMC* drift = new AvalancheMC();
  drift->SetSensor(sensor);
  drift->SetDistanceSteps(2.e-4);
  drift->EnableSignalCalculation();
 
  // Setup HEED
  TrackHeed* track = new TrackHeed();
  track->SetParticle("muon");
  track->SetEnergy(170.e9);
  track->SetSensor(sensor);

  // Start simulation 
  sensor->ClearSignal();
  sensor->NewSignal();
  track->NewTrack(x0, y0, z0, t0, 0., 1., 0.);
  double time;
  double location;
  int count = 0;
  // Loop over all clusters created by muon
  while (track->GetCluster(xc, yc, zc, tc, nc, ec, extra)) {
    count++;
    // Find radial location of cluster
    location = sqrt(xc * xc + yc * yc);
    // Find the creation time of the ions created by the electron cluster
    time = IonTiming(location);
    std::cout << "Cluster: " << count << std::endl;
    std::cout << "Size of Cluster: " << nc << std::endl;
    std::cout << "Location of Cluster: " << location << std::endl;

    // Calculate the signal induced only by the electrons inside the drift tube
    if (location <= rAnode ){
      std::cout << "Cluster is inside the tube" << std::endl;
      // For each electron in the cluster:
      for (int i = 0; i < nc; i++) {
        // Get the number of electron tracks.
        double np = gain * RndmPolya(theta);
        // Scale the effect of the ion induced signal by the number of electron tracks.
        drift->SetIonSignalScalingFactor(np);
        // Drift 1 highly charged ion to estimate the signal created by the avalanche from 1 electron
        drift->DriftIon(xIon, 0., 0., time);
      }
    } else {
      std::cout << "Cluster outside of tube" << std::endl;
    }
  }
  std::cout << "Plotting Signal..." << std::endl;
  ViewSignal* signalView = new ViewSignal();
  signalView->SetSensor(sensor);
  TCanvas* c1 = new TCanvas("c1", "Signal", 21, 28, 500, 527);
  signalView->SetCanvas(c1);
  signalView->PlotSignal("s");
  // Save the plotted driftlines to a pdf
  c1->SaveAs("ApproxSignal.pdf");
  
}

  

