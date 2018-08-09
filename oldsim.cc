//To compile : g++ simulate.cc -o simulate `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm

//#include "DarkPhotons.hh"
#include "Dp_NewSim.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
   double map, ebeam;
   if(argc==3)
   {
      ebeam = std::atof(argv[2]);
      map = std::atof(argv[1]); 
   }
   std::string mname = argv[1];
   std::string ename = argv[2];
   double e_m = 0.00054;
   std::string fname = "oldsim_geant_map_" + mname + "_ebeam_" + ename + ".root";
   TFile *f = new TFile(fname.c_str(),"recreate");
   
   TTree * tree = new TTree("Events", "Tree containing Lorentz Vectors from Geant");
   TLorentzVector * evec = new TLorentzVector();
   TLorentzVector * avec = new TLorentzVector();

   tree->Branch("IncidentParticle","TLorentzVector",evec);
   tree->Branch("APrime","TLorentzVector",avec);

   DarkPhotons* dphoton = new DarkPhotons(map, 0, 1, 28, 14, 2.32, 1);
   momentum pchange;
   double a_mom, a_z, a_t, e_z, e_t, e_x, e_y, e_E, a_y, a_x, a_E;

   for(int i=0;i<1000000;i++)
   {
      pchange = dphoton->SimulateEmission(ebeam);
      a_E = pchange.E0*ebeam;
      a_mom = sqrt(a_E*a_E-map*map);
      a_z = a_mom*cos(pchange.Theta);
      a_t = a_mom*sin(pchange.Theta);
      a_x = a_t*cos(pchange.Phi);
      a_y = a_t*sin(pchange.Phi);
      e_z = sqrt(ebeam*ebeam - e_m*e_m) - a_z;
      e_t = a_t;
      e_x = e_t*cos(pchange.Phi);
      e_y = e_t*sin(pchange.Phi);
      e_E = sqrt(e_x*e_x+e_y*e_y+e_z*e_z+e_m*e_m);
      avec->SetPxPyPzE(a_x,a_y,a_z,a_E);
      evec->SetPxPyPzE(e_x,e_y,e_z,e_E); 
      tree->Fill();
   }

   f->Write();
   f->Close();

   return 0;
}
