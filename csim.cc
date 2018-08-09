//To compile : g++ simulate.cc -o simulate `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm

#include "combsim.hh"
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
   std::string fname = "combsim_geant_map_" + mname + "_ebeam_" + ename + ".root";
   TFile *f = new TFile(fname.c_str(),"recreate");
   
   TTree * tree = new TTree("Events", "Tree containing Lorentz Vectors from Geant");
   TLorentzVector * evec = new TLorentzVector();
   TLorentzVector * avec = new TLorentzVector();

   tree->Branch("IncidentParticle","TLorentzVector",evec);
   tree->Branch("APrime","TLorentzVector",avec);

   DarkPhotons* dphoton = new DarkPhotons(map, 0, 1, 28, 14, 2.32, 1);
   momentum pchange;
   double e_mom, a_z, a_t, e_z, e_t, e_x, e_y, e_E, a_y, a_x, a_E;

   for(int i=0;i<1000000;i++)
   {
      pchange = dphoton->SimulateEmission(ebeam);
      e_E = pchange.E0;
      e_mom = sqrt(e_E*e_E-e_m*e_m);
      e_z = e_mom*cos(pchange.Theta);
      e_t = e_mom*sin(pchange.Theta);
      e_x = e_t*cos(pchange.Phi);
      e_y = e_t*sin(pchange.Phi);
      a_z = sqrt(ebeam*ebeam - map*map) - e_z;
      a_t = e_t;
      a_x = a_t*cos(pchange.Phi);
      a_y = a_t*sin(pchange.Phi);
      a_E = sqrt(a_x*a_x+a_y*a_y+a_z*a_z+map*map);
      avec->SetPxPyPzE(a_x,a_y,a_z,a_E);
      evec->SetPxPyPzE(e_x,e_y,e_z,e_E); 
      tree->Fill();
   }

   f->Write();
   f->Close();

   return 0;
}
