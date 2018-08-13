//To compile : g++ simulate.cc -o simulate `root-config --cflags --glibs` DarkPhotons.cc -lgsl -lgslcblas -lm

#include "DarkPhotons.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
   srand48(std::atof(argv[3]));
   std::string scalefile;
   double map, ebeam;
   if(argc==5)
   {
      ebeam = std::atof(argv[2]);
      map = std::atof(argv[1]); 
      scalefile = argv[4];
   }
   std::string mname = argv[1];
   std::string ename = argv[2];
   std::string jnum = argv[3];
   double e_m = 0.00054;
   std::string fname = "el_sim_geant_map_" + mname + "_ebeam_" + ename + "_" + jnum + ".root";
   TFile *f = new TFile(fname.c_str(),"recreate");
   
   TTree * tree = new TTree("Events", "Tree containing Lorentz Vectors from Geant");
   TLorentzVector * evec = new TLorentzVector();
   TLorentzVector * avec = new TLorentzVector();

   tree->Branch("IncidentParticle","TLorentzVector",evec);
   tree->Branch("APrime","TLorentzVector",avec);

   DarkPhotons* dphoton = new DarkPhotons(map, 0, 1, 28, 14, 2.32, 1, scalefile);
   dphoton->PrepareTable();
   TLorentzVector* pchange = new TLorentzVector();
   double e_mom, a_z, a_t, e_z, e_t, e_x, e_y, a_E, a_y, a_x;
/*
   for(int i=0; i<100; i++)
   {
      printf("Cross section at %d Gev is %e\n", i, dphoton->TotalCrossSectionCalc(i));
   }
*/
   for(int i=0;i<1000000;i++)
   {
      TLorentzVector* pchange = dphoton->SimulateEmission(ebeam);
      evec->SetPxPyPzE(pchange->X(),pchange->Y(),pchange->Z(),pchange->E());
      a_z = sqrt(ebeam*ebeam - e_m*e_m) - evec->Z();
      a_x = evec->X();
      a_y = evec->Y();
      a_E = sqrt(a_x*a_x+a_y*a_y+a_z*a_z+map*map);
      avec->SetPxPyPzE(a_x,a_y,a_z,a_E);
      tree->Fill();
   }

   f->Write();
   f->Close();

   return 0;
}
