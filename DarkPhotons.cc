// This is a class to simulate dark photons A production by electrons in matter eN -> eNA
// To be used in a Geant4 application.
//
// Reference: arXiV:1604.08432 [hep-ph]
//
// The example code to use it is below
//
// In user RunAction constructor:
//   myDarkPhotons = new DarkPhotons(MA, EThresh, SigmaNorm, AA, ZZ, Density, epsilon);
//             where MA              - dark photon mass, GeV
//                   EThresh         - Threshold energy for electron to emit dark photon
//                   SigmaNorm       - Additional cross section factor to provide rare event condition
//                   AA, ZZ, Density - Atomic number, nucleus charge and density of the element (default is lead)
//                   epsilon         - mixing, is used only in the calculation of the final normalization
//
// In user RunAction::BeginOfRunAction:
//   myDarkPhotons->myDarkPhotons->PrepareTable();
//
// In user RunAction::EndOfRunAction:
//   cout << myDarkPhotons->GetNormalization();
//
// In user SteppingAction::UserSteppingAction:
//
//   if( theParticleDefinition == G4Electron::ElectronDefinition() ||
//       theParticleDefinition == G4Positron::PositronDefinition() ) {
//
//     G4double DensityMat = aStep->GetTrack()->GetMaterial()->GetDensity()/(g/cm3);
//
//     if( GetDarkPhotonsPointer()->Emission(ekin, DensityMat, StepLength) ) {
//
//       if(NEmissions) {
//         G4cout << "Attempt to emit more than one A in the event, rejected!!! E = " << ekin << G4endl;
//       }
//       else {
//
//        NEmissions++;
//        //cout << " Ee = " << ekin << " Steplength = " << StepLength << endl;
//        G4double XAcc = GetDarkPhotonsPointer()->SimulateEmission(ekin);
//        track->SetKineticEnergy(ekin*(1.-XAcc)*GeV); // electron loses XAcc of its energy on A emission
//
//
//
#include "DarkPhotons.hh"
//#include "Randomize.hh"

#include <gsl/gsl_math.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>

#include <iostream>
#include <fstream>
#include <sstream>

#define Mel 5.1E-4 // electron mass (GeV)
#define Mmu 0.1056 // muon mass (GeV)
#define alphaEW 1.0/137.0
#define MUp 2.79 // protonMu
#define Mpr 0.938 // proton mass (GeV)
#define max_uint 4294967296.0l
#define GeVtoPb 3.894E+08

#define NPTAB 15


/* +++++++++++++++++++ Util routines: Spline interpolation ++++++++++++ */
/* ++++++++++++++++++++++ C Include Files ++++++++++++++++++++++ */
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <bits/stdc++.h>
#include "TLorentzVector.h"

#define EPSPARINV 1.e-8

double parinv(double x, double a[], double f[], int n)
{
//
//    Interpolation at the point x. Function f(a) is tabulated
//    in arrays a, f with dimension n.
//
  int k1, k2, k3;

  if(n < 3) {std::cerr << "parinv: insufficient number of points" << std::endl; exit(1);}
  if(x < a[0]) {
    double c = fabs(x - a[0]);
    if(c < EPSPARINV*fabs(a[1]-a[0])) return a[0];
    k1 = 0;
  }
  else if(x > a[n-1]) {
    double c = fabs(x - a[n-1]);
    if(c < EPSPARINV*fabs(a[n-1]-a[n-2])) return a[n-1];
    k1 = n-3;
  }
  else {
    k1 = 0;
    k2 = n-1;
    k3 = k2 - k1;
    while(k3 > 1) {
      k3 = k1 + k3/2;
      if( a[k3]-x == 0 ) return f[k3];
      if( a[k3]-x < 0 ) k1 = k3;
      if( a[k3]-x > 0 ) k2 = k3;
      k3 = k2 - k1;
    }
    if(k2 == n-1) k1 = n - 3;
  }
  if(k1 < 0 || k1 > n-3) {std::cerr << "parinv: wrong index found" << std::endl; exit(1);}
  double b1 = a[k1];
  double b2 = a[k1+1];
  double b3 = a[k1+2];
  double b4 = f[k1];
  double b5 = f[k1+1];
  double b6 = f[k1+2];
  return b4 * ((x-b2)*(x-b3))/((b1-b2)*(b1-b3)) +
         b5 * ((x-b1)*(x-b3))/((b2-b1)*(b2-b3)) +
         b6 * ((x-b1)*(x-b2))/((b3-b1)*(b3-b2));
}


DarkPhotons::DarkPhotons(double MAIn, double EThreshIn, double SigmaNormIn, double ANuclIn, double ZNuclIn, double DensityIn,
                         double epsilIn, std::string fname)
:MA(MAIn), EThresh(EThreshIn), SigmaNorm(SigmaNormIn),
ANucl(ANuclIn), ZNucl(ZNuclIn), Density(DensityIn), epsilBench(0.0001), epsil(epsilIn),
AccumulatedProbability(0.)
{
   nptable = NPTAB;
   double epi[NPTAB]={MA+Mmu, MA+Mmu+.1,MA+Mmu+2., MA+Mmu+10., MA+Mmu+20., MA+Mmu+50., MA+Mmu+100., MA+Mmu+150., MA+Mmu+250., MA+Mmu+500., MA+Mmu+800., MA+Mmu+1500., MA+Mmu+2000., MA+Mmu+2500., MA+Mmu+4000.};
   for(int ip=0; ip < nptable; ip++) {ep[ip] = epi[ip];}
   ParseLHE(fname, 5.0); 
   std::sort(energies.begin(), energies.end());
}


DarkPhotons::~DarkPhotons()
{
}

double DsigmaDx(double x, void * pp) 
{
   ParamsForChi* params = (ParamsForChi*)pp;
   
   double beta = sqrt(1 - (params->MMA)*(params->MMA)/(params->EE0)/(params->EE0));
   double num = 1.-x+x*x/3.;
   double denom = (params->MMA)*(params->MMA)*(1.-x)/x+Mel*Mel*x;
   double DsDx = beta*num/denom;

   return DsDx;
}

void DarkPhotons::ParseLHE(std::string fname, double ebeam)
{
   std::ifstream ifile;
   ifile.open(fname.c_str());
   if (!ifile)
   {
      std::cout << "Unable to open LHE file\n";
      exit(1);
   }   
   std::string line;
   int nentries=0;
   mgdata[ebeam];
   printf("Making energy %e.\n",ebeam);
   while(std::getline(ifile, line))
   {
      std::istringstream iss(line);
      int ptype, state;
      double skip, px, py, pz, E, pt, efrac, M;
      if (iss >> ptype >> state >> skip >> skip >> skip >> skip >> px >> py >> pz >> E >> M ) 
      {
         if((ptype==11)&&(state==1))
	 {
	    pt = sqrt(px*px+py*py);
            efrac = (E - M)/(ebeam-M-MA);
	    nentries++;
            mgdata[ebeam].push_back(std::make_pair(efrac,pt));
	 }
      }
   }
   energies.push_back(std::make_pair(ebeam,int(drand48()*nentries)));
   ifile.close();
}

double chi (double t, void * pp) 
{
  ParamsForChi* params = (ParamsForChi*)pp;

/* Reminder II:
   params->AA;
   params->ZZ;
   params->MMA;
   params->EE0;
*/

  double d = 0.164/pow((params->AA),2./3.);
  double ap = 773.0/Mel/pow((params->ZZ),2./3.);
  double a = 111.0/Mel/pow((params->ZZ),1./3.);
  double G2el = (params->ZZ)*(params->ZZ)*a*a*a*a*t*t/(1.0+a*a*t)/(1.0+a*a*t)/(1.0+t/d)/(1.0+t/d);
  double G2in = (params->ZZ)*ap*ap*ap*ap*t*t/(1.0+ap*ap*t)/(1.0+ap*ap*t)/(1.0+t/0.71)/(1.0+t/0.71)
    /(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)
    *(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr)*(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr);
  double G2 = G2el+G2in;
  double ttmin = (params->MMA)*(params->MMA)*(params->MMA)*(params->MMA)/4.0/(params->EE0)/(params->EE0);
  //double ttmin = lowerLimit(x,theta,p);
  double Under = G2*(t-ttmin)/t/t;
//  std::cout << "Under: " << Under << " d: " << d << " AA " << params->AA << " a: " << a << std::endl;
  return Under;
}


double DarkPhotons::TotalCrossSectionCalc(double E0)
{
  double Xmax;
  double sigmaTot;

    if(E0 < 2.*MA) return 0.;

    //begin: chi-formfactor calculation

    gsl_integration_workspace * w
      = gsl_integration_workspace_alloc (1000);

    double result, error;
    double tmin = MA*MA*MA*MA/(4.*E0*E0);
    double tmax = MA*MA;

    gsl_function F;
    ParamsForChi alpha = {1.0, 1.0, 1.0, 1.0};
    F.function = &chi;
    F.params = &alpha;

    alpha.AA = ANucl;
    alpha.ZZ = ZNucl;
    alpha.MMA = MA;
    alpha.EE0 = E0;

    gsl_integration_qags (&F, tmin, tmax, 0, 1e-7, 1000,
                          w, &result, &error);

    //printf ("chi/Z^2 = % .18f\n", result/(ZNucl*ZNucl));
    //printf ("result    = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

    double ChiRes = result;
//    std::cout << "Chi: " << result << " E0: " << E0 << " MA: " << MA << std::endl;

    gsl_integration_workspace_free (w);

    gsl_integration_workspace * s 
       = gsl_integration_workspace_alloc (1000);
    gsl_function G;
    G.function = &DsigmaDx;
    G.params = &alpha;
    double xmin = 0;
    double xmax = 1;
    if((Mel/E0)>(MA/E0)) xmax = 1-Mel/E0;
    else xmax = 1-MA/E0;
    double res, err;

    gsl_integration_qags (&G, xmin, xmax, 0, 1e-7, 1000,
                          s, &res, &err);
    double DsDx = res;

    gsl_integration_workspace_free(s);
//    end: chi-formfactor calculation

    sigmaTot= GeVtoPb*4.*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*DsDx;
    if(sigmaTot < 0.) sigmaTot=0.;

    return sigmaTot;
}

double DsigmaDxmu(double x, void * pp)
{
   ParamsForChi* params = (ParamsForChi*)pp;

   double MMu = 105.658/1000.;
   double beta = sqrt(1- (params->MMA)*(params->MMA)/(params->EE0)/(params->EE0));
   double num = 1.-x+x*x/3.;
   double denom = (params->MMA)*(params->MMA)*(1.-x)/x+MMu*MMu*x;
   double DsDx = beta*num/denom;

   return DsDx;
}

double chimu (double t, void * pp) {
  ParamsForChi* params = (ParamsForChi*)pp;

/* Reminder II:
   params->AA;
   params->ZZ;
   params->MMA;
   params->EE0;
*/

  double d = 0.164/pow((params->AA),2./3.);
  double ap = 773.0/Mel/pow((params->ZZ),2./3.);
  double a = 111.0/Mel/pow((params->ZZ),1./3.);
  double G2el = (params->ZZ)*(params->ZZ)*a*a*a*a*t*t/(1.0+a*a*t)/(1.0+a*a*t)/(1.0+t/d)/(1.0+t/d);
  double G2in = (params->ZZ)*ap*ap*ap*ap*t*t/(1.0+ap*ap*t)/(1.0+ap*ap*t)/(1.0+t/0.71)/(1.0+t/0.71)
    /(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)/(1.0+t/0.71)
    *(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr)*(1.0+t*(MUp*MUp-1.0)/4.0/Mpr/Mpr);
  double G2 = G2el+G2in;
  double ttmin = (params->MMA)*(params->MMA)*(params->MMA)*(params->MMA)/4.0/(params->EE0)/(params->EE0);
  //double ttmin = lowerLimit(x,theta,p);
  double Under = G2*(t-ttmin)/t/t;
//  std::cout << "Under: " << Under << " MMA: " << params->MMA << " EE0: " << params->EE0 << std::endl;
  return Under;
}

double DarkPhotons::TotalMuCrossSectionCalc(double E0)
{
  double Xmin;
  double Xmax;
  double sigmaTot;

    if(E0 < 2.*MA) return 0.;

    Xmin = MA/E0;
    Xmax = 1.0-Xmin;

    //begin: chi-formfactor calculation

    gsl_integration_workspace * w
      = gsl_integration_workspace_alloc (1000);

    double result, error;
    double tmin = MA*MA*MA*MA/(4.*E0*E0);
    double tmax = MA*MA;

    gsl_function F;
    ParamsForChi alpha = {1.0, 1.0, 1.0, 1.0};
    F.function = &chimu;
    F.params = &alpha;

    alpha.AA = ANucl;
    alpha.ZZ = ZNucl;
    alpha.MMA = MA;
    alpha.EE0 = E0;

    gsl_integration_qags (&F, tmin, tmax, 0, 1e-7, 1000,
                          w, &result, &error);

    //printf ("chi/Z^2 = % .18f\n", result/(ZNucl*ZNucl));
    //printf ("result    = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", w->size);

   double ChiRes = result;
   gsl_integration_workspace_free (w);
//    end: chi-formfactor calculation

   gsl_integration_workspace * dxspace = gsl_integration_workspace_alloc (1000);
   gsl_function G;
   G.function = &DsigmaDxmu;
   G.params = &alpha;
   double xmin = 0;
   double xmax = 1;
   if((Mmu/E0)>(MA/E0)) xmax = 1-Mmu/E0;
   else xmax = 1-MA/E0;
   double res, err;

   gsl_integration_qags (&G, xmin, xmax, 0, 1e-7, 1000, dxspace, &res, &err);
   
   double DsDx = res;
   gsl_integration_workspace_free(dxspace);

   sigmaTot = GeVtoPb*4.*alphaEW*alphaEW*alphaEW*epsilBench*epsilBench*ChiRes*DsDx;
   if(sigmaTot < 0.)
   {
      sigmaTot = 0.;
   }

   return sigmaTot;
}
double DarkPhotons::MaxCrossSectionCalc(double E0)
{
  double Xmax, Xmin;
  double sigmaMax;

    if(E0 < 2.*MA) return 0.;

    Xmin = MA/E0;
    Xmax = 0.998;

    double UxthetaMax = MA*MA*(1. - Xmax)/Xmax + Mel*Mel*Xmax;
    double AAMax = (1. - Xmax + Xmax*Xmax/2.0) / (UxthetaMax*UxthetaMax);
    double BBMax = (1. - Xmax)*(1. - Xmax)*MA*MA / (UxthetaMax*UxthetaMax*UxthetaMax*UxthetaMax);
    double CCMax = MA*MA - UxthetaMax*Xmax/(1. - Xmax);
    sigmaMax = Xmax * (AAMax + BBMax*CCMax);

    return sigmaMax;
}


void DarkPhotons::PrepareTable()
{
  for(int ip=0; ip < nptable; ip++) {
    sigmap[ip] = TotalMuCrossSectionCalc(ep[ip]);
    sigmax[ip] = MaxCrossSectionCalc(ep[ip]);
  }
}


double DarkPhotons::GetsigmaTot(double E0)
{
  if(E0<(MA+Mmu))
  {
     return 0;
  }
  double st = parinv(E0, ep, sigmap, nptable);
  if(st<0) {return TotalMuCrossSectionCalc(E0);}
  return parinv(E0, ep, sigmap, nptable);
}


double DarkPhotons::GetsigmaMax(double E0)
{
  return parinv(E0, ep, sigmax,	nptable);
}


bool DarkPhotons::Emission(double E0, double DensityMat, double StepLength)
{
  if(E0 < EThresh) return false;
  if(fabs(DensityMat - Density) > 0.1) return false;
  double prob = SigmaNorm*GetsigmaTot(E0)*StepLength;
  AccumulatedProbability += prob;
  double tmprandom = 0.4;//G4UniformRand();
  if(tmprandom < prob) return true;
  return false;
}

TLorentzVector* DarkPhotons::SimulateEmission(double E0)
{
   TLorentzVector* fParticle = new TLorentzVector;
   double Eout, Theta, Phi, Eta;
   std::pair < double, double > data = GetMadgraphData(E0);
   double XAcc = data.first;
   double Pt = data.second;
   int i = 0;
   while (Pt + Mel > XAcc)
   {
      i++;
      data = GetMadgraphData(E0);
      XAcc = data.first;
      Pt = data.second;
//      printf("Skipped an event\n.");
      if(i>10000)
      {
         printf("Did not manage to simulate.\n");
	 return fParticle;
      }
   }
   double P = sqrt(XAcc*XAcc-Mel*Mel);
   Eout = XAcc;
   Theta = asin(Pt/P);
   Eta = -log(tan(Theta/2.));
   Phi = drand48()*2.*3.14159;
   fParticle->SetPtEtaPhiE(Pt,Eta,Phi,Eout);
   return fParticle;
}

std::pair <double, double> DarkPhotons::GetMadgraphData(double E0)
{
   double samplingE = energies[0].first;
   double efrac = 0;
   double pt = 0;
   int i=0;
   bool pass = false;
   while(!pass)
   {
      i = i+1;
      samplingE = energies[i].first;
      if(E0<=samplingE) {pass=true;}
      if(i>=energies.size()) {pass=true;}
   }
   if(i>0) {i=i-1;}
   
   if(energies[i].second>=mgdata[energies[i].first].size())
   {
      energies[i].second = 0;
   }

   efrac = mgdata[energies[i].first].at(energies[i].second).first;
   pt = mgdata[energies[i].first].at(energies[i].second).second;
   energies[i].second=energies[i].second+1;
   
   double avail = E0 - MA - Mel;
//   printf("Selected energy %e.\n", energies[i].first);
   return std::make_pair(efrac*avail+Mel,pt);
}

TLorentzVector* DarkPhotons::MuSimulateEmission(double E0)
{
   TLorentzVector* fParticle = new TLorentzVector;
   double Eout, Theta, Phi, Eta;
   double ptmax = sqrt(E0*E0-MA*MA);
//   double XAcc = GetMadgraphE(E0);
   double XAcc = 0;
   double width = 1/sqrt((XAcc-Mmu)/(E0-MA-Mmu))/0.8/sqrt(MA)+1.4/MA;
   double integratedPx = 1./width-exp(-width*XAcc)/width;
   for(int i=0;i<10000;i++)
   {
      double PxA = drand48()*integratedPx;
      double PyA = drand48()*integratedPx;
      double Px = -log(1-PxA*width)/width;
      double Py = -log(1-PyA*width)/width;
      double Pt = sqrt(Px*Px+Py*Py);
      if (Pt*Pt+Mmu*Mmu < XAcc*XAcc) 
      {
         double P = sqrt(XAcc*XAcc-Mmu*Mmu);
         Eout = XAcc;
         Theta = asin(Pt/P);
         Eta = -log(tan(Theta/2));
         Phi = drand48()*2*3.14159;
         fParticle->SetPtEtaPhiE(Pt,Eta,Phi,Eout);
         return fParticle;
      }
   }   
   printf ("Did not manage to simulate!. Xacc: %e, E0: %e\n", XAcc, E0);

   return fParticle; // did not manage to simulate
}
