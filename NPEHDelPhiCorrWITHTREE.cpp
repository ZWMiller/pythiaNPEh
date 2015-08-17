//==============================================================================
//  pmainHF2e.cpp
//    
//  This is an example program to study c -> e and b -> e decays
//  in 200 GeV pp collisions with Pythia8.
//
//  The decays are stored in a ROOT tree and written to file.
//
//  Once written most things can be controlled through the runcard,
//  so there's no need to recompile.
//
//  Usage: pmainHF2e  runcard  rootfile
// 
//  Author: Thomas Ullrich  
//  Last update: September 9, 2008
//==============================================================================
#include <cmath>
#include "Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include <vector>
#include "TH2.h"
#define PR(x) std::cout << #x << " = " << (x) << std::endl;
using namespace Pythia8; 

//
//  Forward declarations
//
bool isInAcceptanceE(int, const Event&);  // acceptance filter electron candidate
bool isInAcceptanceH(int, const Event&);  // acceptance filter hadron candidate
int myEvent(Pythia&, double);            // event handler (analyze event)
double deltaPhi(double, double); 
double deltaEta(double, double);
//TH2F* deltaPhiPt = new TH2F("deltaPhiPt","",200,-10,10,200,0,20);
//TH2F* deltaEtaPt = new TH2F("deltaEtaPt","",200,-5,5,200,0,20);
vector<float> dPhiV;
vector<float> dEtaV;

//
// This structure contains all the info we 
// collect. This info is later stored in a tree.
// This is our own business and has nothing to do
// with Pythia directly.
//
struct hf2eDecay_t {
  int orig_id;         // grandmother
  int orig_status;        
    
  int hf_id;         // mother (c/b hadron)
  int hf_status;
  float hf_pt;
  float hf_pz;
  float hf_phi;
  float hf_eta;
  float hf_y;
    
  int e_id;            // electron
  int e_status;
  float e_pt;
  float e_pz;
  float e_phi;
  float e_eta;
  float e_y;
        
  int   q1_id;
  float q1_x;
  int   q2_id;
  float q2_x;
  float Q2fac;
  float alphas;
  float ptHat;
  int   nFinal;
  float pdf1;
  float pdf2;
  int   code;
  float sigmaGen;
  float weight;   // useful for normalization/x-section
  float delPhi;
  float delEta;
 
};

hf2eDecay_t hf2eDecay;

int main(int argc, char* argv[]) {
    
  if (argc != 3) {
    cout << "Usage: " << argv[0] << " runcard  rootfile" << endl;
    return 2;
  }
  char* runcard  = argv[1];
  char* rootfile = argv[2];
  const char* xmlDB    = "/star/u/zbtang/myTools/pythia8142/xmldoc";
    
  //--------------------------------------------------------------
  //  Initialization
  //--------------------------------------------------------------
    
  //
  //  ROOT
  //
  TFile *hfile  = new TFile(rootfile,"RECREATE");
  TTree tree("tree","c -> e decays pp at 200 GeV");
  tree.Branch("hf2eDecay",&hf2eDecay.orig_id, 
	      "orig_id/I:orig_status/I:"
	      "hf_id/I:hf_status/I:hf_pt/F:hf_pz/F:hf_phi/F:hf_eta/F:hf_y/F:"
	      "e_id/I:e_status/I:e_pt/F:e_pz/F:e_phi/F:e_eta/F:e_y/F:"
	      "q1_id/I:q1_x/F:q2_id/I:q2_x/F:"
	      "Q2fac/F:alphas/F:ptHat/F:nFinal/I:pdf1/F:pdf2/F:code/I:sigmaGen/F:weight/F:"
	      "dPhiV/F:dEtaV/F");
    
  //
  //  Create instance of Pythia 
  //
  Pythia pythia(xmlDB); // the init parameters are read from xml files
  // stored in the xmldoc directory. This includes
  // particle data and decay definitions.
    
  //
  // Shorthand for (static) settings
  //
  Settings& settings = pythia.settings;
    
  //
  //  Read in runcard
  //
  pythia.readFile(runcard);  
  cout << "Runcard '" << runcard << "' loaded." << endl;
    
  //
  //  Retrieve number of events and other parameters from the runcard.
  //  We need to deal with those settings ourself. Getting
  //  them through the runcard just avoids recompiling.
  //
  int  maxNumberOfEvents = settings.mode("Main:numberOfEvents");
  int  nList     = settings.mode("Main:numberToList");
  int  nShow     = settings.mode("Main:timesToShow");
  int  maxErrors = settings.mode("Main:timesAllowErrors");
  bool showCS    = settings.flag("Main:showChangedSettings");
  bool showAS    = settings.flag("Main:showAllSettings");
  int  pace = maxNumberOfEvents/nShow;
 
  //
  //  Remark: in this example we do NOT alter the
  //  BRs since they are different for the various charm
  //  hadrons making the normalization at the end rather
  //  cumbersome. In a production version this is what
  //  one probably would implement to save processing time.
  //
    
  //
  //  Initialize Pythia, ready to go
  //
  pythia.init();
    
  //
  // List changed or all data
  //
  if (showCS) settings.listChanged();
  if (showAS) settings.listAll();
    
  //--------------------------------------------------------------
  //  Event loop
  //--------------------------------------------------------------
  int ievent = 0;
  int numberOfElectrons = 0;
  int iErrors = 0;
  int n;
    
  while (ievent < maxNumberOfEvents) {
        
    if (!pythia.next()) {
      if (++iErrors < maxErrors) continue;
      cout << "Error: too many errors in event generation - check your settings & code" << endl;
      break;
    }
    n = myEvent(pythia, maxNumberOfEvents);  // in myEvent we deal with the whole event and return
    // the number of electrons recorded for book keeping
    numberOfElectrons += n; 
    if (n) tree.Fill();   
    ievent++;
    if (ievent%pace == 0) {
      cout << "# of events generated = " << ievent 
	   << ", # of electrons from c/b hadron decays generated so far = " << numberOfElectrons << endl;
    }
        
    // List first few events.
    if (ievent < nList) {
      pythia.info.list(); 
      pythia.process.list();
      pythia.event.list();
    }
  }
    
  //--------------------------------------------------------------
  //  Finish up
  //--------------------------------------------------------------
  pythia.statistics();
  cout << "Writing File" << endl;
  hfile->Write();
    
  return 0;
}

//
//  Event analysis
//
int myEvent(Pythia& pythia, double nMaxEvt)
{
  Event &event = pythia.event;

  int nelectrons = 0;
  int ic = 0;
  for (int i = 0; i < event.size(); i++) {
    if (abs(event[i].id()) == 11) { // event is electron

      //
      //  Check if mother is a c/b hadron
      //
      vector<int> mothers = event.motherList(i);
      if (mothers.size() != 1) {
	cout << "Error: electron has more than one mother. Stop." << endl;
	abort();
      }
      ic = mothers[0];
      int ic_id = abs(event[ic].id());
      int flavor = static_cast<int>(ic_id/pow(10.,static_cast<int>(log10(ic_id))));
      if (flavor != 4 && flavor != 5) continue; // c (b) hadrons start with 4(5)  

      //
      //  Acceptance filter
      //    
      if (!(isInAcceptanceE(i, event))) continue;
            
      nelectrons++;
            
      //
      // Get grandmother (origin of c/b hadron)
      //
      vector<int> grandmothers = event.motherList(ic);
      int iorig = -1;
      switch(grandmothers.size()) {
      case 0: 
	iorig = -1;
	break;
      case 1:
	iorig = grandmothers[0];
	break;
      default:
	iorig = -2;
	break;
      }
            
      // At this point we have J/psi and mother B, the J/psi detectable in
      // STAR.                                                                               // Now Collect (i) all hadrons and (ii) all hadrons from the B.   
      // We require them to be stable, i.e. not decayed.             
      // Also impose pt cut on hadrons as in data.                                                                                  
      vector<int> hadrons;
      vector<int> B_hadrons;
            
      for (int i = 1; i < event.size(); i++) {
        if (event[i].isFinal() && event[i].isCharged() && event[i].pT() > 0.2 && isInAcceptanceH(i, event)) {
	  hadrons.push_back(i);
	  //	  if (event.isAncestor(i, i_B)) B_hadrons.push_back(i); // From Bingchu code, save in case needed later
        }
      }

      //
      //  Store in tuple
      //
            
      // If no origin or more than 1 than id == 0
      // and the status identifies what happens: -1 no mother, 
      // -2 more than 1. This should not happen at all, but ...
      hf2eDecay.orig_id     = iorig >= 0 ? event[iorig].id() : 0;
      hf2eDecay.orig_status = iorig >= 0 ? event[iorig].status() : iorig;
            
      hf2eDecay.hf_id     = event[ic].id();  
      hf2eDecay.hf_status = event[ic].status();
      hf2eDecay.hf_pt     = event[ic].pT();
      hf2eDecay.hf_pz     = event[ic].pz();
      hf2eDecay.hf_phi    = event[ic].phi();
      hf2eDecay.hf_eta    = event[ic].eta();     
      hf2eDecay.hf_y      = event[ic].y();
            
      hf2eDecay.e_id       = event[i].id();     
      hf2eDecay.e_status   = event[i].status();
      hf2eDecay.e_pt    = event[i].pT();
      hf2eDecay.e_pz    = event[i].pz();
      hf2eDecay.e_phi    = event[i].phi();
      hf2eDecay.e_eta      = event[i].eta();     
      hf2eDecay.e_y    = event[i].y();
            
      hf2eDecay.q1_id      = pythia.info.id1();
      hf2eDecay.q1_x       = pythia.info.x1();
      hf2eDecay.q2_id      = pythia.info.id2();
      hf2eDecay.q2_x       = pythia.info.x2();
      hf2eDecay.Q2fac      = pythia.info.Q2Fac();
      hf2eDecay.alphas     = pythia.info.alphaS();
      hf2eDecay.ptHat      = pythia.info.pTHat();
      hf2eDecay.nFinal     = pythia.info.nFinal();
      hf2eDecay.pdf1       = pythia.info.pdf1();
      hf2eDecay.pdf2       = pythia.info.pdf2();
      hf2eDecay.code       = pythia.info.code();
      hf2eDecay.sigmaGen   = pythia.info.sigmaGen();
      hf2eDecay.weight     = pythia.info.sigmaGen()/nMaxEvt; // useful for obtaining x-section
      // vector branches already assigned in branch declaration

      double phi1, phi2;
      double eta1, eta2;
      int hid;
      phi1 = event[i].phi();
      eta1 = event[i].eta();
      for (unsigned int ih=0; ih<hadrons.size(); ih++) {
        hid = hadrons[ih];
        phi2 = event[hid].phi();
	eta2 = event[hid].eta();
        float dphi = deltaPhi(phi1, phi2);
        float deta = deltaEta(eta1, eta2);
        if(event[hid].pT()<0.2) continue;
	dPhiV.push_back(dphi);
	dEtaV.push_back(deta);
	//deltaPhiPt -> Fill(dphi,(float)event[i].pT());
	//deltaEtaPt -> Fill(deta,(float)event[i].pT());
      }                          
    }    

    return nelectrons;
  }
}

//
//  Acceptance filter
//
bool isInAcceptanceE(int i, const Event& event)
{
  // accept all (useful for many studies)
  //  return true;
    
  // limit to STAR TPC/BEMC/ToF acceptance
  double eta = event[i].eta();
  if (fabs(eta) < 0.7)
      return true;
  else
      return false;
}

bool isInAcceptanceH(int i, const Event& event)
{
  // limit to STAR TPC/BEMC/ToF acceptance                                              
  double eta = event[i].eta();
  if (fabs(eta) < 1)
    return true;
  else
    return false;
}

double deltaPhi(double phi1, double phi2)
{
  // move to range [0, 2pi]                                                           
  if (phi1<0) phi1 += 2*M_PI;
  if (phi2<0) phi2 += 2*M_PI;

  // correct difference                                                               
  double delta = phi2-phi1;
  if (delta<-0.5*M_PI) delta += 2*M_PI;
  if (delta> 1.5*M_PI) delta -= 2*M_PI;

  return delta;
}

double deltaEta(double e1, double e2)
{
  double delta = e2-e1;
  return delta;
}
