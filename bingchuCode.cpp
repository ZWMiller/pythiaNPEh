//==============================================================================
//  pmainBJpsiHcorr.cpp
//    
//  Study B -> J/Psi + X / J/Psi hadron correlation 
//  in 200 GeV pp collisions with Pythia8.
//
//  Usage: pmainBJpsiHcorr  runcard  rootfile
// 
//  Author: Thomas Ullrich  
//  Last update: October28, 2008
//  Modified by Zebo Tang
//==============================================================================
#include <ctime>
#include <cmath>
#include <vector>
#include "Pythia.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#define PR(x) std::cout << #x << " = " << (x) << std::endl;
using namespace Pythia8; 

//
//  Forward declarations
//
bool isInAcceptance(int, const Event&);  // acceptance filter
int myEvent(Pythia&, vector<TH2D*> &, vector<TH3D*> &);     // event handler 
int parentB(int, const Event&);
double deltaPhi(double, double);

int main(int argc, char* argv[]) {

  if (argc != 4) {
    cout << "Usage: " << argv[0] << " runcard  rootfile histname" << endl;
    return 2;
  }
  char* runcard  = argv[1];
  char* rootfile = argv[2];
  char* histname = argv[3];
  const char* xmlDB    = "/star/u/zbtang/myTools/pythia8142/xmldoc";

  //--------------------------------------------------------------
  //  Initialization
  //--------------------------------------------------------------

  time_t now = time(0);
  cout << "============================================================================" << endl;
  cout << "Executing program '" << argv[0] << "', start at: " << ctime(&now);
  cout << "Arguments: " << argv[1] << " " << argv[2] << " " << argv[3] << endl;
  cout << "============================================================================" << endl;

  //
  //  ROOT
  //
  char text[64];
  TFile *hfile  = new TFile(rootfile,"RECREATE");
  vector<TH2D*> histos2D;
  sprintf(text,"histos2D%s%d",histname,0);
  histos2D.push_back(new TH2D(text,"J/psi - h", 150,0.,15.,40, -0.5*M_PI, 1.5*M_PI));
  sprintf(text,"histos2D%s%d",histname,1);
  histos2D.push_back(new TH2D(text,"J/psi pt vs y", 150, 0., 15., 60, -3, 3));
  sprintf(text,"histos2D%s%d",histname,2);
  histos2D.push_back(new TH2D(text,"near-side Nch", 150, 0, 15, 50, 0., 50.));
  sprintf(text,"histos2D%s%d",histname,3);
  histos2D.push_back(new TH2D(text,"away-side Nch", 150, 0, 15, 50, 0., 50.));
  sprintf(text,"histos2D%s%d",histname,4);
  histos2D.push_back(new TH2D(text,"near-side pt", 150, 0, 15, 150, 0., 15.));
  sprintf(text,"histos2D%s%d",histname,5);
  histos2D.push_back(new TH2D(text,"away-side pt", 150, 0, 15, 150, 0., 15.));
  sprintf(text,"histos2D%s%d",histname,6);
  histos2D.push_back(new TH2D(text,"near-side m0", 150, 0, 15, 100, 0., 1.));
  sprintf(text,"histos2D%s%d",histname,7);
  histos2D.push_back(new TH2D(text,"away-side m0", 150, 0, 15, 100, 0., 1.));
  sprintf(text,"histos2D%s%d",histname,8);
  histos2D.push_back(new TH2D(text,"pt balance", 150, 0, 15, 100, -10, 10.));
  sprintf(text,"histos2D%s%d",histname,9);
  histos2D.push_back(new TH2D(text,"B daughter pt", 150, 0, 15, 150, 0, 15.));

  vector<TH3D*> histos3D;
  sprintf(text,"histo3D%s%d",histname,0);
  histos3D.push_back(new TH3D(text,"J/psi - h", 150,0.,15.,150,0,15,40, -0.5*M_PI, 1.5*M_PI));
  sprintf(text,"histo3D%s%d",histname,1);
  histos3D.push_back(new TH3D(text,"J/psi - B-->h", 150,0.,15.,150,0,15,40, -0.5*M_PI, 1.5*M_PI));
    //for (unsigned int k=0; k<histos.size(); k++) histos[k]->Sumw2();
    
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
    //  Force J/Psi to decay into ee so we can check if
    //  it would fall into STAR's acceptance.
    //
    pythia.readString("443:onMode = off");
    pythia.readString("443:onIfMatch = 11 -11");

    //
    //  Suppress low pT divergence
    //  (not needed if ptHat cut is on)
    //
    UserHooks *oniumUserHook = new SuppressSmallPT();
    pythia.setUserHooksPtr(oniumUserHook);
    
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
    int numberOfJpsi = 0;
    int iErrors = 0;
    int n;
    
    while (ievent < maxNumberOfEvents) {
        
        if (!pythia.next()) {
            if (++iErrors < maxErrors) continue;
            cout << "Error: too many errors in event generation - check your settings & code" << endl;
            break;
        }
        n = myEvent(pythia, histos2D, histos3D);
        
        numberOfJpsi += n; 
        ievent++;
        if (ievent%pace == 0) {
            cout << "# of events generated = " << ievent 
            << ", # of  numberOfJpsi generated so far = " << numberOfJpsi << endl;
        }
        
        // List first few events.
        if (ievent < nList && n>0) {
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

    now = time(0);
    cout << "============================================================================" << endl;
    cout << "Program finished at: " << ctime(&now);
    cout << "============================================================================" << endl;
    
    return 0;
}

//
//  Event analysis
//
int myEvent(Pythia& pythia, vector<TH2D*> &histos2D, vector<TH3D*> &histos3D)
{    
    Event &event = pythia.event;
    
    //
    // First get the J/psi and the mother B.
    // We skip the case where we have two J/Psi
    //
    int i_B = -1;
    int i_Jpsi = -1;
    int nJpsi = 0;
    for (int i = 1; i < event.size(); i++) {
        if (abs(event[i].id()) == 443 && event[i].pT() > 0) {
            i_B = parentB(i, event);
            if (i_B == -1) {
                cout << "Warning: found J/psi but no B mother" << endl;
                continue;
            }
			//cout<<"jpsi mother = "<<event[i_B].id()<<endl;
            i_Jpsi = i;
            nJpsi++;   // count only if there's a B parent
        }
    }
    
    // 
    // Make sense what we got?
    //
    if (nJpsi == 0) return 0;
    if (nJpsi != 1) {
        cout << "Warning: more than one J/psi (with B parent) in event: n = " << nJpsi << endl;
        return 0;
    }
    if (i_B == -1) {
        cout << "Warning: J/psi but no B parent in event" << endl;
        return 0;
    }
    
    //
    //  Make sure both electrons from J/psi make it in the STAR acceptance
    //
    vector<int> daughters = event.daughterList(i_Jpsi);
    if (daughters.size() != 2) {
        cout << "Error: J/Psi doesn't have 2 daughters n = " << daughters.size() << endl;
        return 0;
    }
    int id1 = daughters[0];
    int id2 = daughters[1];
    if (! (abs(event[id1].id()) == 11 && abs(event[id2].id()) == 11)) {
        cout << "Error: J/Psi didn't decay into e+e-" << endl;
        return 0;
    }
    if (!(isInAcceptance(id1, event) && isInAcceptance(id2, event))) return 0;
    
    //
    // At this point we have J/psi and mother B, the J/psi detectable in
    // STAR.
    //
    // Now Collect (i) all hadrons and (ii) all hadrons from the B.
    // We require them to be stable, i.e. not decayed.
    // Also impose pt cut on hadrons as in data.
    //
    vector<int> hadrons;
    vector<int> B_hadrons;
    
    for (int i = 1; i < event.size(); i++) {
        if (event[i].isFinal() && event[i].isCharged() && event[i].pT() > 0. && isInAcceptance(i, event)) {
            hadrons.push_back(i);
            if (event.isAncestor(i, i_B)) B_hadrons.push_back(i);
        }
    }

    //
    //  Fill histograms
    //

    //histos[2]->Fill(event[i_B].pT(), 1.);
    histos2D[1]->Fill(event[i_Jpsi].pT(), event[i_Jpsi].y());
    Double_t jpsiPt = event[i_Jpsi].pT();
    double phi1, phi2;
    int nnear = 0;
    int naway = 0;
    double ptbalance = jpsiPt;
    int hid;
    phi1 = event[i_Jpsi].phi();
    
    for (unsigned int i=0; i<B_hadrons.size(); i++) {
        hid = B_hadrons[i];
        phi2 = event[hid].phi();
        histos2D[9]->Fill(jpsiPt, event[hid].pT());
	histos3D[1]->Fill(jpsiPt, event[hid].pT(), deltaPhi(phi1, phi2));
    }
   
    for (unsigned int i=0; i<hadrons.size(); i++) {
        hid = hadrons[i];
        phi2 = event[hid].phi();
	double dphi = deltaPhi(phi1, phi2);
	histos3D[0]->Fill(jpsiPt, event[hid].pT(), dphi);
	if(event[hid].pT()<0.5) continue;
	histos2D[0]->Fill(jpsiPt, dphi);
	if( abs(dphi) < 1) {//near side
	  nnear++;
	  ptbalance += event[hid].pT();
	  histos2D[4]->Fill(jpsiPt, event[hid].pT());
	  histos2D[6]->Fill(jpsiPt, event[hid].m0());
	}
	if (abs(dphi-M_PI)<1) { //away side
	  naway++;
	  histos2D[5]->Fill(jpsiPt, event[hid].pT());
	  histos2D[7]->Fill(jpsiPt, event[hid].m0());
	  ptbalance -= event[hid].pT();
	}
    }
    histos2D[2]->Fill(jpsiPt, nnear);
    histos2D[3]->Fill(jpsiPt, naway);
    histos2D[8]->Fill(jpsiPt, ptbalance);
    hadrons.clear();
    
    return 1;    
}

//
// Recursive search for parent B. Returns -1 if not found. 
// Since we start the loop at i=1 we start with the earliest
// B meson in the decay chain. This might be a higher B
// resonance that decays into a B0, B+- etc.
//
int parentB(int ijpsi, const Event& event) 
{    
    for (int k = 1; k < event.size(); k++) {
        if (abs(event[k].id()) > 500 && abs(event[k].id()) < 600) {  // found any B
            if (event.isAncestor(ijpsi, k))  return k;  // is mother
        }
    }
    return -1;
}

//
//  Acceptance filter
//
bool isInAcceptance(int i, const Event& event)
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
