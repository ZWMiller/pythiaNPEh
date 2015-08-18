//==============================================================================
//  pmainHF2e.cpp
//    
//  This is an example program to study c -> e and b -> e decays
//  in 200 GeV pp collisions with Pythia8.
//
//  Creates templated for deltaPhi in b or c, depending on unput cards.
//
//  Usage: pmainHF2e  runcard  rootfile
// 
//  Author: Thomas Ullrich  
//  Last update: September 9, 2008
//  Modified: Z.W. Miller Aug 17, 2015
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
bool isInAcceptanceE(int, const Event&);  // acceptance filter electron candidate
bool isInAcceptanceH(int, const Event&);  // acceptance filter hadron candidate
int myEvent(Pythia&, vector<TH2D*> &, vector<TH3D*>&, double); // event handler (analyze event)
double deltaPhi(double, double); 
double deltaEta(double, double);

int main(int argc, char* argv[]) {
    
  if (argc != 4) {
    cout << "Usage: " << argv[0] << " runcard rootfile histName" << endl;
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
  cout << "============================================================================" \
       << endl;
  cout << "Executing program '" << argv[0] << "', start at: " << ctime(&now);
  cout << "Arguments: " << argv[1] << " " << argv[2] << " " << argv[3] << endl;
  cout << "============================================================================" \
       << endl;
    
  //
  //  ROOT
  //
  char text[64];
  TFile *hfile  = new TFile(rootfile,"RECREATE");
  vector<TH2D*> histos2D;
  sprintf(text,"histos2D%s%d",histname,0);
  histos2D.push_back(new TH2D(text,"NPE - h", 150,0.,15.,200, -0.5*M_PI, 1.5*M_PI));
  sprintf(text,"histos2D%s%d",histname,1);
  histos2D.push_back(new TH2D(text,"NPE pt vs y", 150, 0., 15., 60, -3, 3));
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
  histos3D.push_back(new TH3D(text,"NPE - h", 150,0.,15.,150,0,15,200, -0.5*M_PI, 1.5*M_PI));
  sprintf(text,"histo3D%s%d",histname,1);
  histos3D.push_back(new TH3D(text,"NPE - B-->h", 150,0.,15.,150,0,15,200, -0.5*M_PI, 1.5*M_PI));


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
    n = myEvent(pythia, histos2D, histos3D, maxNumberOfEvents);  // in myEvent we deal with the whole event and return
    // the number of electrons recorded for book keeping
    if(n == 0) continue;
    numberOfElectrons += n; 
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

  now = time(0);
  cout << "============================================================================\
" << endl;
  cout << "Program finished at: " << ctime(&now);
  cout << "============================================================================\
" << endl;
    
  return 0;
}

//
//  Event analysis
//
int myEvent(Pythia& pythia, vector<TH2D*> &histos2D, vector<TH3D*> &histos3D, double nMaxEvt)
{
  Event &event = pythia.event;

  int nelectrons = 0;
  int ic = 0;
  int ie = 0;
  for (int i = 0; i < event.size(); i++) {
    if (abs(event[i].id()) == 11) { // event is electron

      //
      //  Check if mother is a c/b hadron
      //
      vector<int> mothers = event.motherList(i);
      if (mothers.size() > 1) {
	cout << "Error: electron has more than one mother. Stop." << endl;
	//abort();
	return 0;
      }
      ic = mothers[0];
      ie = i;
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
    }
  }          
  // At this point we have J/psi and mother B, the J/psi detectable in
  // STAR.
  // We require them to be stable, i.e. not decayed.
  // Also impose pt cut on hadrons as in data.
  vector<int> hadrons;
  vector<int> B_hadrons;
  
  for (int i = 1; i < event.size(); i++) {
    if (event[i].isFinal() && event[i].isCharged() && event[i].pT() > 0.2 && isInAcceptanceH(i, event) && !(event[i].id()==event[ie].id())) {
      hadrons.push_back(i);
      //	  if (event.isAncestor(i, i_B)) B_hadrons.push_back(i); // From Bingchu code, save in case needed later
    }
  }
  
  //
  //  Fill histograms                                                       
  //

  //histos[2]->Fill(event[i_B].pT(), 1.);                                       
  histos2D[1]->Fill(event[ie].pT(), event[ie].y());
  Double_t npept = event[ie].pT();
  double phi1, phi2;
  int nnear = 0;
  int naway = 0;
  double ptbalance = npept;
  int hid;
  double dphi=999;
  phi1 = event[ie].phi();
  
  for (unsigned int i=0; i<B_hadrons.size(); i++) {
    hid = B_hadrons[i];
    phi2 = event[hid].phi();
    histos2D[9]->Fill(npept, event[hid].pT());
    histos3D[1]->Fill(npept, event[hid].pT(), deltaPhi(phi1, phi2));
  }
  
  for (unsigned int i=0; i<hadrons.size(); i++) {
    hid = hadrons[i];
    phi2 = event[hid].phi();
    if(!(phi1==0) && !(phi2==0))
      dphi = deltaPhi(phi1, phi2);
    histos3D[0]->Fill(npept, event[hid].pT(), dphi);
    if(event[hid].pT()<0.5) continue;
    histos2D[0]->Fill(npept, dphi);
    if( abs(dphi) < 1) {//near side                                                   
      nnear++;
      ptbalance += event[hid].pT();
      histos2D[4]->Fill(npept, event[hid].pT());
      histos2D[6]->Fill(npept, event[hid].m0());
    }
    if (abs(dphi-M_PI)<1) { //away side                                               
      naway++;
      histos2D[5]->Fill(npept, event[hid].pT());
      histos2D[7]->Fill(npept, event[hid].m0());
      ptbalance -= event[hid].pT();
    }
  }
  histos2D[2]->Fill(npept, nnear);
  histos2D[3]->Fill(npept, naway);
  histos2D[8]->Fill(npept, ptbalance);
  hadrons.clear();
  
  return nelectrons;
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
