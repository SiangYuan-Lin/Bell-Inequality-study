/*

  Examines H->WW leptonic decay events for violation of the Bell inequality

*/

#include "hww_results.h"

/// ROOT includes
#include "TH1.h"
#include "TProfile.h"
#include "TH2.h"
#include "TFile.h"
#include "TChain.h"
#include "TMatrixD.h"
#include "TVector.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"

// Tree Reader
#include "ExRootAnalysis/ExRootTreeReader.h"

// DELPHES includes
#include "classes/DelphesClasses.h"

// c++ standard library includes
#include <iomanip>
#include <iostream>
#include <sstream>

//------------------------------------------------------------------------------

struct Plots{
public:
  TH1* hist_lepton_n;
  TH1* hist_lepton_pt;
  TH1* hist_lepton_dphi;
  TH1* hist_met;
  TH1* hist_recoil;
  TH1* hist_dilepton_mass;
  TH1* hist_m_cm;
  TH1* hist_nW;
  TH1* hist_mW;
  TH1* hist_p_cm_W;
  TH1* hist_pt_cm_W;
  TH1* hist_chargeW;
  TH1* hist_mlplusnu;
  TH1* hist_mlminusnubar;
  TH1* hist_nZ;
  TProfile* hist_Cab_elements;
  TProfile* hist_cglmp_elements;
  std::vector<TH1*> hist_cos_theta;
  std::vector<TH1*> hist_cossq_theta;
  std::vector<TH1*> hist_epsilon;
  TH1* hist_sum_cossq_theta;
  TH1* hist_spacelike_invariant;
};

void bookhists(Plots* plots){
  plots->hist_lepton_n = new TH1F("h_lept_mult","Lepton multiplicity",21,-0.5,20.5);
  plots->hist_lepton_pt = new TH1F("h_lept_pt","Lepton pt",50,0,250);
  plots->hist_lepton_dphi = new TH1F("h_lept_dphi","Lepton delta phi",50,0,TMath::Pi());
  plots->hist_met = new TH1F("h_met","Met",50,0,250);
  plots->hist_recoil = new TH1F("h_recoil","Recoil",50,0,250);
  plots->hist_dilepton_mass = new TH1F("h_dilepton_mass","dileptonic invariant mass",50,0,250);
  plots->hist_nW = new TH1F("h_W_mult","W boson multiplicity",6,-0.5,5.5);
  plots->hist_m_cm = new TH1F("h_m_cm","m(llnunu)",150,0,300);;
  plots->hist_mW = new TH1F("h_W_mass","W boson mass",100,0,100);
  plots->hist_chargeW = new TH1F("h_W_charge","W boson charge",3,-1.5,1.5);
  plots->hist_p_cm_W = new TH1F("h_W_p_cm","W boson pt*",100,0,100);
  plots->hist_pt_cm_W = new TH1F("h_W_pt_cm","W boson p*",100,0,100);

  plots->hist_mlplusnu = new TH1F("h_lnu_mass","l+nu mass",100,0,100);
  plots->hist_mlminusnubar = new TH1F("h_lnubar_mass","l-nubar mass",100,0,100);
  plots->hist_nZ = new TH1F("h_Z_mult","Z boson multiplicity",6,-0.5,5.5);
  plots->hist_Cab_elements = new TProfile("h_Cab_elements","elements of Cab",9,0.5,9.5);
  plots->hist_sum_cossq_theta = new TH1F("h_sum_cossq_theta","Sum of cos^2 theta",200,0.,2.5);

  plots->hist_cglmp_elements = new TProfile("h_cglmp_elements","elements of cglmp",9,0.5,9.5);
  
  plots->hist_spacelike_invariant  = new TH1F("h_spacelike_invariant","gamma1gamma2(1+v1v2)",200,0.,10.);
  
  // plots of the cos theta distributions for each lepton and for each axis
  for (int ilept=1; ilept<=2; ilept++){
    for (int icart=0; icart<3; icart++){
      std::ostringstream name;
      name << "cos_theta_l" << ilept << "_b" << icart;
      plots->hist_cos_theta.push_back(new TH1F(name.str().c_str(), name.str().c_str(), 100., -1., 1.));
      std::ostringstream name2;
      name2 << "cossq_theta_l" << ilept << "_b" << icart;
      plots->hist_cossq_theta.push_back(new TH1F(name2.str().c_str(), name2.str().c_str(), 100., -1., 1.));
    }
  }

  // plots of the epsilon_ij distributions
  for (int j=0; j<3; j++){
    for (int i=0; i<3; i++){
      std::ostringstream name;
      name << "epsilon_" << i << "_" << j;
      plots->hist_epsilon.push_back(new TH1F(name.str().c_str(), name.str().c_str(), 100., -1., 1.));
    }
  }
  
  
}

void writehists(Plots*plots){
  plots->hist_lepton_n->Write();
  plots->hist_lepton_pt->Write();
  plots->hist_lepton_dphi->Write();
  plots->hist_met->Write();
  plots->hist_recoil->Write();
  plots->hist_dilepton_mass->Write();
  plots->hist_nW->Write();
  plots->hist_mW->Write();
  plots->hist_p_cm_W->Write();
  plots->hist_pt_cm_W->Write();
  plots->hist_m_cm->Write();
  plots->hist_chargeW->Write();
  plots->hist_mlplusnu->Write();
  plots->hist_mlminusnubar->Write();
  plots->hist_nZ->Write();
  plots->hist_Cab_elements->Write();
  plots->hist_cglmp_elements->Write();
  plots->hist_sum_cossq_theta->Write();
  plots->hist_spacelike_invariant->Write();
  for (int ih=0; ih<plots->hist_cos_theta.size(); ih++){
    plots->hist_cos_theta[ih]->Write();
    plots->hist_cossq_theta[ih]->Write();
  }
  for (int ih=0; ih<plots->hist_epsilon.size(); ih++){
    plots->hist_epsilon[ih]->Write();
  }
  
};


/* Calculate the cglmp Bell-like inequality for qutrits
 using the method of Phys. Rev. A 94, 032119 as adapted
 for Higgs decays by Barr 
 arguments: 
 p_x = cosine of angle of first lepton to x axis
 p_y = cosine of angle of first lepton to y axis
 m_x = cosine of angle of second lepton to x axis
 m_y = cosine of angle of second lepton to y axis
 returns: contribution to <B_cglmp> operator for a single event
*/
double calculate_cglmp(double p_x, double p_y, double m_x, double m_y){
  double value = (8./TMath::Sqrt(3))*(p_x*m_x+p_y*m_y)
    + 25. * (p_x*p_x - p_y*p_y)*(m_x*m_x - m_y*m_y)
    + 100. * p_x * p_y * m_x * m_y;
  return value;
};



double calculate_dphi(double phi1, double phi2){
  double pi = TMath::Pi();
  double d = abs(phi1-phi2);
  if (d>pi) d = 2*pi-d;
  return d;
};

//------------------------------------------------------------------------------

inline bool islepton(const GenParticle *gen){
  const Int_t pid = gen->PID;
  return (abs(pid)==11 || abs(pid)==13);
};

inline bool isWboson(const GenParticle *gen){
  const Int_t pid = gen->PID;
  return (abs(pid)==24);
};

inline bool isZboson(const GenParticle *gen){
  return (gen->PID==23);
};

//------------------------------------------------------------------------------

inline bool isinvisible(const GenParticle *gen){
  const Int_t pid = gen->PID;
  return (abs(pid)==12 || abs(pid)==14 || abs(pid)==16);
};

//------------------------------------------------------------------------------

inline bool isfinal(const GenParticle *gen){
  return gen->Status==1;
};

//------------------------------------------------------------------------------

// The main analysis routine
void truth_bell(const char *inputFile, const char* outputFile, hww::values& inout_values)
{
  using namespace std;
  
  //  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create output histograms
  TFile* output=new TFile(outputFile,"RECREATE");
  cout << "creating output file " << outputFile << endl;
  Plots* myplots = new Plots();
  bookhists(myplots);

  // The choice of basis is important. Here we have the option of various different bases based on the value of this flag:
  // (1) A fixed basis defined by the detector - simple but not very helpful (wrong)
  // (2) A basis defined relative to the beamline (z) and with the x-axis along the missing transverse momentum direction. This seems dodgy, since it breaks the independence of the axes and the measured objects. (i.e. this is wrong)
  // (3) The axes defined around eqn (6) of 2102.11883, but using the W boson rather than the top quark
  // (4) Same as (3) but constructing W+ boson from e+ + nu rather than looking in the event record (since sometimes there is no W in the event record)
  
  // All defined in the CM frame of the H boson  
  const int basis_choice=4;

  // decide which frame to measure the lepton directions in
  // (1) In the CM frame of the higgs boson system (this is wrong)
  // (2) In the CM frame of the direct parent (i.e. Z boson), as recommended in 2102.11883
  const int lepton_frame_choice=2;

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  // The correllation part of the CHSH density matrix, equation (9) of 2102.11883,
  // which may be obtained by integrating over all events
  double correlation_ab[3][3]={0.};
  double correlation_ab_sq[3][3]={0.};
  // we'll later divide by the number of times the quantity above was incrimented

  // the contributions to the CGLMP operator (and its square) 
  double cglmp_contrib[3][3]={0.};
  double cglmp_contrib_sq[3][3]={0.};
  
  int n_correlation=0;

  int n_basis_problems=0;

  const double cross_section = 513.; // fb (ACCORDING TO THE HXSWG 48.6 picobarns * 1.055*10^-2)
  //const double cross_section = 150.; // fb (ACCORDING TO THE HXSWG 48.6 picobarns * 1.055*10^-2)
  
  const double target_lumi   = 140.;   // inverse fb
  const double lumi_scale_factor = cross_section*target_lumi / numberOfEntries;
  std::cout << "x-sec = " << cross_section << " fb " << endl
	    << "target lumi = " << target_lumi << " /fb " << endl
	    << "expected events =" << cross_section*target_lumi <<endl
	    << "generated events = " << numberOfEntries << endl
	    << "scale factor = " << lumi_scale_factor << std::endl;
  
  const double min_W_mass = inout_values.mW_min;
  const double max_W_mass = inout_values.mW_max;
  
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    
    if(entry%50000==0) cout << "Event:\t" << entry <<endl;
    bool printout=entry<2;

    std::vector<GenParticle*> final_state_leptons;

    // MET for the event. NOTE that it's unphysical to know all the components, but let's start there
    TLorentzVector sum_of_invisibles_lab;
    std::vector<GenParticle*> invisibles;

    // Look for W bosons and Z bosons in the lab
    std::vector<GenParticle*> W_bosons;
    std::vector<GenParticle*> Z_bosons;
    
    // loop over all input particles in the event
    for(Int_t i=0; i < branchParticle->GetEntriesFast(); i++)
    {    
     GenParticle *gen = (GenParticle*) branchParticle->At(i);     
     if (islepton(gen) && isfinal(gen)) {
       final_state_leptons.push_back(gen);
       //cout<<"N: "<<i<<", St: "<<gen->Status<<", PID: "<<gen->PID<<", E: "<<gen->E<<", Px: "<<gen->Px<<", Py: "<<gen->Py<<", Pz: "<<gen->Pz<<", M: "<<gen->Mass<<", M1: "<<gen->M1<<", M2: "<<gen->M2<<", D1: "<<gen->D1<<", D2: "<<gen->D2<<endl;
       myplots->hist_lepton_pt->Fill(gen->PT);
     }
     if (isinvisible(gen) && isfinal(gen)){
       sum_of_invisibles_lab += gen->P4();  
       invisibles.push_back(gen);
     }

     // find a W boson if there is one
     if (isWboson(gen)){
       W_bosons.push_back(gen);
     }
     if (isZboson(gen)){
       Z_bosons.push_back(gen);
     }
    } // end loop over particles

    //std::cout << "Number of final state leptons:" << final_state_leptons.size() << std::endl;
    myplots->hist_lepton_n->Fill(final_state_leptons.size());
    myplots->hist_met->Fill(sum_of_invisibles_lab.Pt());

    // Do some checks about the number of W bosons
    int nW=W_bosons.size(); 	
    if (printout) {
      cout << "Number of W bosons: "<< nW << endl;
    }
    
    myplots->hist_nW->Fill(nW);
    myplots->hist_nZ->Fill(Z_bosons.size());
    for (int iW=0; iW<nW; ++iW){
      myplots->hist_mW->Fill(W_bosons[iW]->Mass);
      myplots->hist_chargeW->Fill(W_bosons[iW]->Charge);
    }
    
    // Events with two final state leptons
    if (final_state_leptons.size()==2){
      const int q1 = final_state_leptons[0]->Charge;
      const int q2 = final_state_leptons[1]->Charge;

      // only take opposite-sign pairs
      if ((q1*q2)>=0) continue;

      // only take different-flavour pairs
      if (abs(final_state_leptons[0]->PID)==abs(final_state_leptons[1]->PID)) continue;

      // get the lab lorentz vectors of the two leptons according to their charges.
      TLorentzVector lplus_lab;
      TLorentzVector lminus_lab;
      if (final_state_leptons[0]->Charge>0){
	lplus_lab = final_state_leptons[0]->P4();
	lminus_lab = final_state_leptons[1]->P4();
      } else {
	lplus_lab = final_state_leptons[1]->P4();   AJB corrected Jan 2022. Did not affect H->WW paper
	lminus_lab = final_state_leptons[0]->P4();  AJB corrected Jan 2022. Did not affect H->WW paper
      }
      
      // confirm that there were exactly two invisibles - or print warning
      if (invisibles.size()!=2) cerr << "warning - there were " << invisibles.size()
				     << "invisibles in this event" << endl;
      
      // reconstruct W+ and W- four-vectors in case they aren't in the event record
      // from four-momenta of leptons plus neutrinos
      TLorentzVector WPlusFourVectorReco, WMinusFourVectorReco;
      // first add the leptons
      if (q1>0){
	WPlusFourVectorReco   += final_state_leptons[0]->P4();
	WMinusFourVectorReco  += final_state_leptons[1]->P4();
      }else{
	WPlusFourVectorReco   += final_state_leptons[1]->P4();
	WMinusFourVectorReco  += final_state_leptons[0]->P4();
      }
      // now add the neutrinos
      for (int inu=0; inu<invisibles.size(); inu++){
	GenParticle* nu = invisibles[inu];
	// if it's an electron or muon neutrino then add it to the WPlus:
	if (nu->PID==12 || nu->PID==14) {
	  WPlusFourVectorReco += nu->P4();
	} else if (nu->PID==-12 || nu->PID==-14) {
	  // its an anti-neutrino of e or mu type so add it to the W-
	  WMinusFourVectorReco += nu->P4();
	}	
      }// end loop over neutrinos

      myplots->hist_mlplusnu->Fill(WPlusFourVectorReco.M());
      myplots->hist_mlminusnubar->Fill(WMinusFourVectorReco.M());

      // calculate the Lorentz invariant of the W decay vertices
      double dx = final_state_leptons[0]->X - final_state_leptons[1]->X;
      double dy = final_state_leptons[0]->Y - final_state_leptons[1]->Y;
      double dz = final_state_leptons[0]->Z - final_state_leptons[1]->Z;
      double dt = final_state_leptons[0]->T - final_state_leptons[1]->T;

      if (printout){
	cout << "dx,dy,dz,dt=\t" << dx << " " << dy << " " << dz << " " << dt << std::endl;
      }
      
      // fill with the Lorentz invariant -- convert from mm to fm
      //myplots->hist_spacelike_invariant->Fill((dx*dx+dy*dy+dz*dz-dt*dt)*1.E12);
      
      
      // dilepton four-momentum sum
      const TLorentzVector dilepton_lab = lplus_lab+lminus_lab;
      
      double dphi = calculate_dphi(lplus_lab.Phi(), lminus_lab.Phi());
      myplots->hist_lepton_dphi->Fill(dphi);

      double dilepton_mass=dilepton_lab.M();
      myplots->hist_dilepton_mass->Fill(dilepton_mass);

      // Lorentz vector summing two leptons and invisibles
      const TLorentzVector system_lab = dilepton_lab+sum_of_invisibles_lab;
      myplots->hist_recoil->Fill(system_lab.Pt());

      // Boost to the CM of the dilepton + invisibles - inverse Lorentz Transofrm required, hence -ve sign
      TVector3 boost_lab_to_cm = -system_lab.BoostVector();
      TLorentzVector lplus_cm = lplus_lab;   lplus_cm.Boost(boost_lab_to_cm);
      TLorentzVector lminus_cm = lminus_lab; lminus_cm.Boost(boost_lab_to_cm);
      TLorentzVector system_cm = system_lab; system_cm.Boost(boost_lab_to_cm);      
      
      myplots->hist_m_cm->Fill(system_cm.M());
      
      // Need also to find the beam direction in the CM frame
      TVector3 beam_lab_direction(0,0,1.);
      TLorentzVector beam_lab(beam_lab_direction,1.);
      TLorentzVector beam_cm = beam_lab; beam_cm.Boost(boost_lab_to_cm);
      
      // // Calculate the vectors defined in arXiv:2102.11883
      if (printout) {
	cout << "\nNEW EVENT number " << entry << endl;
       	cout << "l+ \t";    lplus_lab.Print();
	cout << "l- \t";    lminus_lab.Print();
      	cout << "MET\t";   sum_of_invisibles_lab.Print();
	cout << "lab\t";  system_lab.Print();
	cout << "Boost\t"; boost_lab_to_cm.Print();
	cout << "CM\t";    system_cm.Print() ;
	cout << endl;
      }
      
      // make the cuts here...
      if (WPlusFourVectorReco.M() < min_W_mass || WMinusFourVectorReco.M() < min_W_mass) continue;
      if (WPlusFourVectorReco.M() > max_W_mass || WMinusFourVectorReco.M() > max_W_mass) continue;

      // Define the axes in a different way to Eqn (6) of 2102.11883,
      // since we dont have access to the W bosons' directions
      TVector3 basis[3];

      TLorentzVector W_plus_lab, W_plus_cm;
      TLorentzVector W_minus_lab, W_minus_cm;
      
      // choose basis with which to evalate the correlation matrix elements
      if (basis_choice==1){
	// use the cartesian detector axes in the boosted frame
	basis[0] = TVector3(1,0,0);
	basis[1] = TVector3(0,1,0);
	basis[2] = TVector3(0,0,1); // along the beam line
	cerr << "The vector above should be boosted to the in the CM frame" << endl;
	exit(1);
	
      } else if (basis_choice==2){
	// define z along the beam
	// x along the MET direction in the transverse plane
	// y perpendicular to these two
	// all measured in the boosted frame
	TLorentzVector sum_of_invisibles_cm = sum_of_invisibles_lab;
	sum_of_invisibles_cm.Boost(boost_lab_to_cm);
	double phi_invis_cm = sum_of_invisibles_cm.Phi();
	basis[0] = TVector3(TMath::Cos(phi_invis_cm), TMath::Sin(phi_invis_cm),0);
	basis[1] = TVector3(-TMath::Sin(phi_invis_cm), TMath::Cos(phi_invis_cm),0);
	basis[2] = TVector3(0,0,1); // along the beam line
	cerr << "The vector above should be boosted to the in the CM frame" << endl;
	exit(1);
	
      } else if (basis_choice==3 || basis_choice==4){
	// Choose axes as per 2102.11883 but replacing top quark with W boson
	// either from the event record (basis choice 3)
	// or reconstructed W boson (basis choice 4)
	
	if (basis_choice==3 && printout){
	  for (int iW=0; iW<nW; ++iW){
	    cout << "#: " << iW << " ID:" << W_bosons[iW]->PID << " Mass: " <<  W_bosons[iW]->Mass << endl;
	    W_bosons[iW]->P4().Print();
	  }
	}
	// If we're using the event record W and there are none, we're stuck:
	if (basis_choice==3 && nW==0){
	  // no way to do without a W boson present
	  n_basis_problems++;
	  continue;
	}

	// find the direction of the first W boson in the lab 
	// Pick either the first W from the event record
	if (basis_choice==3){
	  TLorentzVector W_truth_lab = W_bosons[0]->P4();
	  int W_charge = W_bosons[0]->Charge;
	  // since there is not necessarily more than one W in the truth, calculate the other from the reco quantities
	  TLorentzVector W_other_lab = system_lab - W_truth_lab;
	  if (W_charge>0){
	    W_plus_lab=W_truth_lab;
	    W_minus_lab=W_other_lab;
	  }
	} else if (basis_choice==4){
	  // or the reconstructed W+ from the leptons and neutrinos	  
	  W_plus_lab = WPlusFourVectorReco;
	  W_minus_lab = WMinusFourVectorReco;
	}
	// Boost to the CM frame
	W_plus_cm = W_plus_lab; W_plus_cm.Boost(boost_lab_to_cm);
	W_minus_cm = W_minus_lab; W_minus_cm.Boost(boost_lab_to_cm);

	myplots->hist_p_cm_W->Fill( W_plus_cm.P() );
	myplots->hist_pt_cm_W->Fill( W_plus_cm.Pt() );

	double function_of_v1v2 = W_plus_cm.Gamma() * W_minus_cm.Gamma()
	  * (1 + W_plus_cm.Beta() + W_minus_cm.Beta() );
	myplots->hist_spacelike_invariant->Fill(function_of_v1v2 );
	
	if (printout){
	  cout << "W bosons in lab + then -" << endl;
	  W_plus_lab.Print();
	  W_minus_lab.Print();
	  cout << "W bosons in cm frame + then -" << endl;
	  W_plus_cm.Print();
	  W_minus_cm.Print();
	}
	
	
	// define khat to be in the direction of the W+ boson,
	TVector3 khat_vector = W_plus_cm.Vect().Unit();
	
	// the quantities eqns (6) and (7) of 2102.11883
	TVector3 phat_vector = beam_cm.Vect().Unit();
	// the quantity called y in eqns (6) and (7) of 2102.11883
	double y_quantity = phat_vector.Dot(khat_vector);
	//double r_quantity = TMath::Sqrt(1-y_quantity*y_quantity);
	TVector3 rhat_vector = ( (phat_vector - y_quantity * khat_vector) ).Unit();
	TVector3 nhat_vector = ( phat_vector.Cross(khat_vector) ).Unit(); 	

	// define the basis relative to these vectors;
	basis[0] = nhat_vector;
	basis[1] = rhat_vector;
	basis[2] = khat_vector;
	if (printout){
	  std::cout << "Basis vectors:"
		    << "\nx=k=" << khat_vector.x() << khat_vector.y() << khat_vector.z();
	  std::cout << "\ny=r=" << rhat_vector.x() << rhat_vector.y() << rhat_vector.z();
	  std::cout << "\nz=n=" << nhat_vector.x() << nhat_vector.y() << nhat_vector.z();
	  std::cout << "\nDot products: k.r: " << khat_vector.Dot(rhat_vector)
		    << " k.n: " <<  khat_vector.Dot(nhat_vector)
		    << " n.r " <<  nhat_vector.Dot(rhat_vector) << std::endl;
	}
      } else {
	cerr << "Invalid choice of basis flag " << basis_choice << endl;
	exit(1);
      } // end of selecting basis vectors
      
      // cosines for angles relative to the cartesian coordinates
      double cos_theta_1[3];
      double cos_theta_2[3];

      // for each cartesian coordinate, calculate the cosine of the angle to each basis vector
      TVector3 lepton_unit_vector_1;
      TVector3 lepton_unit_vector_2;
      // there's a choice of which frame to measure the lepton directions in as follows:
      if (lepton_frame_choice==1){
	// unit vectors in the CM frame of the higgs pointing in the direction of the selected leptons
	lepton_unit_vector_1 = lplus_cm.Vect().Unit();
	lepton_unit_vector_2 = lminus_cm.Vect().Unit();
      } else if (lepton_frame_choice==2) {
	// calculate the unit vectors of the leptons in their respective W+- parent rest frames
	TVector3 boost_cm_to_wplus  = W_plus_cm.BoostVector();
	TVector3 boost_cm_to_wminus = W_minus_cm.BoostVector();
	TLorentzVector lplus_wplus_frame   = lplus_cm;  lplus_wplus_frame.Boost(-boost_cm_to_wplus);
	TLorentzVector lminus_wminus_frame = lminus_cm; lminus_wminus_frame.Boost(-boost_cm_to_wminus);
	lepton_unit_vector_1 = lplus_wplus_frame.Vect().Unit();
	lepton_unit_vector_2 = lminus_wminus_frame.Vect().Unit();
	
	if (printout){
	  cout << "Boost cm->W+ \t"; boost_cm_to_wplus.Print();
	  cout << "Boost cm->W-\t"; boost_cm_to_wminus.Print();
	  cout << "l+ in W+ frame:" << endl;
	  lplus_wplus_frame.Print();
	  cout << "l- in W- frame:" << endl;
	  lminus_wminus_frame.Print();
	  TLorentzVector wplus_wplus_frame = W_plus_cm; wplus_wplus_frame.Boost(-boost_cm_to_wplus);
	  TLorentzVector wminus_wminus_frame = W_minus_cm; wminus_wminus_frame.Boost(-boost_cm_to_wminus);
	  cout << "W+ in W+ frame"<< endl;
	  wplus_wplus_frame.Print();
	  cout << "W- in W- frame"<< endl;
	  wminus_wminus_frame.Print();
	}
      } else{
	cerr << "\nInvalid choice of frame in which to measure lepton direction lepton_frame_choice=" << lepton_frame_choice << "\n" << endl;
	exit(1);
      }
      
      if (printout){
	cout << "Lepton unit vectors" << endl;
	lepton_unit_vector_1.Print();
	lepton_unit_vector_2.Print();
	std::cout << "Basis:" << endl;
	basis[0].Print();
	basis[1].Print();
	basis[2].Print();
      }

      // find sum of cossq theta, for purposes of determining density matrix and checking spin-1
      double sum_cossq_th_l1=0.;
      double sum_cossq_th_l2=0.;

      for (int icart=0; icart<3; icart++){
	// calculate cosines of angles
	cos_theta_1[icart] = lepton_unit_vector_1.Dot(basis[icart]);
	cos_theta_2[icart] = lepton_unit_vector_2.Dot(basis[icart]);

	// fill plots
       	myplots->hist_cos_theta[icart]->Fill(cos_theta_1[icart]);
	myplots->hist_cos_theta[icart+3]->Fill(cos_theta_2[icart]);

	// fill also squares
       	myplots->hist_cossq_theta[icart]->Fill(cos_theta_1[icart]*cos_theta_1[icart]);
	myplots->hist_cossq_theta[icart+3]->Fill(cos_theta_2[icart]*cos_theta_2[icart]);

	// add cartesian componenst to sum of cossqth
	sum_cossq_th_l1+=cos_theta_1[icart]*cos_theta_1[icart];
	sum_cossq_th_l2+=cos_theta_2[icart]*cos_theta_2[icart];
      }
      myplots->hist_sum_cossq_theta->Fill(sum_cossq_th_l1);
      myplots->hist_sum_cossq_theta->Fill(sum_cossq_th_l2);
      
      // now calculate the values of the correllation matrix
      // first calculate the 3x3 quantities in equation (10) of 2102.11883
      double epsilon_ab[3][3];
      if (printout) cout << "epsilons" << endl;
      for (int ia=0;ia<3; ia++){
	for (int ib=0; ib<3; ib++){
	  // calculate the quantity in equation (10)
	  epsilon_ab[ia][ib] = cos_theta_1[ia] * cos_theta_2[ib];
	  if (printout){
	    cout << "\t" << epsilon_ab[ia][ib];
	  }
	  myplots->hist_epsilon[ia+3*ib]->Fill(epsilon_ab[ia][ib]);
	  
	  // ** CHSH ** 
	  // contribution to CHSH integrand is added - factor of 9 from  2102.11883 eqn. (9)
	  // is replaced by factor of 4 for the spin-1 equivalent spin-average of cos(theta)
	  // for a pair of m=l=1 systems
	  double contribution = -4*epsilon_ab[ia][ib];
	  correlation_ab[ia][ib] += contribution;
	  correlation_ab_sq[ia][ib] += contribution*contribution;
	  myplots->hist_Cab_elements->Fill(ia+3*ib+1,contribution);

	  // ** CGLMP **
	  // Work out which angles to take.
	  // For x axis take cosines to y and z axes
	  const double p_x = cos_theta_1[(ia+1)%3];
	  const double p_y = cos_theta_1[(ia+2)%3];
	  const double m_x = cos_theta_2[(ib+1)%3];	  
	  const double m_y = cos_theta_2[(ib+2)%3];
	  
	  double cglmp = calculate_cglmp(p_x, p_y , m_x, m_y);
	  cglmp_contrib[ia][ib] += cglmp;
	  cglmp_contrib_sq[ia][ib] += cglmp*cglmp;
	  myplots->hist_cglmp_elements->Fill(ia+3*ib+1,cglmp);
	}
	if (printout) cout << endl;
      }

      // for normalisation
      n_correlation++;
    }
  }

  cout << endl << "SUMMARY\nW mass selected in the range [" << min_W_mass << "," <<max_W_mass << "]\n"
       <<  "Selected " << n_correlation << " of " << numberOfEntries << " events\n"<< endl;

  cout << "Using basis choice " << basis_choice << " and lepton_frame_choice "  << lepton_frame_choice << " there were " << n_basis_problems << " events for which no basis could be constructed" << endl;
  
  // horrible feature needed for filling a ROOT 3x3 Matrix array
  Double_t Cab_array[9];
  // normalise the integral of the correlation matrix
  cout << "The correlation matrix is:" << endl;
  for (int ia=0;ia<3; ia++){
    for (int ib=0; ib<3; ib++){
      cout << setprecision(5) << std::fixed;
      // Calcualte the correlation matrix and its error
      double element = correlation_ab[ia][ib] / (n_correlation) ;
      double mean_of_squares = correlation_ab_sq[ia][ib] / (n_correlation) ;
      double variance = mean_of_squares - element*element;
      double rms = TMath::Sqrt(variance);
      double error = TMath::Sqrt(variance/n_correlation); // unscaled
      double significance = (element)/(error/TMath::Sqrt(lumi_scale_factor));

      // fill the corresponding component of the array
      Cab_array[ia+ib*3]=element;
      cout << "\t\t" << element << " +- " << error << " unscaled -> ("
	   << significance << " sigma @ "<< target_lumi << "/fb )  rms = " << rms;
    }
    cout << endl;
  }

  // Fill a 3x3 matrix in order to get its eigenvalues
  cout << "Constructing the 3x3 matrix Cab" << endl;
  TMatrixD Cab(3,3);
  Cab.SetMatrixArray(Cab_array);
  Cab.Print();
  // Now form the 3x3 matrix M = C^T C
  // first make the transpose

  //cout << "... and its transpose:" << endl;
  TMatrixD CabT(3,3); CabT.Transpose(Cab);
  //CabT.Print();
  // Now form the product
  cout << "The symmetric matric c^T C is :"<<endl;
  TMatrixD myM(3,3); myM.Mult(CabT,Cab);
  myM.Print();

  cout << "... which has eigenvalues:" <<endl;

  // Find the eigenvalues
  TVectorD eigenvalues;
  myM.EigenVectors(eigenvalues);
  eigenvalues.Print();

  // Find out if we violate the Bell inequality according to eqn(4) of 2102.11883
  double sum_of_two_largest = eigenvalues[0] + eigenvalues[1];
  cout << "The sum of largest two eigenvalues is: " << sum_of_two_largest << endl << endl;
  if (sum_of_two_largest>2) {
    cout << "Opps! This is larger than two so violates Cirel'son QM upper bound of 2 - see Horodecki et al 1995" << endl << endl;
  } else if (sum_of_two_largest>1) {
    cout << "This is larger than unity so violates the generalised Bell inequality" << endl << endl;
  }

  double cglmp_array[9];
  
  std::cout << "elements of CGLMP operator:\n";
  for (int ia=0; ia<3; ia++){
    for (int ib=0; ib<3; ib++){
      cout << setprecision(5) << std::fixed;
      double element = cglmp_contrib[ia][ib] / (n_correlation) ;
      double mean_of_squares = cglmp_contrib_sq[ia][ib] / (n_correlation) ;
      double variance = mean_of_squares - element*element;
      double rms = TMath::Sqrt(variance); 
      double error = TMath::Sqrt(variance/n_correlation);
      double significance = (element-2.)/(error/TMath::Sqrt(lumi_scale_factor));
      cglmp_array[ia+ib*3]=element;
      cout << "\t\t" << element << " +- " << error << " unscaled (" << significance << " sigma for " << target_lumi << "/fb)   rms= " << rms; 
      
      inout_values.I3_cglmp[ia][ib]=element;
      inout_values.signif_cglmp[ia][ib]=significance;
    }
    std::cout << endl;
  }
  std::cout << endl;
  
  std::cout << "diagonal elements of cglmp operator: "
	    << cglmp_contrib[0][0] / n_correlation << " " 
    	    << cglmp_contrib[1][1] / n_correlation << " " 
    	    << cglmp_contrib[2][2] / n_correlation << std::endl;

  TMatrixD CGLMPab(3,3);
  CGLMPab.SetMatrixArray(cglmp_array);
  //CGLMPab.Print();
  
  TMatrixD CGLMP_ab_T(3,3); CGLMP_ab_T.Transpose(Cab);
  cout << "The symmetric matric c^T C is :"<<endl;
  TMatrixD CGLMP_M(3,3); CGLMP_M.Mult(CabT,Cab);
  
  TVectorD eigenvalues_cglmp;
  CGLMP_M.EigenVectors(eigenvalues_cglmp);
  std::cout << "Eigenvalues of C^TC for CGLMP:" << std::endl;
  eigenvalues_cglmp.Print();

  inout_values.fraction_surviving =  (double)n_correlation/(double)numberOfEntries;
  inout_values.number_surviving = (double)n_correlation * lumi_scale_factor; 
  inout_values.Sigma_chsh = sum_of_two_largest;
    
  output->cd();
  writehists(myplots);

  output->Close();

  std::cout << "Plots written to output file " << outputFile << std::endl;
};
