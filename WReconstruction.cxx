#include "TLorentzVector.h"
#include <TMath.h>

std::vector<TLorentzVector> solveForNeutrinoEta(TLorentzVector* lepton, 
                                                TLorentzVector* Whad, 
                                                double nu_eta, 
                                                double mHiggs, 
                                                double mWlep) {

  double nu_cosh = cosh(nu_eta);
  double nu_sinh = sinh(nu_eta);

  double Wmass2   = mWlep*mWlep;
  double Whadmass = Whad->M();
  double Elprime  = lepton->E() * nu_cosh - lepton->Pz() * nu_sinh;
  double Ebprime  = Whad->E()   * nu_cosh - Whad->Pz()   * nu_sinh;

  double A = (lepton->Py() * Ebprime - Whad->Py() * Elprime) / (Whad->Px() * Elprime - lepton->Px() * Ebprime);
  double B = (Elprime * (mHiggs * mHiggs - Wmass2 - Whadmass * Whadmass - 2. * lepton->Dot(*Whad)) - Ebprime * Wmass2) / (2. * (lepton->Px() * Ebprime - Whad->Px() * Elprime));

  double par1 = (lepton->Px() * A + lepton->Py()) / Elprime;
  double C    = A * A + 1. - par1 * par1;
  double par2 = (Wmass2 / 2. + lepton->Px() * B) / Elprime;
  double D    = 2. * (A * B - par2 * par1);
  double F    = B * B - par2 * par2;
  double det  = D * D - 4. * C * F;

  std::vector<TLorentzVector> sol;

  ///-- 0 solutions case --///
  if (det < 0.0){
    return std::move(sol);
  }

  ///-- Only one real solution case --///
  if (det == 0.) {
    double py1 = -D / (2. * C);
    double px1 = A * py1 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pz1 = sqrt(pT2_1) * nu_sinh;

    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    if (!TMath::IsNaN(a1.E()) )
      sol.push_back(a1);
    return std::move(sol);
  }

  ///-- 2 solutions case --///
  if(det > 0){
    double tmp   = sqrt(det) / (2. * C);
    double py1   = -D / (2. * C) + tmp;
    double py2   = -D / (2. * C) - tmp;
    double px1   = A * py1 + B;
    double px2   = A * py2 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pT2_2 = px2 * px2 + py2 * py2;
    double pz1   = sqrt(pT2_1) * nu_sinh;
    double pz2   = sqrt(pT2_2) * nu_sinh;
    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));
    TLorentzVector a2(px2, py2, pz2, sqrt(pT2_2 + pz2 * pz2));

    if (!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
      sol.push_back(a1);
      sol.push_back(a2);
    }
    return std::move(sol);
  }
  
  ///-- Should never reach this point --///
  return std::move(sol);
}


      
      
int WReconstruction(FourMomentum lepton,FourMomentum parton_w_had, FourMomentum nu){

  m_best_NW_weight = -1000;
  particle_reco_w_lep = new TLorentzVector(0.,0.,0.,0.);
  particle_reco_nu= new TLorentzVector(0.,0.,0.,0.);
  //Implementing only the linear mass range variation, for now
  std::vector<double> wmass_points;
  int wmass_linear_nbins = 50;
  double wmass_linear_min   = 0;
  double wmass_linear_max   = 50;
  double wmass_stepsize = (wmass_linear_max - wmass_linear_min)/wmass_linear_nbins;
  for (int i=0; i< wmass_linear_nbins; i++){
      wmass_points.push_back(wmass_linear_min + i*wmass_stepsize);
  }

  // Nu eta linear
  std::vector<double> nueta_linear_points;
  int nueta_linear_nbins = 50;
  double nueta_linear_min   = -3.0;
  double nueta_linear_max   = +3.0;
  double nueta_stepsize = (nueta_linear_max - nueta_linear_min)/nueta_linear_nbins;
  for (int i=0; i< nueta_linear_nbins; i++){
      nueta_linear_points.push_back(nueta_linear_min + i*nueta_stepsize + 0.001);
  }


  TLorentzVector* lepton_lv= new TLorentzVector();
  TLorentzVector* parton_w_had_lv= new TLorentzVector();
  lepton_lv->SetPtEtaPhiM(lepton.pt(), lepton.eta(), convert_phi(lepton.phi()), lepton.mass());
  parton_w_had_lv->SetPtEtaPhiM(parton_w_had.pt(), parton_w_had.eta(), convert_phi(parton_w_had.phi()), parton_w_had.mass());   

  double weight = -99.;
  double weight_for_hist = -99;
  double temp_weight_ex = 0;
  double temp_weight_ey = 0;
  double temp_weight =0;

  for(double w_mass : wmass_points){
      for (double nu_eta : nueta_linear_points){
          std::vector<TLorentzVector> neutrinos;
          neutrinos = solveForNeutrinoEta(lepton_lv , parton_w_had_lv, nu_eta, 125., w_mass);
          
          weight_for_hist = -99;
          temp_weight_ex = 0;
          temp_weight_ey = 0;
          temp_weight =0;

          for (auto neutrino : neutrinos){

              TLorentzVector temp_w_boson = *(lepton_lv) + neutrino;
              TLorentzVector* particle_nu = new TLorentzVector();
              particle_nu->SetPtEtaPhiE(nu.pt(), nu.eta(), convert_phi(nu.phi()), nu.E());
              temp_weight_ex = exp(-pow((neutrino.Px() - particle_nu->Px()),2) / 2*pow(0.5,2));
              temp_weight_ey = exp(-pow((neutrino.Py() - particle_nu->Py()),2) / 2*pow(0.5,2));

              temp_weight    = temp_weight_ex*temp_weight_ey;                 
                      
              if (temp_weight > weight){
                weight = temp_weight;
                particle_reco_w_lep->SetPtEtaPhiM(temp_w_boson.Pt(), temp_w_boson.Eta(), temp_w_boson.Phi(), temp_w_boson.M());
                particle_reco_nu->SetPtEtaPhiM(neutrino.Pt(), neutrino.Eta(), neutrino.Phi(), neutrino.M()); 
                m_best_NW_weight=weight;
              }  

              if (temp_weight > weight_for_hist){
                weight_for_hist = temp_weight;
              }  
              }

              if (weight_for_hist == -99){
                weight_for_hist = 0.0;
              }
              
              
          }

      
}

return 1;

}