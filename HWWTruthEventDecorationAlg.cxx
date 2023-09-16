// PhysicsxAODCode includes
#include "HWWTruthEventDecorationAlg.h"
#include "../Root/TruthTypeUtils.h"
// EDM includes
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthEventContainer.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "HWWProcessIDEnumDef.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>
std::string vector_feature(TLorentzVector P4){ //access vector element 
    std::array<double,9>features={ P4.E(),P4.Px(),P4.Py(),P4.Pz(),P4.M(),P4.Pt(),P4.PseudoRapidity(),P4.Phi(),P4.Rapidity()};
    std::string output = "";
    for (size_t i = 0; i< 9 ;i++){                                              
        std::string val = std::to_string(features[i]);
        output+=val;
        output+=",";
        }
    output+="\n";
    return output;
    }
std::string particle_feature(const xAOD::TruthParticle* particle){ //access pointer                                             
    std::array<double,9>features={ particle->e(),particle->px(),particle->py(),particle->pz(),particle->m(),particle->pt(),particle->eta(),particle->phi(),particle->rapidity()};
    std::string output = "";
    for (size_t i = 0; i< 9 ;i++){                                              
        std::string val = std::to_string(features[i]);
        output+=val;
        output+=",";
        }
    output+="\n";
    return output;
    }
double CGLMP(double px, double mx, double py, double my){//direction cosines of l+(px,py), l-(mx,my)
    double oev = 0;
    oev = (8/sqrt(3))*(px*mx + py*my) + 25*((px*px - py*py)*(mx*mx - my*my)) + 100*(px*py*mx*my);
    return oev;
    }
double B3_Max(double a, double b, double c){
    double tmp = 0;
    double val = 0;
    tmp = std::max(a,b);
    val = std::max(tmp, c);
    return val;
    }
//declare output logs
std::ofstream Wp("WpBoson.csv"); // 1
std::ofstream Wm("WmBoson.csv"); // 2
std::ofstream Lp("LeptonP.csv"); // 3
std::ofstream Lm("LeptonM.csv"); // 4
std::ofstream diLep("diLepton.csv"); // 5
std::ofstream LLep("LeadLepton.csv"); // 6
std::ofstream sLLep("subLeadLepton.csv"); // 7
std::ofstream LeadPair("LeadingPair.csv"); // 8
std::ofstream MissET("MET.csv"); // 9
std::ofstream LpNu("LpNu.csv"); // 10
std::ofstream LmNu("LmNu.csv"); // 11
std::ofstream HBoson("Higgs.csv"); // 12
std::ofstream CGLMP_xi_p_m_xyz("xi_p_m_xyz.csv"); // 13,14,15,16,17,18
/*
std::ofstream CGLMP_xip_xy("xip_xy.csv"); // 13
std::ofstream CGLMP_xip_yz("xip_yz.csv"); // 14
std::ofstream CGLMP_xip_zx("xip_zx.csv"); // 15
std::ofstream CGLMP_xim_xy("xim_xy.csv"); // 16
std::ofstream CGLMP_xim_yz("xim_yz.csv"); // 17
std::ofstream CGLMP_xim_zx("xim_zx.csv"); // 18
*/
std::ofstream CGLMP_BXYZ("Bxy_yz_zx.csv"); // 19,20,21
/*
std::ofstream CGLMP_BXY("Bxy.csv"); // 19
std::ofstream CGLMP_BYZ("Byz.csv"); // 20
std::ofstream CGLMP_BZX("Bzx.csv"); // 21
*/
HWW::TruthEventDecorationAlg::TruthEventDecorationAlg( const std::string& name, ISvcLocator* pSvcLocator ) :
AthAnalysisAlgorithm( name, pSvcLocator ),
m_inElCont(""),
m_inMuCont(""), 
m_inMetCont(""), 
m_metObj(""),
m_VGamORTool(nullptr),
m_VGamORDeco(""),
m_WZProcessID(false)
    {
    declareProperty("InputElectrons",  m_inElCont,         "The name of the input electrons container" );
    declareProperty("InputMuons",      m_inMuCont,         "The name of the input muons container" );
    declareProperty("InputMet",        m_inMetCont,        "The name of the input met container" );
    declareProperty("MissingETObject", m_metObj,           "The name of the input met object" );
    declareProperty("VGamORTool",      m_VGamORTool,       "The VGammaORTool" );
    declareProperty("VGamORDeco",      m_VGamORDeco,       "The name of the VGammaORTool's output decoration" );
    declareProperty("RunWZProcessID",  m_WZProcessID,      "Run the WZ ProcessID algorithm" );
    }

HWW::TruthEventDecorationAlg::~TruthEventDecorationAlg() {
    //close output logs
    Wp.close(); // 1
    Wm.close(); // 2
    Lp.close(); // 3
    Lm.close(); // 4
    diLep.close(); // 5
    LLep.close(); // 6
    sLLep.close(); // 7
    LeadPair.close(); // 8
    MissET.close(); // 9
    LpNu.close(); // 10
    LmNu.close(); // 11
    HBoson.close(); // 12
    CGLMP_xi_p_m_xyz.close(); // 13
    /*
    CGLMP_xip_xy.close(); // 13
    CGLMP_xip_yz.close(); // 14
    CGLMP_xip_zx.close(); // 15
    CGLMP_xim_xy.close(); // 16
    CGLMP_xim_yz.close(); // 17
    CGLMP_xim_zx.close(); // 18
    */
    CGLMP_BXYZ.close(); // 19
    /*
    CGLMP_BXY.close(); // 19
    CGLMP_BYZ.close(); // 20
    CGLMP_BZX.close(); // 21
    */
    }

StatusCode HWW::TruthEventDecorationAlg::initialize(){
    ATH_MSG_DEBUG ("Initializing " << name() << "...");
    return StatusCode::SUCCESS;
    }
StatusCode HWW::TruthEventDecorationAlg::finalize(){
    ATH_MSG_DEBUG ("Finalizing " << name() << "...");
    return StatusCode::SUCCESS;
    }
StatusCode HWW::TruthEventDecorationAlg::execute(){
    const xAOD::EventInfo* eventInfo = nullptr;
    CHECK( evtStore()->retrieve(eventInfo) );
    // Only run HWWTruthEventDecorationAlg.cxx on the MC
    const bool isSim = eventInfo->eventType(xAOD::EventInfo::EventType::IS_SIMULATION);
    if ( !isSim ) {
        ATH_MSG_DEBUG ("It is a data event... nothing to be done...");
        return StatusCode::SUCCESS;
    }
    //auxiliary data to write into PxAOD 
    static SG::AuxElement::Decorator<double> CGLMP_xiP_x("xi_p_n");
    static SG::AuxElement::Decorator<double> CGLMP_xiP_y("xi_p_r");
    static SG::AuxElement::Decorator<double> CGLMP_xiP_z("xi_p_k");
    static SG::AuxElement::Decorator<double> CGLMP_xiM_x("xi_m_n");
    static SG::AuxElement::Decorator<double> CGLMP_xiM_y("xi_m_r");
    static SG::AuxElement::Decorator<double> CGLMP_xiM_z("xi_m_k");
    static SG::AuxElement::Decorator<double> CGLMPLepXYDeco("LepXY");
    static SG::AuxElement::Decorator<double> CGLMPLepYZDeco("LepYZ");
    static SG::AuxElement::Decorator<double> CGLMPLepZXDeco("LepZX");
    static SG::AuxElement::Decorator<double> NuPz("NeutrinoPx");
    static SG::AuxElement::Decorator<double> NuPy("NeutrinoPy");
    static SG::AuxElement::Decorator<double> NuPx("NeutrinoPz");
    static SG::AuxElement::Decorator<double> NuE("NeutrinoE");
    static SG::AuxElement::Decorator<double> antiNuPx("antiNeutrinoPx");
    static SG::AuxElement::Decorator<double> antiNuPy("antiNeutrinoPy");
    static SG::AuxElement::Decorator<double> antiNuPz("antiNeutrinoPz");
    static SG::AuxElement::Decorator<double> antiNuE("antiNeutrinoE");
    const xAOD::TruthParticle* WpBoson(nullptr);// #1 W+
    const xAOD::TruthParticle* WmBoson(nullptr);// #2 W-
    const xAOD::TruthParticle* LeptonP(nullptr);// #3 l+
    const xAOD::TruthParticle* LeptonM(nullptr);// #4 l-
    TLorentzVector llP4; // #5
    const xAOD::TruthParticle* LeadLep = nullptr; // #6
    const xAOD::TruthParticle* subLeadLep = nullptr; // #7
    int LPairLead = 0; // #8
    const xAOD::MissingET* met = nullptr; // 9 
    const xAOD::TruthParticle* LepNu(nullptr);// 10 v
    const xAOD::TruthParticle* anti_LepNu(nullptr);// 11 v_
    const xAOD::TruthParticle* Higgs(nullptr);// 12H
    double XiP[3]; // #13 #14 #15
    double XiM[3]; // #16 #17 #18
    double Bxy = 0.0; // 19
    double Byz = 0.0; // #20
    double Bzx = 0.0; // #21
    size_t NElectrons = 0;
    size_t NMuons = 0;
    // i.e. before any truth PAOD processing if it was scheduled from processTruth=True)
    //ATH_CHECK(evtStore()->retrieve(inCont,inName))
    /* 
    //Alternate method
    const xAOD::TruthEventContainer* truthEvents = nullptr; //the truthElectron/Muon is empty in TruthEventContainer
    ATH_CHECK(evtStore()->retrieve(truthEvents, "TruthEvents"));
    //std::vector<const xAOD::Jet*> jets; 
    std::vector<const xAOD::TruthParticle*> truthElectrons;
    std::vector<const xAOD::TruthParticle*> truthMuons;
    std::vector<const xAOD::TruthParticle*> truthPhotons;
    TruthTypeUtils::fillTruthParticleVectors(truthEvents, truthElectrons, truthMuons, truthPhotons);
     */
    const xAOD::TruthParticleContainer* truthElectrons(nullptr);
    ATH_CHECK(evtStore()->retrieve(truthElectrons, "TruthElectrons"));
    const xAOD::TruthParticleContainer* truthMuons(nullptr);
    ATH_CHECK(evtStore()->retrieve(truthMuons, "TruthMuons"));
    NElectrons = truthElectrons->size();
    NMuons = truthMuons->size();
    //TruthTypeUtils::fillMinPtJetVector(jets); //Get jets that pass min pt selection
    // Retrieve the truth neutrinos
    const xAOD::TruthParticleContainer* truthNeutrinos(nullptr);
    ATH_CHECK(evtStore()->retrieve(truthNeutrinos, "TruthNeutrinos"));
    // Retrieve the truth Bosons
    const xAOD::TruthParticleContainer* truthBosons(nullptr);
    ATH_CHECK(evtStore()->retrieve(truthBosons, "TruthBoson"));
    // Retrieve MET
    const xAOD::MissingETContainer* METcont = nullptr;
    if(!(m_inMetCont.value().empty())){                                              
        ATH_CHECK(evtStore()->retrieve(METcont, m_inMetCont.value()));
        met = (*METcont)[m_metObj.value()];
        }
    bool noTauNu = true;
    bool noHiggs = true;
    bool DF = false;
    bool hasLeadLep = false;
    bool hasSubLeadLep = false;
    std::vector<const xAOD::TruthParticle*> e_Nu;
    std::vector<const xAOD::TruthParticle*> anti_e_Nu;
    std::vector<const xAOD::TruthParticle*> mu_Nu;
    std::vector<const xAOD::TruthParticle*> anti_mu_Nu;
    for(size_t i = 0; i < truthNeutrinos->size(); i++){ //categorize neutrinos
        const xAOD::TruthParticle* part1 = truthNeutrinos->at(i); 
        if(part1->absPdgId() == 16){
            noTauNu = false;
            }
        if(part1->absPdgId() == 12){//electron neutrinos
            if(part1->pdgId() == 12){//electron neutrinos
                e_Nu.push_back(part1);
                }
            else if(part1->pdgId() == -12){//anti electron neutrinos
                anti_e_Nu.push_back(part1);
                }
            }
        if(part1->absPdgId() ==14){ //muon neutrinos
            if(part1->pdgId() == 14){ //muon neutrinos
                mu_Nu.push_back(part1);
                }
            else if(part1->pdgId() == -14){ //anti muon neutrinos
                anti_mu_Nu.push_back(part1);
                }
            }
        }
    ////////////////////////////////////////
    // W+/- lepton neutrino pair search //
    ///////////////////////////////////////
    if(truthBosons->size() < 2){                                              
        ATH_MSG_DEBUG ("No Truth Bosons found.");
        }
    if(noTauNu == true ){ //check and exclude tau processes 
        for (size_t i = 0; i < truthBosons->size(); i++){ //loop over all bosons
            const xAOD::TruthParticle* partA = truthBosons->at(i);
            double tmp = 0; //variable to store invariant mass deviation for comparison
            double delta = 999999; //variable to quantify deviation of invariant mass
            switch(partA->pdgId()){
                case 25:{
                    noHiggs = false;
                    Higgs = partA;
                    break; 
                    }
                case 24:{//Pair W+ boson with l+ v
                    WpBoson = partA;
                    for(size_t j =0; j<partA->nChildren(); j ++){ //loop over all children of W+ boson
                        const xAOD::TruthParticle* partB = partA->child(j);
                        switch(partB->pdgId()){ //check case for e+ or u+
                            case -11:{//case for e+ 
                                for(size_t k = 0; k < e_Nu.size(); k++){//loop over e_Nu container
                                    const xAOD::TruthParticle* partC = e_Nu.at(k); 
                                    tmp = (WpBoson->p4() - partB->p4() - partC->p4()).M(); //calculate difference in invariant mass
                                    if(tmp <= delta){ //if this e+ v pair has a smaller difference
                                        LepNu = partC; //store this electron neutrino as l+ v candidate
                                        delta = tmp; //update value for next comparison
                                        }
                                    }
                                tmp = (WpBoson->p4() - partB->p4() - LepNu->p4()).M(); //calculate difference in invariant mass
                                if(tmp <= delta){ //if this l+ v pair has a smaller difference
                                    LeptonP = partB; //store this lepton as l -v comparison
                                    delta = tmp; //update value for next comparison
                                    }   
                                break;
                                }
                            case -13:{//case for u+
                                for(size_t k = 0; k < mu_Nu.size(); k++){ //loop over mu_Nu container
                                    const xAOD::TruthParticle* partC = mu_Nu.at(k); 
                                    tmp = (WpBoson->p4() - partB->p4() - partC->p4()).M(); //calculate difference in invariant mass
                                    if(tmp <= delta){ //if this u+ v pair has a smaller difference
                                        LepNu = partC; //store this muon neutrino as u+ v candidate
                                        delta = tmp; //update value for next comparison
                                        }
                                    }
                                tmp = (WpBoson->p4() - partB->p4() - LepNu->p4()).M(); //calculate difference in invariant mass
                                if(tmp <= delta){ //if this l+ v pair has a smaller difference
                                    LeptonP = partB; //store this lepton as l -v comparison
                                    delta = tmp; //update value for next comparison
                                    }   
                                break;
                                }
                            default:{
                                ATH_MSG_DEBUG ("W+ children not found");
                                break;
                                }
                            }
                        }
                    delta = 999999;
                    break;
                    }
                case -24:{ //Pair W- boson with l- anti_v
                    WmBoson = partA;
                    for(size_t j =0; j<partA->nChildren(); j ++){ //loop over all children of W- boson
                        const xAOD::TruthParticle* partB = partA->child(j);
                        switch(partB->pdgId()){ //check case for e- or u-
                            case  11:{//case for e-
                                for(size_t k = 0; k < anti_e_Nu.size(); k++){ //loop over anti_e_Nu container
                                    const xAOD::TruthParticle* partC = anti_e_Nu.at(k); 
                                    tmp = (WmBoson->p4() - partB->p4() - partC->p4()).M(); //calculate difference in invariant mass
                                    if(tmp <= delta){ //if this e- anti_v pair has a smaller difference
                                        anti_LepNu = partC; //store this anti electron neutrino as e- anti_v candidate
                                        delta = tmp; //update value for next comparison
                                        }
                                    }
                                tmp = (WmBoson->p4() - partB->p4() - anti_LepNu->p4()).M(); //calculate difference in invariant mass
                                if(tmp <= delta){ //if this l- anti_v pair has a smaller difference 
                                    LeptonM = partB; //store this anti lepton as l- anti_v pair candidate
                                    delta = tmp;
                                    }
                                break;
                                }
                            case  13:{//case for u-
                                for(size_t k = 0; k < anti_mu_Nu.size(); k++){ //loop over anti_mu_Nu container
                                    const xAOD::TruthParticle* partC = anti_mu_Nu.at(k); 
                                    tmp = (WmBoson->p4() - partB->p4() - partC->p4()).M(); //calculate difference in invariant mass
                                    if(tmp <= delta){ //if this u- anti_v pair has a smaller difference
                                        anti_LepNu = partC; //store this anti muon neutrino as u- anti_v candidate
                                        delta = tmp; //update value for next comparision
                                        }
                                    }
                                tmp = (WmBoson->p4() - partB->p4() - anti_LepNu->p4()).M(); //calculate difference in invariant mass
                                if(tmp <= delta){ //if this l- anti_v pair has a smaller difference 
                                    LeptonM = partB; //store this anti lepton as l- anti_v pair candidate
                                    delta = tmp;
                                    }
                                break;
                                }
                            default:{
                                ATH_MSG_DEBUG ("W- children not found");
                                break;
                                }
                            }
                        }
                    delta = 999999;
                    break;
                    }
                default:{
                    ATH_MSG_DEBUG("Truth Boson Missing");
                    break;
                    }
                }
            }
    //////////////////////////////////////////
    // end W+/- lepton neutrino pair search //
    //////////////////////////////////////////
    if(LeptonP->absPdgId() == LeptonM->absPdgId()){                                              
        DF = false;
        }
    else if(LeptonP->absPdgId() != LeptonM->absPdgId()){                                              
        DF = true;
        ///////////////////////////
        // Leading Lepton search //
        ///////////////////////////
        size_t lead = 0; // always use size_t when dealing with memory range
        size_t sublead = 1; 
        double LeadLepPT = 0;
        double subLeadLepPT = 0;
        if(NElectrons + NMuons >=2){ 
            if(NElectrons >= 2 && NMuons >=2){
                const xAOD::TruthParticle* e1 = truthElectrons->at(lead);
                const xAOD::TruthParticle* e2 = truthElectrons->at(sublead);
                const xAOD::TruthParticle* m1 = truthMuons->at(lead);
                const xAOD::TruthParticle* m2 = truthMuons->at(sublead);
                //assume lepton and muon container is pT sorted. store first two elements of each
                std::array<double,4>lepPts={e1->pt()*(-1),m1->pt()*(-1) ,e2->pt()*(-1),m2->pt()*(-1)};//invert ascending sort
                std::array<const xAOD::TruthParticle*,4>leps={e1,m1,e2,m2};
                std::sort(lepPts.begin(), lepPts.end()); //descending sort
                LeadLepPT = lepPts[0]*(-1); //revert value
                subLeadLepPT = lepPts[1]*(-1);
                for(size_t i = 0; i < 4 ; i++){ //match pt value with particle
                    if(LeadLepPT == leps[i]->pt()){
                        LeadLep = leps[i]; //store particle with leading pT
                        }
                    else if(subLeadLepPT == leps[i]->pt()){
                        subLeadLep = leps[i]; //store particle with subleading pT
                        }
                    }
                }
            else if (NMuons == 1 && NElectrons == 1){
                if(truthElectrons->at(lead)->pt() > truthMuons->at(lead)->pt()){
                    LeadLep = truthElectrons->at(lead);
                    subLeadLep = truthMuons->at(lead);
                    }
                else{                                              
                    subLeadLep = truthElectrons->at(lead);
                    LeadLep = truthMuons->at(lead);
                    }
                }  
            else if(NElectrons < 2 ) { 
                if(truthMuons->at(sublead)->pt() < truthElectrons->at(lead)->pt()){
                    if(truthElectrons->at(lead)->pt() > truthMuons->at(lead)->pt()){
                        LeadLep = truthElectrons->at(lead);
                        subLeadLep = truthMuons->at(lead);
                        }
                    else{                                              
                        subLeadLep = truthElectrons->at(lead);
                        LeadLep = truthMuons->at(lead);
                        }
                    }
                else{
                    LeadLep = truthMuons->at(lead);
                    subLeadLep = truthMuons->at(sublead);
                    }
                }
            else if(NMuons < 2 ) { 
                if(truthElectrons->at(sublead)->pt() < truthMuons->at(lead)->pt()){
                    if(truthMuons->at(lead)->pt() > truthElectrons->at(lead)->pt()){
                        subLeadLep = truthElectrons->at(lead);
                        LeadLep = truthMuons->at(lead);
                        }
                    else{                                              
                        LeadLep = truthElectrons->at(lead);
                        subLeadLep = truthMuons->at(lead);
                        }
                    }
                else{                                              
                    LeadLep = truthElectrons->at(lead);
                    subLeadLep = truthElectrons->at(sublead);
                    }    
                }
            if(LeadLep->barcode() == LeptonP->barcode() || LeadLep->barcode() == LeptonM->barcode()){ //match barcode with lepton pair in H-WW final state 
                hasLeadLep = true;
                LPairLead = 1;
                hasSubLeadLep = true;
                LPairLead = 1;
                }
            if(hasLeadLep == true && hasSubLeadLep == true){                                              
                LPairLead = 2;
                }
            }
            
        ///////////////////////////////
        // end Leading Lepton search //
        ///////////////////////////////
        }

        
    NuPx(*eventInfo) = LepNu->p4().X();   
    NuPy(*eventInfo) = LepNu->p4().Y();   
    NuPz(*eventInfo) = LepNu->p4().Z();   
    NuE(*eventInfo)  = LepNu->p4().E();   
    antiNuPx(*eventInfo) =anti_LepNu->p4().X();   
    antiNuPy(*eventInfo) =anti_LepNu->p4().Y();  
    antiNuPz(*eventInfo) =anti_LepNu->p4().Z();  
    antiNuE(*eventInfo)  =anti_LepNu->p4().E(); 
    TLorentzVector WpP4, WmP4, lpP4, lmP4, HiggsP4;
    llP4 = LeptonP->p4() + LeptonM->p4(); // 5
    WpP4 = WpBoson->p4();
    WmP4 = WmBoson->p4();
    lpP4 = LeptonP->p4();
    lmP4 = LeptonM->p4();

    /* 
    double diLeptonPT = sqrt(pow(llP4.Px(),2) + pow(llP4.Py(),2))
    double diLeptonMass = llP4.M();
    double diLeptonPhi = abs(LeptonP.phi() - LeptonM.phi());
    double dPhiMETdiLep = abs(METPhi- diLeptonPhi); 
    double dPhiMETlLep;
    double dPhiMETslLep; 
    */

    /////////////////////////////////////////////
    //Bell inequality study - CGLMP calculation//
    /////////////////////////////////////////////  
    if(noHiggs == false){                                              
        HiggsP4 = Higgs->p4();
        }
    else if(noHiggs == true){                                              
        HiggsP4 = WpP4 + WmP4 ;
        }
    TLorentzVector WpP4_hrf, WmP4_hrf, lpP4_hrf, lmP4_hrf, beam;
    TVector3 p;
    beam.SetVect(TVector3(0,0,1));
    beam.SetE(1); //create a lightlight vector as the beam vector (0,0,1,1) (px,py,pz,E)
    TVector3 HiggsBoost = HiggsP4.BoostVector();
    beam.Boost(-HiggsBoost);//beam vector as seen in higgs rest frame
    p.SetXYZ(beam.Px(), beam.Py(), beam.Pz());
    p = p.Unit();
    //Boost W+ W- l+ l- to Higgs Rest frame
    WpP4_hrf = WpP4;
    WpP4_hrf.Boost(-HiggsBoost);
    WmP4_hrf = WmP4;
    WmP4_hrf.Boost(-HiggsBoost);
    lpP4_hrf = lpP4;
    lpP4_hrf.Boost(-HiggsBoost);
    lmP4_hrf = lmP4;
    lmP4_hrf.Boost(-HiggsBoost);
    //construct {n, r, k}  basis -> {x, y, z}
    TVector3 k,r,n;
    k.SetXYZ(WpP4_hrf.Px(), WpP4_hrf.Py(), WpP4_hrf.Pz());
    k = k.Unit(); 
    double_t y = p.Dot(k);
    r = (p - y*k).Unit(); //corresponds to y
    n = (p.Cross(k).Unit()); //corresponds to x
    //boost to W+ W- rest frame
    TVector3 WpBoost, WmBoost;//boost vector // to k
    WpBoost = WpP4_hrf.BoostVector();
    WmBoost = WmP4_hrf.BoostVector();
    lpP4_hrf.Boost(-WpBoost);
    lmP4_hrf.Boost(-WmBoost);
    TVector3 lp_Wprf, lm_Wmrf;
    lp_Wprf.SetXYZ(lpP4_hrf.Px(), lpP4_hrf.Py(), lpP4_hrf.Pz());
    lm_Wmrf.SetXYZ(lmP4_hrf.Px(), lmP4_hrf.Py(), lmP4_hrf.Pz());
    //{n, r, k} projection
    //calculate directional cosines
    XiP[0] = cos(lp_Wprf.Angle(n));XiP[1] = cos(lp_Wprf.Angle(r));XiP[2] = cos(lp_Wprf.Angle(k));
    XiM[0] = cos(lm_Wmrf.Angle(n));XiM[1] = cos(lm_Wmrf.Angle(r));XiM[2] = cos(lm_Wmrf.Angle(k));
    Bxy = CGLMP(XiP[0], XiM[0], XiP[1], XiM[1]);
    Byz = CGLMP(XiP[1], XiM[1], XiP[2], XiM[2]);
    Bzx = CGLMP(XiP[2], XiM[2], XiP[0], XiM[0]);
    ///////////////////////////
    // end CGLMP Calculation //
    ///////////////////////////
    }
    ////////////////////
    // csv log output //
    ////////////////////
    if(noTauNu == true && DF == true ){ 
        Wp<<particle_feature(WpBoson); // #1
        Wm<<particle_feature(WmBoson); // #2
        Lp<<particle_feature(LeptonP); // #3
        Lm<<particle_feature(LeptonM); // #4
        diLep<<vector_feature(llP4); // #5
        LLep<<particle_feature(LeadLep); // #6
        sLLep<<particle_feature(subLeadLep); // #7
        LeadPair<<LPairLead<<",\n"; // #8
        MissET<<std::to_string(met->met())+","+std::to_string(met->mpx())+","+std::to_string(met->mpy())+","+std::to_string(met->phi())+",\n"; // #9
        LpNu<<particle_feature(LepNu); // #10
        LmNu<<particle_feature(anti_LepNu); // #11
        HBoson<<particle_feature(Higgs); // #12
        CGLMP_xi_p_m_xyz << XiP[0] << ","<< XiP[1] << ","<< XiP[2] << ","<< XiM[0] << ","<< XiM[1] << ","<< XiM[2] << ",\n"; // #13,14,15,16,17,18
        /*
        CGLMP_xip_xy << XiP[0] << ",\n"; // #13
        CGLMP_xip_yz << XiP[1] << ",\n"; // #14
        CGLMP_xip_zx << XiP[2] << ",\n"; // #15
        CGLMP_xim_xy << XiM[0] << ",\n"; // #16
        CGLMP_xim_yz << XiM[1] << ",\n"; // #17
        CGLMP_xim_zx << XiM[2] << ",\n"; // #18
        */
        CGLMP_BXYZ << Bxy << ","<< Byz << ","<< Bzx << ",\n"; // #19,20,21
        /*
        CGLMP_BXY << Bxy << ",\n"; // #19
        CGLMP_BYZ << Byz << ",\n"; // #20
        CGLMP_BZX << Bzx << ",\n"; // #21
        */
        CGLMP_xiP_x(*eventInfo) = XiP[0];
        CGLMP_xiP_y(*eventInfo) = XiP[1];
        CGLMP_xiP_z(*eventInfo) = XiP[2];
        CGLMP_xiM_x(*eventInfo) = XiM[0];
        CGLMP_xiM_y(*eventInfo) = XiM[1];
        CGLMP_xiM_z(*eventInfo) = XiM[2];
        CGLMPLepXYDeco(*eventInfo) = Bxy;
        CGLMPLepYZDeco(*eventInfo) = Byz;
        CGLMPLepZXDeco(*eventInfo) = Bzx;
        }
    else if(noTauNu == false || DF == false ){ // csv log null output                                              
        std::string NotANumber = "-99999,-99999,-99999,-99999,-99999,-99999,-99999,-99999,-99999,-99999,\n";
        double XP[3] = {-99999,-99999,-99999};
        double XM[3] = {-99999,-99999,-99999};
        Bxy = -99999;
        Byz = -99999;
        Bzx = -99999;
        Wp   << NotANumber; // #1
        Wm   << NotANumber; // #2
        Lp   << NotANumber; // #3
        Lm   << NotANumber; // #4
        diLep<< NotANumber; // #5
        LLep << NotANumber; // #6
        sLLep<< NotANumber; // #7
        LeadPair<< "-99999,\n"; // #8
        MissET<< "-99999,-99999,-99999,-99999,\n"; // #9
        LpNu<< NotANumber; // #10
        LmNu<< NotANumber; // #11
        HBoson<<NotANumber; // #12
        CGLMP_xi_p_m_xyz << XP[0] << ","<< XP[1] << ","<< XP[2] << ","<< XM[0] << ","<< XM[1] << ","<< XM[2] << ",\n"; // #13,14,15,16,17,18
        /*
        CGLMP_xip_xy << XiP[0] << ","; // #13
        CGLMP_xip_yz << XiP[1] << ","; // #14
        CGLMP_xip_zx << XiP[2] << ","; // #15
        CGLMP_xim_xy << XiM[0] << ","; // #16
        CGLMP_xim_yz << XiM[1] << ","; // #17
        CGLMP_xim_zx << XiM[2] << ","; // #18
        */
        CGLMP_BXYZ << Bxy << ","<< Byz << ","<< Bzx << ",\n"; // #19,20,21
        /*
        CGLMP_BXY << Bxy << ","; // #19
        CGLMP_BYZ << Byz << ","; // #20
        CGLMP_BZX << Bzx << ","; // #21
        */
        CGLMP_xiP_x(*eventInfo) = XP[0];
        CGLMP_xiP_y(*eventInfo) = XP[1];
        CGLMP_xiP_z(*eventInfo) = XP[2];
        CGLMP_xiM_x(*eventInfo) = XM[0];
        CGLMP_xiM_y(*eventInfo) = XM[1];
        CGLMP_xiM_z(*eventInfo) = XM[2];
        CGLMPLepXYDeco(*eventInfo) = Bxy;
        CGLMPLepYZDeco(*eventInfo) = Byz;
        CGLMPLepZXDeco(*eventInfo) = Bzx;
        }
        
        ////////////////////////
        // end csv log output //
        ////////////////////////
    /*	//Mass cut
    for (int i=0; i<=5;i++){
        if(std::min(WpP4.M(),WmP4.M()) > 100001){
            XY[i] += Bxy; YZ[i] += Byz; ZX[i] += Bzx;
            Bxyz_N[i] ++;
            Bell_CGLMP << "("<<i*10<<" GeV cut) Bxy = " << Bxy << std::endl;
            Bell_CGLMP << "("<<i*10<<" GeV cut) Byz = " << Byz << std::endl;
            Bell_CGLMP << "("<<i*10<<" GeV cut) Bzx = " << Bzx << std::endl;
            Bell_CGLMP << "("<<i*10<<" GeV cut) Bxy average: " << XY[i]/Bxyz_N[i] << " Byz average: " << YZ[i]/Bxyz_N[i] << " Bzx average: " << ZX[i]/Bxyz_N[i] << " Counts: " << Bxyz_N[i]<< std::endl;
            XY2[i] += Bxy*Bxy; YZ2[i] += Byz*Byz; ZX2[i] += Bzx*Bzx;
            Bell_CGLMP << "("<<i*10<<" GeV cut) Bxy squared average: " << XY2[i]/Bxyz_N[i] << " Byz squared average: " << YZ2[i]/Bxyz_N[i] << " Bzx squared average: " << ZX2[i]/Bxyz_N[i] <<  std::endl;
            Bell_CGLMP << "("<<i*10<<" GeV cut) Bxy mean squared error: " << sqrt((XY2[i]/Bxyz_N[i]) - (XY[i]/Bxyz_N[i])*(XY[i]/Bxyz_N[i])) << " Byz mean squared error: " << sqrt((YZ2[i]/Bxyz_N[i]) - (YZ[i]/Bxyz_N[i])*(YZ[i]/Bxyz_N[i])) << " Bzx mean squared error: " << sqrt((ZX2[i]/Bxyz_N[i]) - (ZX[i]/Bxyz_N[i])*(ZX[i]/Bxyz_N[i])) << std::endl;
            }
        }
    */
    return StatusCode::SUCCESS;
}
