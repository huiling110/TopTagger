// -*- C++ -*-
//
// Package:    TopTagger/TopTagger
// Class:      SHOTProducerforBSM
//
/**\class SHOTProducerforBSM SHOTProducerforBSM.cc TopTagger/TopTagger/plugins/SHOTProducerforBSM.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Nathaniel Pastika
//         Created:  Thu, 09 Nov 2017 21:29:56 GMT
//
//


// system include files
#include <memory>
#include <vector>
#include <string>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

//#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/ESHandle.h"

#include "TLorentzVector.h"

#include "TopTagger/TopTagger/interface/TopTagger.h"
#include "TopTagger/TopTagger/interface/TopTaggerResults.h"
#include "TopTagger/TopTagger/interface/TopObject.h"
#include "TopTagger/TopTagger/interface/Constituent.h"
#include "TopTagger/TopTagger/interface/TopObjLite.h"

//this include is necessary to handle exceptions thrown by the top tagger code
#include "TopTagger/CfgParser/interface/TTException.h"

class SHOTProducerforBSM : public edm::stream::EDProducer<>
{
public:
    explicit SHOTProducerforBSM(const edm::ParameterSet&);
    ~SHOTProducerforBSM();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
    template<typename LEP>
        //?not familar with template
    bool lepJetPassdRMatch(const int& jetIndex, const LEP& lepton, const edm::Handle<std::vector<pat::Jet> >& jets)
    {
        //Get the index of the jet closest to the lepton
        if(jetIndex < 0 || jetIndex >= static_cast<int>(jets->size()))
        {
            //There is no match
            return false;
        }

        //Check if it falls within the dR criterion
        const pat::Jet& jet = (*jets)[jetIndex];
        double dR = deltaR(lepton.p4(), jet.p4());
        return dR < leptonJetDr_;
    }
    std::unique_ptr<std::vector<TopObjLite>> getFinalTopObjLiteVec(const TopTaggerResults& ttr);
    std::unique_ptr<std::vector<TopObjLite>> getCandidateTopObjLiteVec(const TopTaggerResults& ttr);

    virtual void beginStream(edm::StreamID) override;
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
    virtual void endStream() override;

    // helper function to calculate subjet qg input vars
    void compute(const reco::Jet * jet, bool isReco, double& totalMult_, double& ptD_, double& axis1_, double& axis2_);

    // ----------member data ---------------------------
    edm::EDGetTokenT<std::vector<pat::Jet> > JetTok_;
    edm::EDGetTokenT<std::vector<pat::Muon> > muonTok_;
    edm::EDGetTokenT<std::vector<pat::Electron> > elecTok_;

    std::string elecIDBitName_, qgTaggerKey_, deepCSVBJetTags_, bTagKeyString_, taggerCfgFile_, elecIsoName_, muonIsoName_;
    double ak4ptCut_, leptonJetDr_, discriminatorCut_, muonPtCut_, elecPtCut_, muonIsoCut_, elecIsoCut_;
    bool doLeptonCleaning_, saveAllTopCandidates_;
    std::string elecIDFlagName_;
    reco::Muon::Selector muonIDFlag_;
//    int elecIDFlag_;

    TopTagger tt;
};


SHOTProducerforBSM::SHOTProducerforBSM(const edm::ParameterSet& iConfig)
{
    //register vector of top objects
    produces<std::vector<TopObjLite>>();
    //?where is produces defined?

    //now do what ever other initialization is needed
    edm::InputTag jetSrc = iConfig.getParameter<edm::InputTag>("ak4JetSrc");
    edm::InputTag muonSrc = iConfig.getParameter<edm::InputTag>("muonSrc");
    edm::InputTag elecSrc = iConfig.getParameter<edm::InputTag>("elecSrc");
    //so I think the difference betwen the InputTag and regular C++ type is that InputTag is a module name.

    ak4ptCut_ = iConfig.getParameter<double>("ak4PtCut");
    muonPtCut_ = iConfig.getParameter<double>("muonPtCut");
    elecPtCut_ = iConfig.getParameter<double>("elecPtCut");

    muonIsoName_ = iConfig.getParameter<std::string>("muonIsoName");
    elecIsoName_ = iConfig.getParameter<std::string>("elecIsoName");

    muonIsoCut_ = iConfig.getParameter<double>("muonIsoCut");
    elecIsoCut_ = iConfig.getParameter<double>("elecIsoCut");

    std::string muonIDFlagName = iConfig.getParameter<std::string>("muonIDFlag");
    if     (muonIDFlagName.compare("CutBasedIdLoose")  == 0) muonIDFlag_ = reco::Muon::CutBasedIdLoose;
    else if(muonIDFlagName.compare("CutBasedIdMedium") == 0) muonIDFlag_ = reco::Muon::CutBasedIdMedium;
    else if(muonIDFlagName.compare("CutBasedIdTight")  == 0) muonIDFlag_ = reco::Muon::CutBasedIdTight;
    else
    {
        muonIDFlag_ = reco::Muon::CutBasedIdLoose;
        edm::LogWarning("SHOTProducerforBSM") << "Warning!!! Unrecognized muon ID flag\"" << muonIDFlagName << "\", defaulting to loose.  (Options are CutBasedIdLoose, CutBasedIdMedium, CutBasedIdTight)" << std::endl;
    }

  //  elecIDBitName_ = iConfig.getParameter<std::string>("elecIDBitFieldName");//"VIDNestedWPBitmap"

//    std::string elecIDFlagName = iConfig.getParameter<std::string>("elecIDFlag");//"CutBasedIdVeto"
      elecIDFlagName_ = iConfig.getParameter<std::string>("elecIDFlag");//"CutBasedIdVeto"
//    if     (elecIDFlagName.compare("CutBasedIdVeto")   == 0) elecIDFlag_ = 1;
//    else if(elecIDFlagName.compare("CutBasedIdLoose")  == 0) elecIDFlag_ = 2;
//    else if(elecIDFlagName.compare("CutBasedIdMedium") == 0) elecIDFlag_ = 3;
//    else if(elecIDFlagName.compare("CutBasedIdTight")  == 0) elecIDFlag_ = 4;
//    else
//    {
//        elecIDFlag_ = 1;
//        edm::LogWarning("SHOTProducerforBSM") << "Warning!!! Unrecognized muon ID flag \"" << elecIDFlagName << "\", defaulting to veto.  (Options are CutBasedIdVeto, CutBasedIdLoose, CutBasedIdMedium, CutBasedIdTight)" << std::endl;
//    }

    leptonJetDr_ = iConfig.getParameter<double>("leptonJetDr");//0.2 from cfi, the delta R of jet and lepton

    doLeptonCleaning_  = iConfig.getParameter<bool>("doLeptonCleaning");

    qgTaggerKey_ = iConfig.getParameter<std::string>("qgTaggerKey");
    deepCSVBJetTags_ = iConfig.getParameter<std::string>("deepCSVBJetTags");
    bTagKeyString_ = iConfig.getParameter<std::string>("bTagKeyString");//'pfCombinedInclusiveSecondaryVertexV2BJetTags'

    taggerCfgFile_ = iConfig.getParameter<edm::FileInPath>("taggerCfgFile").fullPath();//"TopTagger/TopTagger/data/TopTaggerCfg-DeepResolved_DeepCSV_GR_noDisc_Release_v1.0.0/TopTagger.cfg"

    discriminatorCut_ = iConfig.getParameter<double>("discriminatorCut");

    saveAllTopCandidates_  = iConfig.getParameter<bool>("saveAllTopCandidates");

    JetTok_ = consumes<std::vector<pat::Jet> >(jetSrc);
    muonTok_ = consumes<std::vector<pat::Muon>>(muonSrc);
    elecTok_ = consumes<std::vector<pat::Electron>>(elecSrc);//'slimmedElectrons'

    //configure the top tagger
    try
    {
        //For working directory use cfg file location
        size_t splitLocation = taggerCfgFile_.rfind("/");
        std::string workingDir = taggerCfgFile_.substr(0, splitLocation);
        std::string configName = taggerCfgFile_.substr(splitLocation + 1);
        tt.setWorkingDirectory(workingDir);//tt:TopTagger class
        tt.setCfgFile(configName);
    }
    catch(const TTException& e)
        //?
    {
        //Convert the TTException into a cms::Exception
        throw cms::Exception(e.getFileName() + ":" + std::to_string(e.getLineNumber()) + ", in function \"" + e.getFunctionName() + "\" -- " + e.getMessage());
    }
}


SHOTProducerforBSM::~SHOTProducerforBSM()
{

}


//
// member functions
//
std::unique_ptr<std::vector<TopObjLite>> SHOTProducerforBSM::getFinalTopObjLiteVec(const TopTaggerResults& ttr)
{/*{{{*/
    //get reconstructed top
    const std::vector<TopObject*>& tops = ttr.getTops();

    //Translate TopObject to TopObjLite and save to event
    std::unique_ptr<std::vector<TopObjLite>> liteTops(new std::vector<TopObjLite>());
    for(auto const * const top : tops)
    {
        if(top->getDiscriminator() > discriminatorCut_)
        {
            liteTops->emplace_back(*top);
        }
    }
    return liteTops;
}/*}}}*/

std::unique_ptr<std::vector<TopObjLite>> SHOTProducerforBSM::getCandidateTopObjLiteVec(const TopTaggerResults& ttr)
{/*{{{*/
    //get top candidates
    const std::vector<TopObject>& tops = ttr.getTopCandidates();

    //Translate TopObject to TopObjLite and save to event
    std::unique_ptr<std::vector<TopObjLite>> liteTops(new std::vector<TopObjLite>());
    for(const auto& top : tops)
    {
        if(top.getDiscriminator() > discriminatorCut_)
        {
            liteTops->emplace_back(top);
        }
    }
    return liteTops;
}/*}}}*/


// ------------ method called to produce the data  ------------
void SHOTProducerforBSM::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    //Get jet collection
    edm::Handle<std::vector<pat::Jet> > jets;
    iEvent.getByToken(JetTok_, jets);

    std::set<int> jetToClean;
    if(doLeptonCleaning_)//True
    {
        //Get lepton collections
        edm::Handle<std::vector<pat::Muon> > muons;
        iEvent.getByToken(muonTok_, muons);
        edm::Handle<std::vector<pat::Electron> > elecs;
        iEvent.getByToken(elecTok_, elecs);

        //remove leptons that match to jets with a dR cone
        for(const pat::Muon& muon : *muons)
        {
            const int jetIndex = static_cast<int>(muon.userCand("jet").key());//here
            if(muon.pt() > muonPtCut_ &&
               jetIndex >= 0 && jetIndex < static_cast<int>(static_cast<int>(jets->size())) &&  //is matched to an existing jet
               muon.passed(muonIDFlag_) && muon.userFloat(muonIsoName_)/muon.pt() < muonIsoCut_)  //passes ID and Iso  //here
            {
                if(lepJetPassdRMatch(jetIndex, muon, jets)) jetToClean.insert(jetIndex);
            }
        }
        for(const pat::Electron& elec : *elecs)
        {
            const int jetIndex = static_cast<int>(elec.userCand("jet").key());//here
            //?elec.userCand("jet").key()?
            if(elec.pt() > elecPtCut_ &&
               jetIndex >= 0 && jetIndex < static_cast<int>(jets->size()))  //matched to an existing jet
            {
                //recalculate electron ID without pfRelIso
                //?why have to calculate ID without Iso?
          /*      const int NCUTS = 10;
                const int BITSTRIDE = 3;
                const int BITMASK = 0x7;
                const int ISOBITMASK = 070000000;  //note to the curious, 0 before an integer is octal, so 070 = 0x38 = 56, so this corrosponds to the three pfRelIso bits
                int elecCutBits = elec.userInt(elecIDBitName_);//"VIDNestedWPBitmap"
                int cutBits = elecCutBits | ISOBITMASK; // the | masks the iso cut
                //?
                int elecID = 07; // start with the largest 3 bit number
                for(int i = 0; i < NCUTS; ++i)
                {
                    elecID = std::min(elecID, cutBits & BITMASK);
                    cutBits = cutBits >> BITSTRIDE;
                }*/

//                if(elecID >= elecIDFlag_ && elec.userFloat(elecIsoName_)/elec.pt() < elecIsoCut_)
                if(elec.electronID(elecIDFlagName_) && elec.userFloat(elecIsoName_)/elec.pt() < elecIsoCut_)
                {
                    if(lepJetPassdRMatch(jetIndex, elec, jets)) jetToClean.insert(jetIndex);
                }
            }
        }
    }

    //container holding input jet info for top tagger
    std::vector<Constituent> constituents;

    //initialize iJet such that it is incrememted to 0 upon start of loop
    for(int iJet = 0; iJet < static_cast<int>(jets->size()); ++iJet)
    {
        const pat::Jet& jet = (*jets)[iJet];

        //Apply pt cut on jets
        if(jet.pt() < ak4ptCut_) continue;

        //Apply lepton cleaning
        if(doLeptonCleaning_ && jetToClean.count(iJet)) continue;//cout():Count elements with a specific value

        TLorentzVector perJetLVec(jet.p4().X(), jet.p4().Y(), jet.p4().Z(), jet.p4().T());

        double qgPtD = jet.userFloat("qgptD");
        double qgAxis1 = jet.userFloat("qgAxis1");
        double qgAxis2 = jet.userFloat("qgAxis2");
        double qgMult = static_cast<double>(jet.userInt("qgMult"));
        double deepCSVb = jet.bDiscriminator((deepCSVBJetTags_+":probb").c_str());
        double deepCSVc = jet.bDiscriminator((deepCSVBJetTags_+":probc").c_str());
        double deepCSVl = jet.bDiscriminator((deepCSVBJetTags_+":probudsg").c_str());
        double deepCSVbb = jet.bDiscriminator((deepCSVBJetTags_+":probbb").c_str());
        double deepCSVcc = jet.bDiscriminator((deepCSVBJetTags_+":probcc").c_str());
        double btag = jet.bDiscriminator(bTagKeyString_.c_str());
        double chargedHadronEnergyFraction = jet.chargedHadronEnergyFraction();
        double neutralHadronEnergyFraction = jet.neutralHadronEnergyFraction();
        double chargedEmEnergyFraction = jet.chargedEmEnergyFraction();
        double neutralEmEnergyFraction = jet.neutralEmEnergyFraction();
        double muonEnergyFraction = jet.muonEnergyFraction();
        double photonEnergyFraction = jet.photonEnergyFraction();
        double electronEnergyFraction = jet.electronEnergyFraction();
        double recoJetsHFHadronEnergyFraction = jet.HFHadronEnergyFraction();
        double recoJetsHFEMEnergyFraction = jet.HFEMEnergyFraction();
        double chargedHadronMultiplicity = jet.chargedHadronMultiplicity();
        double neutralHadronMultiplicity = jet.neutralHadronMultiplicity();
        double photonMultiplicity = jet.photonMultiplicity();
        double electronMultiplicity = jet.electronMultiplicity();
        double muonMultiplicity = jet.muonMultiplicity();

        constituents.emplace_back(perJetLVec, btag, 0.0);
        constituents.back().setIndex(iJet);
        constituents.back().setExtraVar("qgMult"                              , qgMult);
        constituents.back().setExtraVar("qgPtD"                               , qgPtD);
        constituents.back().setExtraVar("qgAxis1"                             , qgAxis1);
        constituents.back().setExtraVar("qgAxis2"                             , qgAxis2);
        constituents.back().setExtraVar("recoJetschargedHadronEnergyFraction" , chargedHadronEnergyFraction);
        constituents.back().setExtraVar("recoJetschargedEmEnergyFraction"     , chargedEmEnergyFraction);
        constituents.back().setExtraVar("recoJetsneutralEmEnergyFraction"     , neutralEmEnergyFraction);
        constituents.back().setExtraVar("recoJetsmuonEnergyFraction"          , muonEnergyFraction);
        constituents.back().setExtraVar("recoJetsHFHadronEnergyFraction"      , recoJetsHFHadronEnergyFraction);
        constituents.back().setExtraVar("recoJetsHFEMEnergyFraction"          , recoJetsHFEMEnergyFraction);
        constituents.back().setExtraVar("recoJetsneutralEnergyFraction"       , neutralHadronEnergyFraction);
        constituents.back().setExtraVar("PhotonEnergyFraction"                , photonEnergyFraction);
        constituents.back().setExtraVar("ElectronEnergyFraction"              , electronEnergyFraction);
        constituents.back().setExtraVar("ChargedHadronMultiplicity"           , chargedHadronMultiplicity);
        constituents.back().setExtraVar("NeutralHadronMultiplicity"           , neutralHadronMultiplicity);
        constituents.back().setExtraVar("PhotonMultiplicity"                  , photonMultiplicity);
        constituents.back().setExtraVar("ElectronMultiplicity"                , electronMultiplicity);
        constituents.back().setExtraVar("MuonMultiplicity"                    , muonMultiplicity);
        constituents.back().setExtraVar("DeepCSVb"                            , deepCSVb);
        constituents.back().setExtraVar("DeepCSVc"                            , deepCSVc);
        constituents.back().setExtraVar("DeepCSVl"                            , deepCSVl);
        constituents.back().setExtraVar("DeepCSVbb"                           , deepCSVbb);
        constituents.back().setExtraVar("DeepCSVcc"                           , deepCSVcc);
    }

    //run top tagger
    try
    {
        tt.runTagger(std::move(constituents));
    }
    catch(const TTException& e)
    {
        //Convert the TTException into a cms::Exception
        throw cms::Exception(e.getFileName() + ":" + std::to_string(e.getLineNumber()) + ", in function \"" + e.getFunctionName() + "\" -- " + e.getMessage());
    }

    //retrieve the top tagger results object
    const TopTaggerResults& ttr = tt.getResults();

    //save tops or top candidates
    if(saveAllTopCandidates_) iEvent.put(std::move( getCandidateTopObjLiteVec(ttr) ));
    else                      iEvent.put(std::move( getFinalTopObjLiteVec(ttr)     ));
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void SHOTProducerforBSM::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void SHOTProducerforBSM::endStream()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
//  After this has been created, when scram build is run a cfi file will be automatically generated
void SHOTProducerforBSM::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SHOTProducerforBSM);
