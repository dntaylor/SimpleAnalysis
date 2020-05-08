#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/MuonSimInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Association.h"

#include "TTree.h"

class MuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
 public:
    MuonAnalyzer(edm::ParameterSet const& iParameterSet);
 
    void beginJob() override {}
    void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
    void endJob() override {}

    // tree
    TTree* tree_;

    // tokens
    edm::EDGetTokenT<edm::View<reco::Muon> > muonsToken_;
    edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
    edm::EDGetTokenT<edm::View<reco::Vertex> > verticesToken_;
    edm::EDGetTokenT<edm::View<reco::PFCandidate> > pfCandsToken_;

    // MuonSimHit
    std::map<std::string, edm::EDGetTokenT<edm::ValueMap<reco::MuonSimInfo> > > muonSimInfoTokens_;
    std::map<std::string, edm::EDGetTokenT<edm::Association<reco::GenParticleCollection> > > muonSimInfoGenPrimaryTokens_;

    // branches
    std::vector<std::string> branchNames_ = { "muon_pt", "muon_eta", "muon_phi", "muon_mass", "muon_charge",
      "muon_isGlobalMuon", "muon_isStandAloneMuon", "muon_isTrackerMuon", "muon_isPFMuon",
      "muon_muonTime_nDof", "muon_muonTime_time", "muon_muonTime_timeError",
      "muon_rpcTime_nDof", "muon_rpcTime_time", "muon_rpcTime_timeError",
      "muon_pf_time", "muon_pf_timeError",
      "muon_innerTrack_pt", "muon_innerTrack_eta", "muon_innerTrack_phi",
      "muon_outerTrack_pt", "muon_outerTrack_eta", "muon_outerTrack_phi",
      "muon_InTimeMuon",
      "gen_pt", "gen_eta", "gen_phi", "gen_mass", "gen_pdgId",
    };
    std::map<std::string, std::vector<float> > branches_;
};

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet &iConfig) :
    muonsToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonsSrc"))),
    genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genSrc"))),
    verticesToken_(consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("verticesSrc"))),
    pfCandsToken_(consumes<edm::View<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("pfSrc")))
{

  // load the various muon simhit collections
  std::vector<edm::ParameterSet> msiConfig = iConfig.getParameter<std::vector<edm::ParameterSet> >("MuonSimInfoConfiguration");
  for (auto pset : msiConfig) {
    std::string label = pset.getParameter<std::string>("label");
    muonSimInfoTokens_[label] = consumes<edm::ValueMap<reco::MuonSimInfo> >(pset.getParameter<edm::InputTag>("muonSimInfo"));
    muonSimInfoGenPrimaryTokens_[label] = consumes<edm::Association<reco::GenParticleCollection> >(pset.getParameter<edm::InputTag>("muonSimInfoGenPrimary"));
    branchNames_.push_back("muon_muonSimInfo_"+label+"_primaryClass");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_extendedClass");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_flavour");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_pdgId");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_pt");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_eta");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_phi");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_mass");
    branchNames_.push_back("muon_muonSimInfo_"+label+"_tpAssoQuality");
    //branchNames_.push_back("muon_muonSimInfo_"+label+"_genIndex");
  }

  // create TTree
  usesResource("TFileService");
  edm::Service<TFileService> FS;
  tree_ = FS->make<TTree>("MuonTree", "MuonTree");

  for (auto b : branchNames_) {
    branches_[b] = std::vector<float>();
    tree_->Branch(b.c_str(), &branches_[b]);
  }

}

void MuonAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
    // Get the collections
    edm::Handle<edm::View<reco::Muon> > muons;
    iEvent.getByToken(muonsToken_, muons);

    edm::Handle<edm::View<reco::GenParticle> > genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    edm::Handle<edm::View<reco::Vertex> > vertices;
    iEvent.getByToken(verticesToken_, vertices);

    edm::Handle<edm::View<reco::PFCandidate> > pfCands;
    iEvent.getByToken(pfCandsToken_, pfCands);

    std::map<std::string, edm::Handle<edm::ValueMap<reco::MuonSimInfo> > > muonSimInfos;
    for (auto keyVal : muonSimInfoTokens_) {
      muonSimInfos[keyVal.first] = edm::Handle<edm::ValueMap<reco::MuonSimInfo> >();
      iEvent.getByToken(muonSimInfoTokens_[keyVal.first], muonSimInfos[keyVal.first]);
    }

    std::map<std::string, edm::Handle<edm::Association<reco::GenParticleCollection> > > muonSimMatches;
    for (auto keyVal : muonSimInfoGenPrimaryTokens_) {
      muonSimMatches[keyVal.first] = edm::Handle<edm::Association<reco::GenParticleCollection> >();
      iEvent.getByToken(muonSimInfoGenPrimaryTokens_[keyVal.first], muonSimMatches[keyVal.first]);
    }

    // clear the tree
    for (auto b : branchNames_) {
      branches_[b].clear();
    }

    // fill a muon ref to pf map
    std::map<unsigned int, unsigned int> muonRefToPF;
    unsigned int idpf = 0;
    for (const auto pf: *pfCands) {
      if (abs(pf.pdgId()) == 13 ) {
        muonRefToPF[pf.muonRef().key()] = idpf;
      }
      ++idpf;
    }

    // store gen
    unsigned int idx = 0;
    unsigned int jdx = 0;
    std::map<unsigned int, unsigned int> genIdxToOutIdxMap;
    for (const auto gen: *genParticles) {
      if (gen.status()!=1 || gen.pt()<5) {
        idx++;
        continue;
      }
      branches_["gen_pt"].push_back(gen.pt());
      branches_["gen_eta"].push_back(gen.eta());
      branches_["gen_phi"].push_back(gen.phi());
      branches_["gen_mass"].push_back(gen.mass());
      branches_["gen_pdgId"].push_back(gen.pdgId());

      genIdxToOutIdxMap[idx] = jdx;
      std::cout << "gen " << idx << " pos " << jdx << std::endl;

      idx++;
      jdx++;
    }

    // store muon
    idx = 0;
    for (const auto muon: *muons) {

      // only store global muons for now
      if (muon.pt()<5) {
        ++idx;
        continue;
      }

      edm::RefToBase<reco::Muon> muonRefTmp = muons->refAt(idx);
      reco::CandidateBaseRef muonBaseRef(muonRefTmp);

      branches_["muon_pt"].push_back(muon.pt());
      branches_["muon_eta"].push_back(muon.eta());
      branches_["muon_phi"].push_back(muon.phi());
      branches_["muon_mass"].push_back(muon.mass());
      branches_["muon_charge"].push_back(muon.charge());

      branches_["muon_isGlobalMuon"].push_back(muon.isGlobalMuon());
      branches_["muon_isStandAloneMuon"].push_back(muon.isStandAloneMuon());
      branches_["muon_isTrackerMuon"].push_back(muon.isTrackerMuon());
      branches_["muon_isPFMuon"].push_back(muon.isPFMuon());

      branches_["muon_muonTime_nDof"].push_back(muon.time().nDof);
      branches_["muon_muonTime_time"].push_back(muon.time().timeAtIpInOut);
      branches_["muon_muonTime_timeError"].push_back(muon.time().timeAtIpInOutErr);
      branches_["muon_rpcTime_nDof"].push_back(muon.rpcTime().nDof);
      branches_["muon_rpcTime_time"].push_back(muon.rpcTime().timeAtIpInOut);
      branches_["muon_rpcTime_timeError"].push_back(muon.rpcTime().timeAtIpInOutErr);

      // reference to pf
      auto iter = muonRefToPF.find(idx);
      if (iter != muonRefToPF.end()) {
        const auto& pfMuon = pfCands->at(iter->second);
        branches_["muon_pf_time"].push_back(pfMuon.time());
        branches_["muon_pf_timeError"].push_back(pfMuon.timeError());
        muonRefToPF.erase(iter);
      } else {
        branches_["muon_pf_time"].push_back(0);
        branches_["muon_pf_timeError"].push_back(0);
      }

      // innerTrack
      if (muon.innerTrack().isNonnull()) {
        branches_["muon_innerTrack_pt"].push_back(muon.innerTrack()->pt());
        branches_["muon_innerTrack_eta"].push_back(muon.innerTrack()->eta());
        branches_["muon_innerTrack_phi"].push_back(muon.innerTrack()->phi());
      } else {
        branches_["muon_innerTrack_pt"].push_back(0);
        branches_["muon_innerTrack_eta"].push_back(0);
        branches_["muon_innerTrack_phi"].push_back(0);
      }

      // outerTrack
      if (muon.outerTrack().isNonnull()) {
        branches_["muon_outerTrack_pt"].push_back(muon.outerTrack()->pt());
        branches_["muon_outerTrack_eta"].push_back(muon.outerTrack()->eta());
        branches_["muon_outerTrack_phi"].push_back(muon.outerTrack()->phi());
      } else {
        branches_["muon_outerTrack_pt"].push_back(0);
        branches_["muon_outerTrack_eta"].push_back(0);
        branches_["muon_outerTrack_phi"].push_back(0);
      }

      // selectors
      branches_["muon_InTimeMuon"].push_back(muon.passed(reco::Muon::InTimeMuon));

      // the muon sim hit info
      for (auto keyVal : muonSimInfos) {
        const auto& msi = (*keyVal.second)[muonBaseRef];
        const auto& gp = (*muonSimMatches[keyVal.first])[muonBaseRef];
        branches_["muon_muonSimInfo_"+keyVal.first+"_primaryClass"].push_back(msi.primaryClass);
        branches_["muon_muonSimInfo_"+keyVal.first+"_extendedClass"].push_back(msi.extendedClass);
        branches_["muon_muonSimInfo_"+keyVal.first+"_flavour"].push_back(msi.flavour);
        branches_["muon_muonSimInfo_"+keyVal.first+"_pdgId"].push_back(msi.pdgId);
        branches_["muon_muonSimInfo_"+keyVal.first+"_pt"].push_back(msi.p4.pt());
        branches_["muon_muonSimInfo_"+keyVal.first+"_eta"].push_back(msi.p4.eta());
        branches_["muon_muonSimInfo_"+keyVal.first+"_phi"].push_back(msi.p4.phi());
        branches_["muon_muonSimInfo_"+keyVal.first+"_mass"].push_back(msi.p4.mass());
        branches_["muon_muonSimInfo_"+keyVal.first+"_tpAssoQuality"].push_back(msi.tpAssoQuality);
        // not available, so forget about it
        //branches_["muon_muonSimInfo_"+keyVal.first+"_genIndex"].push_back(-1);
      }

      idx++;
    }


    // Fill the tree
    tree_->Fill();

}

DEFINE_FWK_MODULE(MuonAnalyzer);
