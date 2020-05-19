#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"

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
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMuon/GlobalTrackingTools/interface/GlobalMuonTrackMatcher.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "TTree.h"

class MuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
 public:
  MuonAnalyzer(edm::ParameterSet const& iParameterSet);

  void beginJob() override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
  void endJob() override {}

  typedef std::pair<const Trajectory*, reco::TrackRef> TrackCand;

  // tree
  TTree* tree_;
  TTree* staTree_;

  // tokens
  edm::EDGetTokenT<edm::View<reco::Muon> > muonsToken_;
  edm::EDGetTokenT<edm::View<reco::GenParticle> > genParticlesToken_;
  edm::EDGetTokenT<edm::View<reco::Vertex> > verticesToken_;
  edm::EDGetTokenT<edm::View<reco::PFCandidate> > pfCandsToken_;
  edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> standAloneToken_;
  edm::EDGetTokenT<reco::TrackCollection> globalToken_;

  typedef edm::Association<reco::GenParticleCollection> GenParticleMatch;
  edm::EDGetTokenT<GenParticleMatch> trackGenMatchToken_;
  edm::EDGetTokenT<GenParticleMatch> standAloneGenMatchToken_;
  edm::EDGetTokenT<GenParticleMatch> globalGenMatchToken_;

  // muon reco stuff
  std::unique_ptr<MuonServiceProxy> theService;
  GlobalMuonTrackMatcher* theGlbMatcher;

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
    "muon_innerTrack_t0", "muon_innerTrack_t0Error",
    "muon_innerTrack_gen_pt", "muon_innerTrack_gen_eta", "muon_innerTrack_gen_phi",
    "muon_innerTrack_gen_mass", "muon_innerTrack_gen_pdgId", "muon_innerTrack_gen_index",
    "muon_outerTrack_pt", "muon_outerTrack_eta", "muon_outerTrack_phi",
    "muon_outerTrack_t0", "muon_outerTrack_t0Error",
    "muon_outerTrack_gen_pt", "muon_outerTrack_gen_eta", "muon_outerTrack_gen_phi",
    "muon_outerTrack_gen_mass", "muon_outerTrack_gen_pdgId", "muon_outerTrack_gen_index",
    "muon_InTimeMuon",
    "muon_combinedQuality_chi2LocalMomentum",
    "muon_combinedQuality_localDistance",
    "muon_combinedQuality_globalDeltaEtaPhi",
    "muon_combinedQuality_chi2LocalPosition",
    "muon_combinedQuality_chi2Time",
    "gen_pt", "gen_eta", "gen_phi", "gen_mass", "gen_pdgId",
    "gen_isPromptFinalState",
    "gen_isPromptDecayed",
    "gen_isLastCopy",
    "gen_isHardProcess",
    "gen_isDirectPromptTauDecayProductFinalState",
    "gen_fromHardProcessFinalState",
    "gen_fromHardProcessDecayed",
  };
  std::vector<std::string> staBranchNames_ = {
    "standalone_pt", "standalone_eta", "standalone_phi", "standalone_charge",
    "standalone_t0", "standalone_t0Error",
    "standalone_gen_pt", "standalone_gen_eta", "standalone_gen_phi", "standalone_gen_mass", "standalone_gen_pdgId",
    //"standalone_muonSimInfo_primaryClass",
    //"standalone_muonSimInfo_extendedClass",
    //"standalone_muonSimInfo_flavour",
    //"standalone_muonSimInfo_pdgId",
    //"standalone_muonSimInfo_pt",
    //"standalone_muonSimInfo_eta",
    //"standalone_muonSimInfo_phi",
    //"standalone_muonSimInfo_mass",
    //"standalone_muonSimInfo_tpAssoQuality",
  };
  std::vector<std::string> staTkBranchNames_ = {
    "standalone_track_pt", "standalone_track_eta", "standalone_track_phi", "standalone_track_charge",
    "standalone_track_t0", "standalone_track_t0Error",
    "standalone_track_chi2LocalMomentum",
    "standalone_track_localDistance",
    "standalone_track_globalDeltaEtaPhi",
    "standalone_track_chi2LocalPosition",
    "standalone_track_chi2Time",
    "standalone_track_genMatch",
    //"standalone_track_gen_pt",
    //"standalone_track_gen_eta",
    //"standalone_track_gen_phi",
    //"standalone_track_gen_mass",
    //"standalone_track_gen_pdgId",
  };
  std::map<std::string, std::vector<float> > branches_;
  std::map<std::string, float> staBranches_;
  std::map<std::string, std::vector<float> > staTkBranches_;
};

MuonAnalyzer::MuonAnalyzer(const edm::ParameterSet &iConfig) :
    muonsToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonsSrc"))),
    genParticlesToken_(consumes<edm::View<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genSrc"))),
    verticesToken_(consumes<edm::View<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("verticesSrc"))),
    pfCandsToken_(consumes<edm::View<reco::PFCandidate> >(iConfig.getParameter<edm::InputTag>("pfSrc"))),
    tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("trackSrc"))),
    standAloneToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("standAloneSrc"))),
    globalToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("globalSrc"))),
    trackGenMatchToken_(consumes<reco::GenParticleMatch>(iConfig.getParameter<edm::InputTag>("trackMCMatch"))),
    standAloneGenMatchToken_(consumes<reco::GenParticleMatch>(iConfig.getParameter<edm::InputTag>("standAloneMCMatch"))),
    globalGenMatchToken_(consumes<reco::GenParticleMatch>(iConfig.getParameter<edm::InputTag>("globalMCMatch")))
{

  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService = std::make_unique<MuonServiceProxy>(serviceParameters, consumesCollector());

  edm::ParameterSet trackMatcherPSet = iConfig.getParameter<edm::ParameterSet>("GlobalMuonTrackMatcher");
  theGlbMatcher = new GlobalMuonTrackMatcher(trackMatcherPSet, theService.get());

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
    branchNames_.push_back("muon_muonSimInfo_"+label+"_genIndex");
  }

  // create TTree
  usesResource("TFileService");
  edm::Service<TFileService> FS;
  tree_ = FS->make<TTree>("MuonTree", "MuonTree");
  staTree_ = FS->make<TTree>("StandAloneTree", "StandAloneTree");

  for (auto b : branchNames_) {
    branches_[b] = std::vector<float>();
    tree_->Branch(b.c_str(), &branches_[b]);
  }

  for (auto b : staBranchNames_) {
    staBranches_[b] = 0.;
    staTree_->Branch(b.c_str(), &staBranches_[b]);
  }

  for (auto b : staTkBranchNames_) {
    staTkBranches_[b] = std::vector<float>();
    staTree_->Branch(b.c_str(), &staTkBranches_[b]);
  }

}

void MuonAnalyzer::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {

  theService->update(iSetup);

  // Get the collections
  edm::Handle<edm::View<reco::Muon> > muons;
  iEvent.getByToken(muonsToken_, muons);

  edm::Handle<edm::View<reco::GenParticle> > genParticles;
  iEvent.getByToken(genParticlesToken_, genParticles);

  edm::Handle<edm::View<reco::Vertex> > vertices;
  iEvent.getByToken(verticesToken_, vertices);

  edm::Handle<edm::View<reco::PFCandidate> > pfCands;
  iEvent.getByToken(pfCandsToken_, pfCands);

  edm::Handle<reco::TrackCollection> tracks;
  iEvent.getByToken(tracksToken_, tracks);

  edm::Handle<reco::TrackCollection> standAloneTracks;
  iEvent.getByToken(standAloneToken_, standAloneTracks);

  edm::Handle<reco::TrackCollection> globalTracks;
  iEvent.getByToken(globalToken_, globalTracks);

  edm::Handle<GenParticleMatch> trackGenMatch;
  iEvent.getByToken(trackGenMatchToken_, trackGenMatch);

  edm::Handle<GenParticleMatch> standAloneGenMatch;
  iEvent.getByToken(standAloneGenMatchToken_, standAloneGenMatch);

  edm::Handle<GenParticleMatch> globalGenMatch;
  iEvent.getByToken(globalGenMatchToken_, globalGenMatch);

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
    branches_["gen_isPromptFinalState"].push_back(gen.isPromptFinalState());
    branches_["gen_isPromptDecayed"].push_back(gen.isPromptDecayed());
    branches_["gen_isLastCopy"].push_back(gen.isLastCopy());
    branches_["gen_isHardProcess"].push_back(gen.isHardProcess());
    branches_["gen_isDirectPromptTauDecayProductFinalState"].push_back(gen.isDirectPromptTauDecayProductFinalState());
    branches_["gen_fromHardProcessFinalState"].push_back(gen.fromHardProcessFinalState());
    branches_["gen_fromHardProcessDecayed"].push_back(gen.fromHardProcessDecayed());

    genIdxToOutIdxMap[idx] = jdx;

    idx++;
    jdx++;
  }

  // store standalone and tracks matches
  // in separate tree
  // keep, but redo it iterating through muon collection instead of tracks
  double chi2, d, dist, Rpos, time_chi2;
  unsigned int sdx = -1;
  unsigned int tdx;
  for (const auto sta: *standAloneTracks) {
    sdx++;

    if (sta.p()<5)
      continue;

    staBranches_["standalone_pt"] = sta.pt();
    staBranches_["standalone_eta"] = sta.eta();
    staBranches_["standalone_phi"] = sta.phi();
    staBranches_["standalone_charge"] = sta.charge();
    staBranches_["standalone_t0"] = sta.t0();
    staBranches_["standalone_t0Error"] = sta.t0Error();
    
    reco::TrackRef staRef(standAloneTracks, sdx);
    // gen match
    reco::GenParticleRef gp = (*standAloneGenMatch)[staRef];
    if (gp.isNonnull()) {
      staBranches_["standalone_gen_pt"] = gp->pt();
      staBranches_["standalone_gen_eta"] = gp->eta();
      staBranches_["standalone_gen_phi"] = gp->phi();
      staBranches_["standalone_gen_mass"] = gp->mass();
      staBranches_["standalone_gen_pdgId"] = gp->pdgId();
    } else {
      staBranches_["standalone_gen_pt"] = 0;
      staBranches_["standalone_gen_eta"] = 0;
      staBranches_["standalone_gen_phi"] = 0;
      staBranches_["standalone_gen_mass"] = 0;
      staBranches_["standalone_gen_pdgId"] = 0;
    }

    for (auto b : staTkBranchNames_) {
      staTkBranches_[b].clear();
    }
    tdx = -1;
    for (const auto tk: *tracks) {
      tdx++;

      if (reco::deltaR(sta, tk)>0.4 || tk.p()<2)
        continue;

      staTkBranches_["standalone_track_pt"].push_back(tk.pt());
      staTkBranches_["standalone_track_eta"].push_back(tk.eta());
      staTkBranches_["standalone_track_phi"].push_back(tk.phi());
      staTkBranches_["standalone_track_charge"].push_back(tk.charge());
      staTkBranches_["standalone_track_t0"].push_back(tk.t0());
      staTkBranches_["standalone_track_t0Error"].push_back(tk.t0Error());

      reco::TrackRef tkRef(tracks, tdx);
      TrackCand staCand = TrackCand((Trajectory*)nullptr, staRef);
      TrackCand tkCand = TrackCand((Trajectory*)nullptr, tkRef);
      chi2 = theGlbMatcher->match(staCand, tkCand, 0, 0);
      d = theGlbMatcher->match(staCand, tkCand, 1, 0);
      Rpos = theGlbMatcher->match(staCand, tkCand, 2, 0);
      dist = theGlbMatcher->match(staCand, tkCand, 3, 0);
      time_chi2 = theGlbMatcher->match(staCand, tkCand, 4, 0);
      staTkBranches_["standalone_track_chi2LocalMomentum"].push_back(chi2);
      staTkBranches_["standalone_track_localDistance"].push_back(d);
      staTkBranches_["standalone_track_chi2LocalPosition"].push_back(dist);
      staTkBranches_["standalone_track_globalDeltaEtaPhi"].push_back(Rpos);
      staTkBranches_["standalone_track_chi2Time"].push_back(time_chi2);

      // gen match
      reco::GenParticleRef gpt = (*trackGenMatch)[tkRef];
      if (gpt.isNonnull() && gp.isNonnull()) {
        bool match = (gp->pt()>0);
        match = (match && std::abs(gp->pt()-gpt->pt())<1e-3);
        match = (match && std::abs(gp->eta()-gpt->eta())<1e-3);
        match = (match && std::abs(gp->phi()-gpt->phi())<1e-3);
        match = (match && std::abs(gp->mass()-gpt->mass())<1e-3);
        staTkBranches_["standalone_track_genMatch"].push_back(match);
      } else {
        staTkBranches_["standalone_track_genMatch"].push_back(0);
      }
    }

    // Fill the standalone tree
    staTree_->Fill();
  }

  // store muon
  idx = 0;
  for (const auto muon: *muons) {

    //if (muon.pt()<5 || !(muon.isTrackerMuon() || muon.isStandAloneMuon() || muon.isGlobalMuon())) {
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
    if (muon.isQualityValid()) {
      branches_["muon_combinedQuality_chi2LocalMomentum"].push_back(muon.combinedQuality().chi2LocalMomentum);
      branches_["muon_combinedQuality_localDistance"].push_back(muon.combinedQuality().localDistance);
      branches_["muon_combinedQuality_globalDeltaEtaPhi"].push_back(muon.combinedQuality().globalDeltaEtaPhi);
      branches_["muon_combinedQuality_chi2LocalPosition"].push_back(muon.combinedQuality().chi2LocalPosition);
      branches_["muon_combinedQuality_chi2Time"].push_back(muon.combinedQuality().chi2Time);
    } else {
      branches_["muon_combinedQuality_chi2LocalMomentum"].push_back(0);
      branches_["muon_combinedQuality_localDistance"].push_back(0);
      branches_["muon_combinedQuality_globalDeltaEtaPhi"].push_back(0);
      branches_["muon_combinedQuality_chi2LocalPosition"].push_back(0);
      branches_["muon_combinedQuality_chi2Time"].push_back(0);
    }

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
      branches_["muon_innerTrack_t0"].push_back(muon.innerTrack()->t0());
      branches_["muon_innerTrack_t0Error"].push_back(muon.innerTrack()->t0Error());
      reco::GenParticleRef gp = (*trackGenMatch)[muon.innerTrack()];
      if (gp.isNonnull()) {
        branches_["muon_innerTrack_gen_pt"].push_back(gp->pt());
        branches_["muon_innerTrack_gen_eta"].push_back(gp->eta());
        branches_["muon_innerTrack_gen_phi"].push_back(gp->phi());
        branches_["muon_innerTrack_gen_mass"].push_back(gp->mass());
        branches_["muon_innerTrack_gen_pdgId"].push_back(gp->pdgId());
        branches_["muon_innerTrack_gen_index"].push_back(genIdxToOutIdxMap[gp.key()]);
      } else {
        branches_["muon_innerTrack_gen_pt"].push_back(0);
        branches_["muon_innerTrack_gen_eta"].push_back(0);
        branches_["muon_innerTrack_gen_phi"].push_back(0);
        branches_["muon_innerTrack_gen_mass"].push_back(0);
        branches_["muon_innerTrack_gen_pdgId"].push_back(0);
        branches_["muon_innerTrack_gen_index"].push_back(-1);
      }
    } else {
      branches_["muon_innerTrack_pt"].push_back(0);
      branches_["muon_innerTrack_eta"].push_back(0);
      branches_["muon_innerTrack_phi"].push_back(0);
      branches_["muon_innerTrack_t0"].push_back(0);
      branches_["muon_innerTrack_t0Error"].push_back(0);
      branches_["muon_innerTrack_gen_pt"].push_back(0);
      branches_["muon_innerTrack_gen_eta"].push_back(0);
      branches_["muon_innerTrack_gen_phi"].push_back(0);
      branches_["muon_innerTrack_gen_mass"].push_back(0);
      branches_["muon_innerTrack_gen_pdgId"].push_back(0);
      branches_["muon_innerTrack_gen_index"].push_back(-1);
    }

    // outerTrack
    if (muon.outerTrack().isNonnull()) {
      branches_["muon_outerTrack_pt"].push_back(muon.outerTrack()->pt());
      branches_["muon_outerTrack_eta"].push_back(muon.outerTrack()->eta());
      branches_["muon_outerTrack_phi"].push_back(muon.outerTrack()->phi());
      branches_["muon_outerTrack_t0"].push_back(muon.outerTrack()->t0());
      branches_["muon_outerTrack_t0Error"].push_back(muon.outerTrack()->t0Error());
      reco::GenParticleRef gp = (*standAloneGenMatch)[muon.outerTrack()];
      if (gp.isNonnull()) {
        branches_["muon_outerTrack_gen_pt"].push_back(gp->pt());
        branches_["muon_outerTrack_gen_eta"].push_back(gp->eta());
        branches_["muon_outerTrack_gen_phi"].push_back(gp->phi());
        branches_["muon_outerTrack_gen_mass"].push_back(gp->mass());
        branches_["muon_outerTrack_gen_pdgId"].push_back(gp->pdgId());
        branches_["muon_outerTrack_gen_index"].push_back(genIdxToOutIdxMap[gp.key()]);
      } else {
        branches_["muon_outerTrack_gen_pt"].push_back(0);
        branches_["muon_outerTrack_gen_eta"].push_back(0);
        branches_["muon_outerTrack_gen_phi"].push_back(0);
        branches_["muon_outerTrack_gen_mass"].push_back(0);
        branches_["muon_outerTrack_gen_pdgId"].push_back(0);
        branches_["muon_outerTrack_gen_index"].push_back(-1);
      }
    } else {
      branches_["muon_outerTrack_pt"].push_back(0);
      branches_["muon_outerTrack_eta"].push_back(0);
      branches_["muon_outerTrack_phi"].push_back(0);
      branches_["muon_outerTrack_t0"].push_back(0);
      branches_["muon_outerTrack_t0Error"].push_back(0);
      branches_["muon_outerTrack_gen_pt"].push_back(0);
      branches_["muon_outerTrack_gen_eta"].push_back(0);
      branches_["muon_outerTrack_gen_phi"].push_back(0);
      branches_["muon_outerTrack_gen_mass"].push_back(0);
      branches_["muon_outerTrack_gen_pdgId"].push_back(0);
      branches_["muon_outerTrack_gen_index"].push_back(-1);
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
      //branches_["muon_muonSimInfo_"+keyVal.first+"_genIndex"].push_back(msi.);
      //
      // manually find gen
      unsigned int gdx = 0;
      unsigned int matched;
      bool match = false;
      for (const auto gen : *genParticles) {

        match = msi.p4.pt()>0;
        match = (match && std::abs(msi.p4.pt()-gen.pt())<1e-3);
        match = (match && std::abs(msi.p4.eta()-gen.eta())<1e-3);
        match = (match && std::abs(msi.p4.phi()-gen.phi())<1e-3);
        match = (match && std::abs(msi.p4.mass()-gen.mass())<1e-3);
        match = (match && std::abs(msi.pdgId-gen.pdgId())<1e-3);
        if (match) {
          matched = genIdxToOutIdxMap[gdx];
          break;
        }

        gdx++;
      }
      if (match) {
        branches_["muon_muonSimInfo_"+keyVal.first+"_genIndex"].push_back(matched);
      } else {
        branches_["muon_muonSimInfo_"+keyVal.first+"_genIndex"].push_back(-1);
      }
    }

    idx++;
  }


  // Fill the tree
  tree_->Fill();

}

DEFINE_FWK_MODULE(MuonAnalyzer);
