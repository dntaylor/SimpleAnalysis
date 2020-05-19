/** \class NewMCTrackMatcher
 *
 * \author Luca Lista, INFN
 *
 *
 */
#include "DataFormats/Common/interface/Association.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "DataFormats/RecoCandidate/interface/TrackAssociation.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

namespace edm {
  class ParameterSet;
}

using namespace edm;
using namespace std;
using namespace reco;

class NewMCTrackMatcher : public edm::stream::EDProducer<> {
public:
  /// constructor
  NewMCTrackMatcher(const edm::ParameterSet &);

private:
  void produce(edm::Event &evt, const edm::EventSetup &es) override;
  edm::EDGetTokenT<reco::SimToRecoCollection> associatorMap_s2r_;
  edm::EDGetTokenT<reco::RecoToSimCollection> associatorMap_r2s_;
  edm::EDGetTokenT<edm::View<reco::Track> > tracks_;
  edm::EDGetTokenT<GenParticleCollection> genParticles_;
  edm::EDGetTokenT<TrackingParticleCollection> trackingParticles_;
  typedef edm::Association<reco::GenParticleCollection> GenParticleMatch;
  typedef edm::Association<edm::View<reco::Track> > TrackMatch;
};

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

NewMCTrackMatcher::NewMCTrackMatcher(const ParameterSet &p)
    : associatorMap_s2r_(consumes<reco::SimToRecoCollection>(p.getParameter<InputTag>("associator"))),
      associatorMap_r2s_(consumes<reco::RecoToSimCollection>(p.getParameter<InputTag>("associator"))),
      tracks_(consumes<edm::View<reco::Track> >(p.getParameter<InputTag>("tracks"))),
      genParticles_(consumes<GenParticleCollection>(p.getParameter<InputTag>("genParticles"))),
      trackingParticles_(consumes<TrackingParticleCollection>(p.getParameter<InputTag>("trackingParticles"))) {
  produces<GenParticleMatch>();
  //produces<TrackMatch>();
}

void NewMCTrackMatcher::produce(Event &evt, const EventSetup &es) {
  Handle<reco::RecoToSimCollection> r2s;
  evt.getByToken(associatorMap_r2s_, r2s);

  Handle<reco::SimToRecoCollection> s2r;
  evt.getByToken(associatorMap_s2r_, s2r);

  Handle<View<Track> > tracks;
  evt.getByToken(tracks_, tracks);

  Handle<TrackingParticleCollection> trackingParticles;
  evt.getByToken(trackingParticles_, trackingParticles);

  Handle<GenParticleCollection> genParticles;
  evt.getByToken(genParticles_, genParticles);

  unique_ptr<GenParticleMatch> matchGen(new GenParticleMatch(GenParticleRefProd(genParticles)));
  //unique_ptr<TrackMatch> matchTrack(new TrackMatch(edm::RefProd<edm::View<reco::Track> >(tracks)));

  GenParticleMatch::Filler fillerGen(*matchGen);
  //TrackMatch::Filler fillerTrack(*matchTrack);

  size_t n = tracks->size();
  vector<int> indices(n, -1);
  for (size_t i = 0; i < n; ++i) {
    RefToBase<Track> track(tracks, i);
    RecoToSimCollection::const_iterator f = r2s->find(track);
    if (f != r2s->end()) {
      TrackingParticleRef tp = f->val.front().first;
      TrackingParticle::genp_iterator j, b = tp->genParticle_begin(), e = tp->genParticle_end();
      for (j = b; j != e; ++j) {
        const reco::GenParticle *p = j->get();
        if (p->status() == 1) {
          indices[i] = j->key();
          break;
        }
      }
    }
  }
  fillerGen.insert(tracks, indices.begin(), indices.end());
  fillerGen.fill();

  //size_t nt = trackingParticles->size();
  //size_t ng = genParticles->size();
  //vector<int> indicesGen(ng, -1);
  //for (size_t i = 0; i < nt; ++i) {
  //  TrackingParticleRef tp(trackingParticles, i);
  //  SimToRecoCollection::const_iterator f = s2r->find(tp);
  //  if (f != s2r->end()) {
  //    const RefToBase<Track> t = f->val.front().first;
  //    TrackingParticle::genp_iterator j, b = tp->genParticle_begin(), e = tp->genParticle_end();
  //    for (j = b; j != e; ++j) {
  //      const reco::GenParticle *p = j->get();
  //      if (p->status() == 1) {
  //        indicesGen[j->key()] = t.key();
  //        break;
  //      }
  //    }
  //  }
  //}
  //fillerTrack.insert(genParticles, indicesGen.begin(), indicesGen.end());
  //fillerTrack.fill();

  evt.put(std::move(matchGen));
  //evt.put(std::move(matchTrack));
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(NewMCTrackMatcher);
