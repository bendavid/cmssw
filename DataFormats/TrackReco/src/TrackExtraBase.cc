#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

using namespace reco;

unsigned int TrackExtraBase::firstRecHit() const {
  //return stored index only if original collection is available
  if (m_hitCollection.isAvailable()) {
    return m_firstHit;
  }
  else if (m_hitCollection.isThinnedAvailable(m_firstHit, m_hitCollection.productGetter())) {
    unsigned int thinnedKey;
    m_hitCollection.getThinnedProductPtr(typeid(TrackingRecHit), thinnedKey, m_hitCollection.productGetter());
    return thinnedKey;
  }
  return m_firstHit;
}
