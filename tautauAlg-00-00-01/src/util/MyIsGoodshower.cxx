#include "tautauAlg/util/MyIsGoodshower.h"
bool isGoodShower(EvtRecTrack *recTrack) {

  IDataProviderSvc *eventSvc = NULL;
  Gaudi::svcLocator()->service("EventDataSvc", eventSvc);

  SmartDataPtr<EvtRecEvent> recEvt(eventSvc, EventModel::EvtRec::EvtRecEvent);
  SmartDataPtr<EvtRecTrackCol> recTrkCol(eventSvc,
                                         EventModel::EvtRec::EvtRecTrackCol);

  if (!recTrack->isEmcShowerValid())
    return false;
  RecEmcShower *emcTrk = recTrack->emcShower();
  double eraw = emcTrk->energy();
  double phi = emcTrk->phi();
  double the = emcTrk->theta();
  HepLorentzVector shP4(eraw * sin(the) * cos(phi), eraw * sin(the) * sin(phi),
                        eraw * cos(the), eraw);
  double cosThetaSh = shP4.vect().cosTheta();
  if (!((fabs(cosThetaSh) < 0.8 && shP4.e() > 0.025) ||
        ((fabs(cosThetaSh) > 0.86) && (fabs(cosThetaSh) < 0.92) &&
         (shP4.e() > 0.05))))
    return false;
  double time = emcTrk->time();
  if (time < 0 || time > 14)
    return false;

  Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
  double dang = 200.;

  for (int j = 0; j < recEvt->totalCharged(); j++) {
    EvtRecTrackIterator jtTrk = recTrkCol->begin() + j;
    if (!(*jtTrk)->isExtTrackValid())
      continue;
    RecExtTrack *extTrk = (*jtTrk)->extTrack();
    if (extTrk->emcVolumeNumber() == -1)
      continue;
    Hep3Vector extpos = extTrk->emcPosition();
    double angd1 = extpos.angle(emcpos);
    if (angd1 < dang)
      dang = angd1;
  }

  // 2020.1.11, revise it
  if (dang > 200)
    return false;
  dang = dang * 180 / (CLHEP::pi);
  if (fabs(dang) <= 10)
    return false;

  return true;
}
