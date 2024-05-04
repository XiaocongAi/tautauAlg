#include "tautauAlg/util/MyIsGoodtrack.h"

bool isGoodTrack(EvtRecTrack *trk, Hep3Vector xorigin) {
  if (!(trk->isMdcTrackValid()))
    return false;

  RecMdcTrack *mdcTrk = trk->mdcTrack();
  HepVector a = mdcTrk->helix();
  HepSymMatrix Ea = mdcTrk->err();
  HepPoint3D pivot(0., 0., 0.);
  HepPoint3D IP(xorigin[0], xorigin[1], xorigin[2]);
  VFHelix helixp(pivot, a, Ea);
  helixp.pivot(IP);
  HepVector vec = helixp.a();
  double vrl = -100.;
  vrl = vec[0];
  double vzl = -100.;
  vzl = vec[3];
  double costheta = -100.;
  costheta = cos(mdcTrk->theta());

  if (mdcTrk->charge() == 0)
    return false;

  if (fabs(vrl) < 1 && fabs(vzl) < 10 && fabs(costheta) < 0.93)
    return true;
  return false;
}
