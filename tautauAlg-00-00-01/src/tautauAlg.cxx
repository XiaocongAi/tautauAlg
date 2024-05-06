//-----------------------------
// e+ e- -> tau+ tau-,
//          |    |-> e- nu anti-nu
//          mu+ nu anti_nu
//----------------------------

#include "tautauAlg/util/MyConst.h"
#include "tautauAlg/util/MyUtil.h"
#include "tautauAlg/util/MyIsGoodshower.h"
#include "tautauAlg/util/MyIsGoodtrack.h"
#include "tautauAlg/util/MyInitIP.h"
#include "tautauAlg/util/MyInfoShower.h"
#include "tautauAlg/util/MyInfoCharge.h"
#include "tautauAlg/util/MyGamConv.h"
#include "tautauAlg/util/Mygamee.h"
#include "tautauAlg/tautauAlg.h"

#include "GaudiKernel/IHistogramSvc.h"

using namespace Event;

typedef std::vector<int> Vint;
typedef std::vector<HepLorentzVector> Vp4;
int Ncut[10];

//***************************Declare******************************
tautauAlg::tautauAlg(const std::string &name, ISvcLocator *pSvcLocator)
    : Algorithm(name, pSvcLocator) {
  m_nEvtDisp = 100;
  m_nEvt = 0;
  declareProperty("writeGenOnly", m_writeGenOnly = false);
  declareProperty("Ecms", m_ecms = 4.23);
  declareProperty("testMC", m_testMC = 1);
  declareProperty("gamNumCut", m_gamNumCut = 1);
  declareProperty("pi0Chi2Cut", m_pi0Chi2Cut = 100);
  declareProperty("pi0NumCut", m_pi0NumCut = 1);
}

//***************************Fill Tree******************************
StatusCode tautauAlg::initialize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in initialize()" << endmsg;

  StatusCode status;
  // add m_cutflow
  status = service("THistSvc", m_histSvc);
  if (status.isFailure()) {
    log << MSG::INFO << "Unable to retrieve pointer to THistSvc" << endreq;
    return status;
  }
  m_cutflow = new TH1F("cutflow", "cutflow", 6, -0.5, 5.5); //-0.5——5.5,bin = 6
  // m_recEmcHitMap= new TH2F("recEmcHitMap", "recEmcHitMap", 43, 0, 43, 119, 0,
  // 119); //see
  // DetectorDescription/Identifier/Identifier-00-02-17/Identifier/EmcID.h
  status = m_histSvc->regHist("/HIST/cutflow", m_cutflow);
  // status = m_histSvc->regHist("/HIST/recEmcHitMap", m_recEmcHitMap);

  NTuplePtr nt1(ntupleSvc(), "FILE1/gen");
  if (nt1)
    m_tuple1 = nt1;
  else {
    m_tuple1 = ntupleSvc()->book("FILE1/gen", CLID_ColumnWiseTuple,
                                 "         theory verfication");
    if (m_tuple1) {
      m_tuple1->addItem("within_kine_bounds", m_within_kine_bounds);
      m_tuple1->addItem("pion_px", m_pi_px);
      m_tuple1->addItem("pion_py", m_pi_py);
      m_tuple1->addItem("pion_pz", m_pi_pz);
      m_tuple1->addItem("pion_e", m_pi_e);

      m_tuple1->addItem("gam_px", m_gam_px);
      m_tuple1->addItem("gam_py", m_gam_py);
      m_tuple1->addItem("gam_pz", m_gam_pz);
      m_tuple1->addItem("gam_e", m_gam_e);

      m_tuple1->addItem("tau_px", m_tau_px);
      m_tuple1->addItem("tau_py", m_tau_py);
      m_tuple1->addItem("tau_pz", m_tau_pz);
      m_tuple1->addItem("tau_e", m_tau_e);

      m_tuple1->addItem("nu_px", m_nu_px);
      m_tuple1->addItem("nu_py", m_nu_py);
      m_tuple1->addItem("nu_pz", m_nu_pz);
      m_tuple1->addItem("nu_e", m_nu_e);

      m_tuple1->addItem("x", m_x);
      m_tuple1->addItem("y", m_y);
      m_tuple1->addItem("t", m_t);
      m_tuple1->addItem("s", m_s);
      m_tuple1->addItem("E", m_E);
      m_tuple1->addItem("z", m_z);

      m_tuple1->addItem("f_IB", m_f_IB);
      m_tuple1->addItem("f_VV", m_f_VV);
      m_tuple1->addItem("f_AV", m_f_AV);
      m_tuple1->addItem("f_AA", m_f_AA);
      m_tuple1->addItem("f_IB_V", m_f_IB_V);
      m_tuple1->addItem("f_IB_A", m_f_IB_A);
      m_tuple1->addItem("Gam_IB", m_Gam_IB);
      m_tuple1->addItem("Gam_SD", m_Gam_SD);
      m_tuple1->addItem("Gam_int", m_Gam_int);
      m_tuple1->addItem("ngamfromTau", m_ngamfromTau, 0, 1000000);
      m_tuple1->addIndexedItem("gam_motherIndexMC", m_ngamfromTau,
                               m_gam_motherIndex);
      m_tuple1->addItem("nefrommu", m_nefrommu);
      m_tuple1->addItem("nefromtau", m_nefromtau);
    } else {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple1)
          << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt2(ntupleSvc(), "FILE1/mctruth");
  if (nt2)
    m_tuple2 = nt2;
  else {
    m_tuple2 = ntupleSvc()->book("FILE1/mctruth", CLID_ColumnWiseTuple,
                                 "mctruth/all particle");
    if (m_tuple2) {
      m_tuple2->addItem("run", m_run);
      m_tuple2->addItem("event", m_event);
      m_tuple2->addItem("indexmc", m_idxmc, 0, 100);
      m_tuple2->addIndexedItem("pdgid", m_idxmc, m_pdgid);
      m_tuple2->addIndexedItem("motheridx", m_idxmc, m_motheridx);
      m_tuple2->addItem("info_tru", 7, 8, m_info_tru);
    } else {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple2)
          << endmsg;
      return StatusCode::FAILURE;
    }
  }

  NTuplePtr nt3(ntupleSvc(), "FILE1/rec");
  if (nt3)
    m_tuple3 = nt3;
  else {
    m_tuple3 = ntupleSvc()->book("FILE1/rec", CLID_ColumnWiseTuple,
                                 " variables at reco-level");
    if (m_tuple3) {
      m_tuple3->addItem("rec_run", m_rec_run);
      m_tuple3->addItem("rec_event", m_rec_event);
      m_tuple3->addItem("rec_nchrg", m_rec_nchrg);
      m_tuple3->addItem("rec_nneu", m_rec_nneu);

      m_tuple3->addItem("indexmc", m_rec_idxmc, 0, 100);
      m_tuple3->addIndexedItem("pdgid", m_rec_idxmc, m_rec_pdgid);
      m_tuple3->addIndexedItem("motheridx", m_rec_idxmc, m_rec_motheridx);

      // MDC information
      m_tuple3->addItem("charge_1", m_charge_1);
      m_tuple3->addItem("charge_2", m_charge_2);
      m_tuple3->addItem("pmdc_1", m_pmdc_1);
      m_tuple3->addItem("pmdc_2", m_pmdc_2);
      m_tuple3->addItem("ptmdc_1", m_ptmdc_1);
      m_tuple3->addItem("ptmdc_2", m_ptmdc_2);
      m_tuple3->addItem("pkal_1", 5, m_pkal_1);
      m_tuple3->addItem("pkal_2", 5, m_pkal_2);
      m_tuple3->addItem("ptkal_1", 5, m_ptkal_1);
      m_tuple3->addItem("ptkal_2", 5, m_ptkal_2);
      m_tuple3->addItem("p3_1", 3, m_p3_1);
      m_tuple3->addItem("p3_2", 3, m_p3_2);
      m_tuple3->addItem("theta_1", m_theta_1);
      m_tuple3->addItem("theta_2", m_theta_2);

      // PID information
      m_tuple3->addItem("probe_1", m_probe_1);
      m_tuple3->addItem("probmu_1", m_probmu_1);
      m_tuple3->addItem("probpi_1", m_probpi_1);
      m_tuple3->addItem("probk_1", m_probk_1);
      m_tuple3->addItem("probp_1", m_probp_1);
      m_tuple3->addItem("probe_2", m_probe_2);
      m_tuple3->addItem("probmu_2", m_probmu_2);
      m_tuple3->addItem("probpi_2", m_probpi_2);
      m_tuple3->addItem("probk_2", m_probk_2);
      m_tuple3->addItem("probp_2", m_probp_2);

      m_tuple3->addItem("nTofInfo_1", m_nTofInfo_1, 0, 10);
      m_tuple3->addItem("nTofInfo_2", m_nTofInfo_2, 0, 10);
      m_tuple3->addItem("fgtof_1", m_fgtof_1);
      m_tuple3->addItem("counter_1", m_counter_1);
      m_tuple3->addItem("isbarrel_1", m_isbarrel_1);
      m_tuple3->addItem("layertof_1", m_layertof_1);
      m_tuple3->addItem("iscluster_1", m_iscluster_1);
      m_tuple3->addItem("tof_1", m_tof_1);
      m_tuple3->addItem("texe_1", m_texe_1);
      m_tuple3->addItem("texmu_1", m_texmu_1);
      m_tuple3->addItem("texpi_1", m_texpi_1);
      m_tuple3->addItem("texk_1", m_texk_1);
      m_tuple3->addItem("texp_1", m_texp_1);
      m_tuple3->addItem("dte_1", m_dte_1);
      m_tuple3->addItem("dtmu_1", m_dtmu_1);
      m_tuple3->addItem("dtpi_1", m_dtpi_1);
      m_tuple3->addItem("dtk_1", m_dtk_1);
      m_tuple3->addItem("dtp_1", m_dtp_1);
      m_tuple3->addItem("fgtof_2", m_fgtof_2);
      m_tuple3->addItem("counter_2", m_counter_2);
      m_tuple3->addItem("isbarrel_2", m_isbarrel_2);
      m_tuple3->addItem("layertof_2", m_layertof_2);
      m_tuple3->addItem("iscluster_2", m_iscluster_2);
      m_tuple3->addItem("tof_2", m_tof_2);
      m_tuple3->addItem("texe_2", m_texe_2);
      m_tuple3->addItem("texmu_2", m_texmu_2);
      m_tuple3->addItem("texpi_2", m_texpi_2);
      m_tuple3->addItem("texk_2", m_texk_2);
      m_tuple3->addItem("texp_2", m_texp_2);
      m_tuple3->addItem("dte_2", m_dte_2);
      m_tuple3->addItem("dtmu_2", m_dtmu_2);
      m_tuple3->addItem("dtpi_2", m_dtpi_2);
      m_tuple3->addItem("dtk_2", m_dtk_2);
      m_tuple3->addItem("dtp_2", m_dtp_2);

      // MUC information
      m_tuple3->addItem("maxhitsinlay_1", m_maxhitsinlay_1);
      m_tuple3->addItem("numhits_1", m_numhits_1);
      m_tuple3->addItem("numlayers_1", m_numlayers_1);
      m_tuple3->addItem("depth_1", m_depth_1);
      m_tuple3->addItem("mucchi2_1", m_mucchi2_1);
      m_tuple3->addItem("maxhitsinlay_2", m_maxhitsinlay_2);
      m_tuple3->addItem("numhits_2", m_numhits_2);
      m_tuple3->addItem("numlayers_2", m_numlayers_2);
      m_tuple3->addItem("depth_2", m_depth_2);
      m_tuple3->addItem("mucchi2_2", m_mucchi2_2);

      // EMC information
      m_tuple3->addItem("evp_1", m_evp_1);
      m_tuple3->addItem("ene_1", m_ene_1);
      m_tuple3->addItem("evp_2", m_evp_2);
      m_tuple3->addItem("ene_2", m_ene_2);
      m_tuple3->addItem("nTrk_EMC", m_nTrk_EMC);

      // dedx information
      m_tuple3->addItem("chi_e_1", m_chi_e_1);
      m_tuple3->addItem("chi_mu_1", m_chi_mu_1);
      m_tuple3->addItem("chi_pi_1", m_chi_pi_1);
      m_tuple3->addItem("chi_k_1", m_chi_k_1);
      m_tuple3->addItem("chi_p_1", m_chi_p_1);
      m_tuple3->addItem("chi_e_2", m_chi_e_2);
      m_tuple3->addItem("chi_mu_2", m_chi_mu_2);
      m_tuple3->addItem("chi_pi_2", m_chi_pi_2);
      m_tuple3->addItem("chi_k_2", m_chi_k_2);
      m_tuple3->addItem("chi_p_2", m_chi_p_2);

      m_tuple3->addItem("mTrks", 5, m_mTrks);
      m_tuple3->addItem("m2Trks", 5, m_m2Trks);
      m_tuple3->addItem("pkal_tot_px", 5, m_pkal_tot_px);
      m_tuple3->addItem("pkal_tot_py", 5, m_pkal_tot_py);
      m_tuple3->addItem("pkal_tot_pz", 5, m_pkal_tot_pz);
      m_tuple3->addItem("pkal_tot_e", 5, m_pkal_tot_e);
      m_tuple3->addItem("ang_mdc_acol", m_ang_mdc_acol);
      m_tuple3->addItem("ang_mdc_acop", m_ang_mdc_acop);
      m_tuple3->addItem("the_add", m_the_add);
      m_tuple3->addItem("phi_diff", m_phi_diff);
      m_tuple3->addItem("ptem", m_PTEM);
      m_tuple3->addItem("ee_acop", m_ee_acop);
      m_tuple3->addItem("emu_acop", m_emu_acop);
      m_tuple3->addItem("mue_acop", m_mue_acop);
      m_tuple3->addItem("epi_acop", m_epi_acop);
      m_tuple3->addItem("pie_acop", m_pie_acop);
      m_tuple3->addItem("ek_acop", m_ek_acop);
      m_tuple3->addItem("ke_acop", m_ke_acop);
      m_tuple3->addItem("mumu_acop", m_mumu_acop);
      m_tuple3->addItem("mupi_acop", m_mupi_acop);
      m_tuple3->addItem("pimu_acop", m_pimu_acop);
      m_tuple3->addItem("pipi_acop", m_pipi_acop);
      m_tuple3->addItem("muk_acop", m_muk_acop);
      m_tuple3->addItem("kmu_acop", m_kmu_acop);
      m_tuple3->addItem("pik_acop", m_pik_acop);
      m_tuple3->addItem("kpi_acop", m_kpi_acop);
      m_tuple3->addItem("kk_acop", m_kk_acop);
      m_tuple3->addItem("erho_11_acop", m_erho_11_acop);
      m_tuple3->addItem("erho_12_acop", m_erho_12_acop);
      m_tuple3->addItem("murho_11_acop", m_murho_11_acop);
      m_tuple3->addItem("murho_12_acop", m_murho_12_acop);
      m_tuple3->addItem("rhorho_11_acop", m_rhorho_11_acop);
      m_tuple3->addItem("rhorho_12_acop", m_rhorho_12_acop);
      m_tuple3->addItem("ee_PTEM", m_ee_PTEM);
      m_tuple3->addItem("emu_PTEM", m_emu_PTEM);
      m_tuple3->addItem("mue_PTEM", m_mue_PTEM);
      m_tuple3->addItem("epi_PTEM", m_epi_PTEM);
      m_tuple3->addItem("pie_PTEM", m_pie_PTEM);
      m_tuple3->addItem("ek_PTEM", m_ek_PTEM);
      m_tuple3->addItem("ke_PTEM", m_ke_PTEM);
      m_tuple3->addItem("mumu_PTEM", m_mumu_PTEM);
      m_tuple3->addItem("mupi_PTEM", m_mupi_PTEM);
      m_tuple3->addItem("pimu_PTEM", m_pimu_PTEM);
      m_tuple3->addItem("pipi_PTEM", m_pipi_PTEM);
      m_tuple3->addItem("muk_PTEM", m_muk_PTEM);

      m_tuple3->addItem("kmu_PTEM", m_kmu_PTEM);
      m_tuple3->addItem("pik_PTEM", m_pik_PTEM);
      m_tuple3->addItem("kpi_PTEM", m_kpi_PTEM);
      m_tuple3->addItem("kk_PTEM", m_kk_PTEM);
      m_tuple3->addItem("dlt_mpi0", m_dlt_mpi0);
      m_tuple3->addItem("eg_11", m_eg_11);
      m_tuple3->addItem("eg_12", m_eg_12);
      m_tuple3->addItem("eg_21", m_eg_21);
      m_tuple3->addItem("eg_22", m_eg_22);
      m_tuple3->addItem("mpi0_1", m_mpi0_1);
      m_tuple3->addItem("mpi0_2", m_mpi0_2);
      m_tuple3->addItem("mrho_11", m_mrho_11);
      m_tuple3->addItem("mrho_12", m_mrho_12);
      m_tuple3->addItem("mrho_21", m_mrho_21);
      m_tuple3->addItem("mrho_22", m_mrho_22);
      m_tuple3->addItem("prho_11", m_prho_11);
      m_tuple3->addItem("prho_12", m_prho_12);
      m_tuple3->addItem("prho_21", m_prho_21);
      m_tuple3->addItem("prho_22", m_prho_22);
      m_tuple3->addItem("theta_pi_11", m_theta_pi_11);
      m_tuple3->addItem("theta_pi_12", m_theta_pi_12);
      m_tuple3->addItem("theta_pi_21", m_theta_pi_21);
      m_tuple3->addItem("theta_pi_22", m_theta_pi_22);
      m_tuple3->addItem("ppi_cms_11", m_ppi_cms_11);
      m_tuple3->addItem("ppi_cms_12", m_ppi_cms_12);
      m_tuple3->addItem("ppi_cms_21", m_ppi_cms_21);
      m_tuple3->addItem("ppi_cms_22", m_ppi_cms_22);
      m_tuple3->addItem("theta_m_11", m_theta_m_11);
      m_tuple3->addItem("theta_m_12", m_theta_m_12);
      m_tuple3->addItem("theta_m_21", m_theta_m_21);
      m_tuple3->addItem("theta_m_22", m_theta_m_22);
      //*****************

      // angle between two charged tracks
      m_tuple3->addItem("ee_angle", m_ee_angle);
      m_tuple3->addItem("emu_angle", m_emu_angle);
      m_tuple3->addItem("mue_angle", m_mue_angle);
      m_tuple3->addItem("epi_angle", m_epi_angle);
      m_tuple3->addItem("pie_angle", m_pie_angle);
      m_tuple3->addItem("ek_angle", m_ek_angle);
      m_tuple3->addItem("ke_angle", m_ke_angle);
      m_tuple3->addItem("mumu_angle", m_mumu_angle);
      m_tuple3->addItem("mupi_angle", m_mupi_angle);
      m_tuple3->addItem("pimu_angle", m_pimu_angle);
      m_tuple3->addItem("pipi_angle", m_pipi_angle);
      m_tuple3->addItem("muk_angle", m_muk_angle);
      m_tuple3->addItem("kmu_angle", m_kmu_angle);
      m_tuple3->addItem("pik_angle", m_pik_angle);
      m_tuple3->addItem("kpi_angle", m_kpi_angle);
      m_tuple3->addItem("kk_angle", m_kk_angle);
      m_tuple3->addItem("erho_11_angle", m_erho_11_angle);
      m_tuple3->addItem("erho_12_angle", m_erho_12_angle);
      m_tuple3->addItem("murho_11_angle", m_murho_11_angle);
      m_tuple3->addItem("murho_12_angle", m_murho_12_angle);
      m_tuple3->addItem("rhorho_11_angle", m_rhorho_11_angle);
      m_tuple3->addItem("rhorho_12_angle", m_rhorho_12_angle);

      m_tuple3->addItem("miss_m2_ee", m_miss_m2_ee);
      m_tuple3->addItem("miss_m2_emu", m_miss_m2_emu);
      m_tuple3->addItem("miss_m2_mue", m_miss_m2_mue);
      m_tuple3->addItem("miss_m2_epi", m_miss_m2_epi);
      m_tuple3->addItem("miss_m2_pie", m_miss_m2_pie);
      m_tuple3->addItem("miss_m2_ek", m_miss_m2_ek);
      m_tuple3->addItem("miss_m2_ke", m_miss_m2_ke);
      m_tuple3->addItem("miss_m2_mumu", m_miss_m2_mumu);
      m_tuple3->addItem("miss_m2_mupi", m_miss_m2_mupi);
      m_tuple3->addItem("miss_m2_pimu", m_miss_m2_pimu);
      m_tuple3->addItem("miss_m2_pipi", m_miss_m2_pipi);
      m_tuple3->addItem("miss_m2_muk", m_miss_m2_muk);
      m_tuple3->addItem("miss_m2_kmu", m_miss_m2_kmu);
      m_tuple3->addItem("miss_m2_pik", m_miss_m2_pik);
      m_tuple3->addItem("miss_m2_kpi", m_miss_m2_kpi);
      m_tuple3->addItem("miss_m2_kk", m_miss_m2_kk);
      m_tuple3->addItem("miss_m2_erho_11", m_miss_m2_erho_11);
      m_tuple3->addItem("miss_m2_erho_12", m_miss_m2_erho_12);
      m_tuple3->addItem("miss_m2_murho_11", m_miss_m2_murho_11);
      m_tuple3->addItem("miss_m2_murho_12", m_miss_m2_murho_12);
      m_tuple3->addItem("miss_m2_rhorho_11", m_miss_m2_rhorho_11);
      m_tuple3->addItem("miss_m2_rhorho_12", m_miss_m2_rhorho_12);

      m_tuple3->addItem("nsig_shower", m_nsig_shower, 0, 100);
      m_tuple3->addIndexedItem("shower_info", m_nsig_shower, 19, m_shower_info);
      m_tuple3->addItem("tot_egam", m_esum);
      m_tuple3->addItem("max_egam", m_egam_max);
      m_tuple3->addItem("max_gam_id", m_id_gam_max);
      m_tuple3->addItem("gamcon", 26, m_gamcon);
      m_tuple3->addItem("info_elp", 26, m_elp_info);
      m_tuple3->addItem("info_mum", 26, m_mum_info);
      m_tuple3->addItem("info_eta_gamee", 16, m_eta_gamee);
      m_tuple3->addItem("info_pi0_gamee", 16, m_pi0_gamee);
      m_tuple3->addItem("prob_elp", 6, m_elp_prob);
      m_tuple3->addItem("prob_mum", 6, m_mum_prob);

      m_tuple3->addItem("nelpFSR", m_nleppFSR, 0, 100);
      m_tuple3->addItem("elpFSRArray", m_nleppFSR, m_eleppFSRArray);
      m_tuple3->addItem("dtheelpFSRArray", m_nleppFSR, m_dtheleppFSRArray);

      m_tuple3->addItem("info_gam_fromelp", 8, m_gam_fromelp);
      m_tuple3->addItem("info_elp_corfsr", 8, m_elp_corfsr);
      m_tuple3->addItem("info_elp_mum", 3, m_info_elp_mum);
      m_tuple3->addItem("info_p4_emu", 3, 8, m_info_emu);

      // For topology
      m_tuple3->addItem("nTrackMC", m_nTrack, 0, 100);
      m_tuple3->addIndexedItem("trackIDMC", m_nTrack, ma_trackID);
      m_tuple3->addIndexedItem("trackIndexMC", m_nTrack, ma_trackIndex);
      m_tuple3->addIndexedItem("motherIDMC", m_nTrack, ma_motherID);
      m_tuple3->addIndexedItem("motherIndexMC", m_nTrack, ma_motherIndex);
      m_tuple3->addIndexedItem("fromGenerato3MC", m_nTrack, ma_fromGenerator);
      m_tuple3->addIndexedItem("primaryParticleMC", m_nTrack,
                               ma_primaryParticle);
      m_tuple3->addIndexedItem("pxMC", m_nTrack, ma_px);
      m_tuple3->addIndexedItem("pyMC", m_nTrack, ma_py);
      m_tuple3->addIndexedItem("pzMC", m_nTrack, ma_pz);
      m_tuple3->addIndexedItem("eMC", m_nTrack, ma_e);
      m_tuple3->addItem("numGam", m_nGam, 0, 100);

      // Kalman1CFit for reconstruct pi0  2024.04.10
      m_tuple3->addItem("numpi0", m_npi0, 0, 50);
      m_tuple3->addIndexedItem("m_pi0before", m_npi0, m_m_pi0before);
      m_tuple3->addIndexedItem("m_pi0fit", m_npi0, m_m_pi0fit);
      m_tuple3->addIndexedItem("chi_pi0fit", m_npi0, m_chi_pi0fit);
      m_tuple3->addIndexedItem("i1gam", m_npi0, m_i1gam_pi0);
      m_tuple3->addIndexedItem("i2gam", m_npi0, m_i2gam_pi0);

      // track and gamma angle /invariant mass
      m_tuple3->addIndexedItem("gam_eraw", m_nGam, m_eraw);
      m_tuple3->addIndexedItem("gam_phi", m_nGam, m_phi);
      m_tuple3->addIndexedItem("gam_the", m_nGam, m_the);

      m_tuple3->addIndexedItem("dang_pos_e", m_nGam, m_dang_pos_e);
      m_tuple3->addIndexedItem("dthe_pos_e", m_nGam, m_dthe_pos_e);
      m_tuple3->addIndexedItem("dphi_pos_e", m_nGam, m_dphi_pos_e);
      m_tuple3->addIndexedItem("dang_neg_e", m_nGam, m_dang_neg_e);
      m_tuple3->addIndexedItem("dthe_neg_e", m_nGam, m_dthe_neg_e);
      m_tuple3->addIndexedItem("dphi_neg_e", m_nGam, m_dphi_neg_e);
      m_tuple3->addIndexedItem("dang_pos_mu", m_nGam, m_dang_pos_mu);
      m_tuple3->addIndexedItem("dthe_pos_mu", m_nGam, m_dthe_pos_mu);
      m_tuple3->addIndexedItem("dphi_pos_mu", m_nGam, m_dphi_pos_mu);
      m_tuple3->addIndexedItem("dang_neg_mu", m_nGam, m_dang_neg_mu);
      m_tuple3->addIndexedItem("dthe_neg_mu", m_nGam, m_dthe_neg_mu);
      m_tuple3->addIndexedItem("dphi_neg_mu", m_nGam, m_dphi_neg_mu);
      m_tuple3->addIndexedItem("dang_pos_pi", m_nGam, m_dang_pos_pi);
      m_tuple3->addIndexedItem("dthe_pos_pi", m_nGam, m_dthe_pos_pi);
      m_tuple3->addIndexedItem("dphi_pos_pi", m_nGam, m_dphi_pos_pi);
      m_tuple3->addIndexedItem("dang_neg_pi", m_nGam, m_dang_neg_pi);
      m_tuple3->addIndexedItem("dthe_neg_pi", m_nGam, m_dthe_neg_pi);
      m_tuple3->addIndexedItem("dphi_neg_pi", m_nGam, m_dphi_neg_pi);
      m_tuple3->addIndexedItem("dang_pos_k", m_nGam, m_dang_pos_k);
      m_tuple3->addIndexedItem("dthe_pos_k", m_nGam, m_dthe_pos_k);
      m_tuple3->addIndexedItem("dphi_pos_k", m_nGam, m_dphi_pos_k);
      m_tuple3->addIndexedItem("dang_neg_k", m_nGam, m_dang_neg_k);
      m_tuple3->addIndexedItem("dthe_neg_k", m_nGam, m_dthe_neg_k);
      m_tuple3->addIndexedItem("dphi_neg_k", m_nGam, m_dphi_neg_k);

      m_tuple3->addItem("pos_chg_px_e", m_px_pos_e);
      m_tuple3->addItem("pos_chg_py_e", m_py_pos_e);
      m_tuple3->addItem("pos_chg_pz_e", m_pz_pos_e);
      m_tuple3->addItem("pos_chg_px_mu", m_px_pos_mu);
      m_tuple3->addItem("pos_chg_py_mu", m_py_pos_mu);
      m_tuple3->addItem("pos_chg_pz_mu", m_pz_pos_mu);
      m_tuple3->addItem("pos_chg_px_pi", m_px_pos_pi);
      m_tuple3->addItem("pos_chg_py_pi", m_py_pos_pi);
      m_tuple3->addItem("pos_chg_pz_pi", m_pz_pos_pi);
      m_tuple3->addItem("pos_chg_px_k", m_px_pos_k);
      m_tuple3->addItem("pos_chg_py_k", m_py_pos_k);
      m_tuple3->addItem("pos_chg_pz_k", m_pz_pos_k);
      m_tuple3->addItem("pos_chg_p", m_p_pos);
      m_tuple3->addItem("pos_chg_pt", m_pt_pos);

      m_tuple3->addItem("neg_chg_px_e", m_px_neg_e);
      m_tuple3->addItem("neg_chg_py_e", m_py_neg_e);
      m_tuple3->addItem("neg_chg_pz_e", m_pz_neg_e);
      m_tuple3->addItem("neg_chg_px_mu", m_px_neg_mu);
      m_tuple3->addItem("neg_chg_py_mu", m_py_neg_mu);
      m_tuple3->addItem("neg_chg_pz_mu", m_pz_neg_mu);
      m_tuple3->addItem("neg_chg_px_pi", m_px_neg_pi);
      m_tuple3->addItem("neg_chg_py_pi", m_py_neg_pi);
      m_tuple3->addItem("neg_chg_pz_pi", m_pz_neg_pi);
      m_tuple3->addItem("neg_chg_px_k", m_px_neg_k);
      m_tuple3->addItem("neg_chg_py_k", m_py_neg_k);
      m_tuple3->addItem("neg_chg_pz_k", m_pz_neg_k);
      m_tuple3->addItem("neg_chg_p", m_p_neg);
      m_tuple3->addItem("neg_chg_pt", m_pt_neg);

      m_tuple3->addItem("cang_ee", m_cang_ee);
      m_tuple3->addItem("cthe_ee", m_cthe_ee);
      m_tuple3->addItem("cphi_ee", m_cphi_ee);
      m_tuple3->addItem("costhe_ee", m_costheta_ee);
      m_tuple3->addItem("cang_emu", m_cang_emu);
      m_tuple3->addItem("cthe_emu", m_cthe_emu);
      m_tuple3->addItem("cphi_emu", m_cphi_emu);
      m_tuple3->addItem("costhe_emu", m_costheta_emu);
      m_tuple3->addItem("cang_epi", m_cang_epi);
      m_tuple3->addItem("cthe_epi", m_cthe_epi);
      m_tuple3->addItem("cphi_epi", m_cphi_epi);
      m_tuple3->addItem("costhe_epi", m_costheta_epi);
      m_tuple3->addItem("cang_ek", m_cang_ek);
      m_tuple3->addItem("cthe_ek", m_cthe_ek);
      m_tuple3->addItem("cphi_ek", m_cphi_ek);
      m_tuple3->addItem("costhe_ek", m_costheta_ek);
      m_tuple3->addItem("cang_mue", m_cang_mue);
      m_tuple3->addItem("cthe_mue", m_cthe_mue);
      m_tuple3->addItem("cphi_mue", m_cphi_mue);
      m_tuple3->addItem("costhe_mue", m_costheta_mue);
      m_tuple3->addItem("cang_mumu", m_cang_mumu);
      m_tuple3->addItem("cthe_mumu", m_cthe_mumu);
      m_tuple3->addItem("cphi_mumu", m_cphi_mumu);
      m_tuple3->addItem("costhe_mumu", m_costheta_mumu);
      m_tuple3->addItem("cang_mupi", m_cang_mupi);
      m_tuple3->addItem("cthe_mupi", m_cthe_mupi);
      m_tuple3->addItem("cphi_mupi", m_cphi_mupi);
      m_tuple3->addItem("costhe_mupi", m_costheta_mupi);
      m_tuple3->addItem("cang_muk", m_cang_muk);
      m_tuple3->addItem("cthe_muk", m_cthe_muk);
      m_tuple3->addItem("cphi_muk", m_cphi_muk);
      m_tuple3->addItem("costhe_muk", m_costheta_muk);
      m_tuple3->addItem("cang_pie", m_cang_pie);
      m_tuple3->addItem("cthe_pie", m_cthe_pie);
      m_tuple3->addItem("cphi_pie", m_cphi_pie);
      m_tuple3->addItem("costhe_pie", m_costheta_pie);
      m_tuple3->addItem("cang_pimu", m_cang_pimu);
      m_tuple3->addItem("cthe_pimu", m_cthe_pimu);
      m_tuple3->addItem("cphi_pimu", m_cphi_pimu);
      m_tuple3->addItem("costhe_pimu", m_costheta_pimu);
      m_tuple3->addItem("cang_pipi", m_cang_pipi);
      m_tuple3->addItem("cthe_pipi", m_cthe_pipi);
      m_tuple3->addItem("cphi_pipi", m_cphi_pipi);
      m_tuple3->addItem("costhe_pipi", m_costheta_pipi);
      m_tuple3->addItem("cang_pik", m_cang_pik);
      m_tuple3->addItem("cthe_pik", m_cthe_pik);
      m_tuple3->addItem("cphi_pik", m_cphi_pik);
      m_tuple3->addItem("costhe_pik", m_costheta_pik);
      m_tuple3->addItem("cang_ke", m_cang_ke);
      m_tuple3->addItem("cthe_ke", m_cthe_ke);
      m_tuple3->addItem("cphi_ke", m_cphi_ke);
      m_tuple3->addItem("costhe_ke", m_costheta_ke);
      m_tuple3->addItem("cang_kmu", m_cang_kmu);
      m_tuple3->addItem("cthe_kmu", m_cthe_kmu);
      m_tuple3->addItem("cphi_kmu", m_cphi_kmu);
      m_tuple3->addItem("costhe_kmu", m_costheta_kmu);
      m_tuple3->addItem("cang_kpi", m_cang_kpi);
      m_tuple3->addItem("cthe_kpi", m_cthe_kpi);
      m_tuple3->addItem("cphi_kpi", m_cphi_kpi);
      m_tuple3->addItem("costhe_kpi", m_costheta_kpi);
      m_tuple3->addItem("cang_kk", m_cang_kk);
      m_tuple3->addItem("cthe_kk", m_cthe_kk);
      m_tuple3->addItem("cphi_kk", m_cphi_kk);
      m_tuple3->addItem("costhe_kk", m_costheta_kk);

      m_tuple3->addIndexedItem("egam_invam_pos", m_nGam, m_m_egam_pos);
      m_tuple3->addIndexedItem("mugam_invam_pos", m_nGam, m_m_mugam_pos);
      m_tuple3->addIndexedItem("pigam_invam_pos", m_nGam, m_m_pigam_pos);
      m_tuple3->addIndexedItem("kgam_invam_pos", m_nGam, m_m_kgam_pos);
      m_tuple3->addIndexedItem("egam_invam_neg", m_nGam, m_m_egam_neg);
      m_tuple3->addIndexedItem("mugam_invam_neg", m_nGam, m_m_mugam_neg);
      m_tuple3->addIndexedItem("pigam_invam_neg", m_nGam, m_m_pigam_neg);
      m_tuple3->addIndexedItem("kgam_invam_neg", m_nGam, m_m_kgam_neg);

      m_tuple3->addItem("nRecEmcHits", m_nRecEmcHits, 0, 1000);
      m_tuple3->addIndexedItem("emcHit_id_theta", m_nRecEmcHits,
                               m_emc_id_theta);
      m_tuple3->addIndexedItem("emcHit_id_phi", m_nRecEmcHits, m_emc_id_phi);
      m_tuple3->addIndexedItem("emcHit_energy", m_nRecEmcHits, m_emc_energy);
      m_tuple3->addIndexedItem("emcHit_bc", m_nRecEmcHits, m_emc_bc);
      m_tuple3->addIndexedItem("emcHit_tdc", m_nRecEmcHits, m_emc_tdc);
      m_tuple3->addIndexedItem("emcHit_x", m_nRecEmcHits, m_emc_x);
      m_tuple3->addIndexedItem("emcHit_y", m_nRecEmcHits, m_emc_y);
      m_tuple3->addIndexedItem("emcHit_z", m_nRecEmcHits, m_emc_z);
      m_tuple3->addIndexedItem("emcHit_pos_theta", m_nRecEmcHits,
                               m_emc_pos_theta);
      m_tuple3->addIndexedItem("emcHit_pos_phi", m_nRecEmcHits, m_emc_pos_phi);

    } else {
      log << MSG::ERROR << "    Cannot book N-tuple:" << long(m_tuple3)
          << endmsg;
      return StatusCode::FAILURE;
    }
  }

  log << MSG::INFO << "successfully return from initialize()" << endmsg;
  return StatusCode::SUCCESS;
}

//***************************Execution Event
// Filtering******************************
StatusCode tautauAlg::execute() {
  MsgStream log(msgSvc(), name());

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),
                                               "/Event/EventHeader");
  if (!eventHeader) {
    log << MSG::FATAL << "Could not find Event Header" << endreq;

    return (StatusCode::FAILURE);
  }
  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),
                                        "/Event/EvtRec/EvtRecEvent");
  if (!evtRecEvent) {
    log << MSG::FATAL << "Could not find EvtRecEvent" << endreq;
    return StatusCode::FAILURE;
  }
  SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),
                                            "/Event/EvtRec/EvtRecTrackCol");
  if (!evtRecTrkCol) {
    log << MSG::FATAL << "Could not find EvtRecTrackCol" << endreq;
    return StatusCode::FAILURE;
  }

  EvtRecTrackIterator track_begin = evtRecTrkCol->begin();
  EvtRecTrackIterator charged_begin = track_begin;
  EvtRecTrackIterator charged_end = track_begin + evtRecEvent->totalCharged();

  EvtRecTrackIterator neutral_begin = track_begin + evtRecEvent->totalCharged();
  EvtRecTrackIterator neutral_end = track_begin + evtRecEvent->totalTracks();

  int runNo = eventHeader->runNumber();
  int event = eventHeader->eventNumber();

  m_cutflow->Fill(0);
  Ncut[0]++;

  m_run = runNo;
  m_event = event;
  m_rec_run = runNo;
  m_rec_event = event;
  m_rec_nchrg = evtRecEvent->totalCharged();
  m_rec_nneu = evtRecEvent->totalNeutral();

  if (0 == (m_nEvt % m_nEvtDisp))
    std::cout << "Run num : " << runNo << ",  Event " << m_nEvt
              << ",   Event number online " << event << endl;
  m_nEvt++;

  // 2020.03.27 for topology
  ////////////////////////////////////////////////////////////////////////////////////wu
  /// begin
  //
  if (eventHeader->runNumber() < 0 && m_testMC == 1) {
    int m_numParticle = 0;
    int gamfromTau = 0;
    double pow_tau = 3.1572315;
    double pow_r_2 = 0.0061699104;
    int nefrommu = 0;
    int nefromtau = 0;
    m_pi_e = 0.0;
    m_gam_e = 0.0;
    m_tau_e = 0.0;
    m_nu_e = 0.0;
    // MC information
    SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
                                                     "/Event/MC/McParticleCol");
    if (!mcParticleCol) {
      std::cout << "Could not retrieve McParticelCol" << std::endl;
      return StatusCode::FAILURE;
    }
    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();

    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      HepLorentzVector lvMC = (*iter_mc)->initialFourMomentum();
      // what do you mean by ma??
      ma_trackID[m_numParticle] = (*iter_mc)->particleProperty();
      ma_motherID[m_numParticle] = ((*iter_mc)->mother()).particleProperty();
      ma_trackIndex[m_numParticle] = (*iter_mc)->trackIndex();
      ma_motherIndex[m_numParticle] = ((*iter_mc)->mother()).trackIndex();
      ma_fromGenerator[m_numParticle] = (*iter_mc)->decayFromGenerator();
      ma_primaryParticle[m_numParticle] = (*iter_mc)->primaryParticle();
      ma_px[m_numParticle] = lvMC.px();
      ma_py[m_numParticle] = lvMC.py();
      ma_pz[m_numParticle] = lvMC.pz();
      ma_e[m_numParticle] = lvMC.e();

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == -211 &&
          ((*iter_mc)->mother()).particleProperty() == 15 &&
          lvMC.e() > m_pi_e) {
        m_pi_px = lvMC.px();
        m_pi_py = lvMC.py();
        m_pi_pz = lvMC.pz();
        m_pi_e = lvMC.e();
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == 22 && lvMC.e() > m_gam_e) {
        m_gam_px = lvMC.px();
        m_gam_py = lvMC.py();
        m_gam_pz = lvMC.pz();
        m_gam_e = lvMC.e();
        m_gam_motherIndex[gamfromTau] = ((*iter_mc)->mother()).trackIndex();
        gamfromTau += 1;
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == 16 &&
          ((*iter_mc)->mother()).particleProperty() == 15 &&
          lvMC.e() > m_nu_e) {
        m_nu_px = lvMC.px();
        m_nu_py = lvMC.py();
        m_nu_pz = lvMC.pz();
        m_nu_e = lvMC.e();
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == 15 && lvMC.e() > m_tau_e) {
        m_tau_px = lvMC.px();
        m_tau_py = lvMC.py();
        m_tau_pz = lvMC.pz();
        m_tau_e = lvMC.e();
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == -11 &&
          ((*iter_mc)->mother()).particleProperty() == -13) {
        nefrommu += 1;
      }

      if ((*iter_mc)->decayFromGenerator() &&
          (*iter_mc)->particleProperty() == -11 &&
          ((*iter_mc)->mother()).particleProperty() == -15) {
        nefromtau += 1;
      }

      m_numParticle += 1;
    } // end of loop of all particles

    HepLorentzVector p4Pi, p4Nu, p4Gam, p4Tau;
    setPxPyPzE(p4Pi, m_pi_px, m_pi_py, m_pi_pz, m_pi_e);
    setPxPyPzE(p4Nu, m_nu_px, m_nu_py, m_nu_pz, m_nu_e);
    setPxPyPzE(p4Gam, m_gam_px, m_gam_py, m_gam_pz, m_gam_e);
    setPxPyPzE(p4Tau, m_tau_px, m_tau_py, m_tau_pz, m_tau_e);

    m_x = 2.0 * p4Tau.dot(p4Gam) / pow_tau;
    m_y = 2.0 * p4Tau.dot(p4Pi) / pow_tau;
    m_z = m_x + m_y - 1.0;
    m_t = (p4Pi + p4Gam).m2();
    m_s = pow_tau * m_z; // m_s is the same as m_t

    m_within_kine_bounds = false;
    if (m_x >= 0 && m_x <= (1.0 - pow_r_2) &&
        m_y >= (1.0 - m_x + pow_r_2 / (1.0 - m_x)) && m_y <= (1 + pow_r_2)) {
      m_within_kine_bounds = true;
    }
    m_ngamfromTau = gamfromTau;
    m_nefrommu = nefrommu;
    m_nefromtau = nefromtau;
    // m_nTrack = m_numParticle;
  } // if this is MC and testMC is requested
    ///////////////////////////////////////////////////////////////////////////////////

  //***************************Access Truth
  // information******************************
  // 2020.3.9, add the MC truth information
  if (eventHeader->runNumber() < 0) {
    mctruth();
    m_tuple1->write();
    m_tuple2->write();
  }

  if (m_writeGenOnly) {
    return StatusCode::SUCCESS;
  }

  // reconstruct 2 good charged tracks
  Hep3Vector xorigin(0, 0, 0);
  VertexParameter m_privxpar;
  MyInitIP myInitIP;
  myInitIP.InitIP(xorigin, m_privxpar);

  Vint iGood;
  iGood.clear();
  Vint icharge;
  icharge.clear();
  int charge = 0;
  for (int i = 0; i < evtRecEvent->totalCharged(); i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (isGoodTrack(*itTrk, xorigin)) { // 2020.3.9, not use mdcKalTrack
      iGood.push_back(i);
      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
      charge += mdcTrk->charge();
      icharge.push_back(mdcTrk->charge());
    }
  }
  int nGood = iGood.size();
  if (!(nGood == 2 && charge == 0))
    return StatusCode::SUCCESS;
  m_cutflow->Fill(1);
  Ncut[1]++;

  // make sure the first track is positive, the second one is negetive;
  if (icharge[0] == -1) {
    int temidx = iGood[0];
    iGood[0] = iGood[1];
    iGood[1] = temidx;
  }

  //***************************Select Good Photons******************************
  Vint iGam;
  Vint shelf;
  iGam.clear();
  Vp4 pGam;
  pGam.clear();

  int igamma = 0;
  double esum = 0.;
  // The energy and id of gam with maximum energy
  double tmp_e_gam_max = 0.; // 2020.04.13 modify from -100 to 0 !!
  int tmp_gam_max_id = -100;

  for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks();
       i++) {
    if (i >= evtRecTrkCol->size()) {
      break;
    }
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + i;
    if (!isGoodShower(*itTrk)) {
      continue;
    }
    RecEmcShower *aShower = NULL;
    aShower = (*itTrk)->emcShower();
    esum += aShower->energy();

    double tmp_e_gam = 0.;
    tmp_e_gam = aShower->energy();
    if (tmp_e_gam > tmp_e_gam_max) {
      tmp_e_gam_max = tmp_e_gam;
      tmp_gam_max_id = igamma;
    }
    MyInfoShower(aShower, m_shower_info, igamma);
    igamma++;
    iGam.push_back(i);
  }

  int nGam = iGam.size();
  //  std::cout<< "nGam = " << nGam << std::endl;
  if (nGam > 100 || nGam < m_gamNumCut) {
    return StatusCode::SUCCESS;
  }
  m_cutflow->Fill(2);
  Ncut[2]++;
  m_nGam = nGam;

  int count = 0;
  m_nsig_shower = igamma;
  m_esum = esum;
  m_egam_max = tmp_e_gam_max;
  m_id_gam_max = tmp_gam_max_id;
  // 2020.03.16, save gamma track
  for (int i = 0; i < nGam; i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGam[i];
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector ptrk;
    ptrk.setPx(eraw * sin(the) * cos(phi));
    ptrk.setPy(eraw * sin(the) * sin(phi));
    ptrk.setPz(eraw * cos(the));
    ptrk.setE(eraw);
    // ptrk = ptrk.boost(-0.011,0,0);// boost to cms
    pGam.push_back(ptrk);
    shelf.push_back(count);
    count += 1;
  }

  int eventFlag = -1;
  // 2 or 3 good photons
  if (nGam > 1 && nGam < 4) {
    eventFlag = 1;
  }
  // 4 or 5 good photons
  if (nGam > 3 && nGam < 6) {
    eventFlag = 2;
  }

  //------------------ Reconstruct pi0 using Kalman1Cfit-------------------
  int tmp_npi0 = 0;
  int onenumber = 0;
  int twonumber = 0;
  KalmanKinematicFit *kmfit1c = KalmanKinematicFit::instance();
  while (shelf.size() > 1) {
    double final_pi0chi2 = 9999;
    bool find_it = false;
    for (int i = 0; i < shelf.size(); i++) {
      for (int j = i + 1; j < shelf.size(); j++) {
        double mpi0 = (pGam[shelf[i]] + pGam[shelf[j]]).m();
        HepLorentzVector ppi0fit;
        RecEmcShower *g1Trk =
            (*(evtRecTrkCol->begin() + iGam[shelf[i]]))->emcShower();
        RecEmcShower *g2Trk =
            (*(evtRecTrkCol->begin() + iGam[shelf[j]]))->emcShower();
        kmfit1c->init();
        kmfit1c->setBeamPosition(xorigin);
        kmfit1c->AddTrack(0, 0.0, g1Trk);
        kmfit1c->AddTrack(1, 0.0, g2Trk);
        kmfit1c->AddResonance(0, 0.1349768, 0, 1);
        kmfit1c->setChisqCut(9999, 0.05);
        bool oksq = kmfit1c->Fit();
        double pi0chi2 = kmfit1c->chisq();
        if (!oksq || pi0chi2 > m_pi0Chi2Cut)
          continue;
        if (pi0chi2 < final_pi0chi2) {
          final_pi0chi2 = pi0chi2;
          onenumber = i;
          twonumber = j;
          find_it = true;
          ppi0fit = kmfit1c->pfit(0) + kmfit1c->pfit(1);
          m_m_pi0before[tmp_npi0] = mpi0;
          m_m_pi0fit[tmp_npi0] = ppi0fit.m();
          m_chi_pi0fit[tmp_npi0] = pi0chi2;
          m_i1gam_pi0[tmp_npi0] = iGam[shelf[i]];
          m_i2gam_pi0[tmp_npi0] = iGam[shelf[j]];
        }
      }
    }
    tmp_npi0 += 1;
    if (shelf.size() == 2 || !find_it)
      break;
    shelf.erase(shelf.begin() + onenumber);
    shelf.erase(shelf.begin() + twonumber - 1);
  }
  m_npi0 = tmp_npi0;

  if (tmp_npi0 > m_pi0NumCut) {
    return StatusCode::SUCCESS;
  }
  m_cutflow->Fill(3);
  Ncut[3]++;

  //***************************Save Information and
  //PID******************************
  // 2020.03.14-15, save information of charged tracks
  int npi_pid = 0;
  int nele_pid = 0;
  int nmu_pid = 0;
  double theta_1;
  double theta_2;
  double phi_1;
  double phi_2;
  HepLorentzVector P4_1[5];
  HepLorentzVector P4_2[5];
  HepLorentzVector cms_P4_1[5];
  HepLorentzVector cms_P4_2[5];
  Hep3Vector P3_1;
  Hep3Vector P3_2;
  RecMdcKalTrack::PidType pidtype[5] = {
    RecMdcKalTrack::electron, RecMdcKalTrack::muon,  RecMdcKalTrack::pion,
    RecMdcKalTrack::kaon,     RecMdcKalTrack::proton
  };
  Vint itrkWithValidEmcShower;
  for (int i = 0; i < nGood; i++) {
    EvtRecTrackIterator itTrk = evtRecTrkCol->begin() + iGood[i];
    if (!(*itTrk)->isMdcTrackValid()) {
      throw std::runtime_error("This should not happen. The selected good "
                               "charged tracks are already required to have "
                               "valid MdcTrack.");
    }
    RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
    if ((*itTrk)->isEmcShowerValid()) {
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      double pctrk = mdcTrk->p();
      double thetamdc = mdcTrk->theta();
      double eemc = emcTrk->energy();
      // The E/p
      double evp = eemc / pctrk;
      itrkWithValidEmcShower.push_back(iGood[i]);
      if (i < 1) {
        // Mdc information for the first track
        RecMdcKalTrack *mdckaltrk = (*itTrk)->mdcKalTrack();
        for (int j = 0; j < 5; j++) {
          RecMdcKalTrack::setPidType(pidtype[j]);
          HepLorentzVector ptrk;
          ptrk.setPx(mdckaltrk->px());
          ptrk.setPy(mdckaltrk->py());
          ptrk.setPz(mdckaltrk->pz());
          double p3 = ptrk.mag();
          ptrk.setE(sqrt(p3 * p3 + xmass[j] * xmass[j]));
          P4_1[j] = ptrk;
        }

        m_charge_1 = mdckaltrk->charge();
        theta_1 = mdcTrk->theta();
        phi_1 = mdcTrk->phi();

        m_pmdc_1 = mdcTrk->p();
        m_ptmdc_1 = mdcTrk->pxy();

        for (int j = 0; j < 5; j++) {
          m_pkal_1[j] = P4_1[j].rho();
          m_ptkal_1[j] = P4_1[j].perp();
        }

        m_theta_1 = thetamdc;
        m_p3_1[0] = mdcTrk->px();
        m_p3_1[1] = mdcTrk->py();
        m_p3_1[2] = mdcTrk->pz();
        P3_1.set(mdcTrk->px(), mdcTrk->py(), mdcTrk->pz());

        // PID information for the first track
        ParticleID *pid = ParticleID::instance();
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut(4);
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() |
                       pid->useEmc());
        pid->identify(pid->onlyProton() | pid->onlyKaon() | pid->onlyPion() |
                      pid->onlyElectron() | pid->onlyMuon());
        pid->calculate();
        if (pid->IsPidInfoValid()) {
          m_probe_1 = pid->probElectron();
          m_probmu_1 = pid->probMuon();
          m_probpi_1 = pid->probPion();
          m_probk_1 = pid->probKaon();
          m_probp_1 = pid->probProton();
        } else {
          m_probe_1 = -1.0;
          m_probmu_1 = -1.0;
          m_probpi_1 = -1.0;
          m_probk_1 = -1.0;
          m_probp_1 = -1.0;
        }

        double tmp_probe_1 = m_probe_1;
        double tmp_probmu_1 = m_probmu_1;
        double tmp_probpi_1 = m_probpi_1;
        if (tmp_probpi_1 > 0.001 && tmp_probpi_1 > tmp_probe_1 &&
            tmp_probpi_1 > tmp_probmu_1) {
          npi_pid++;
        }
        if (tmp_probe_1 > 0.001 && tmp_probe_1 > tmp_probpi_1 &&
            tmp_probe_1 > tmp_probmu_1) {
          nele_pid++;
        }
        if (tmp_probmu_1 > 0.001 && tmp_probmu_1 > tmp_probpi_1 &&
            tmp_probmu_1 > tmp_probe_1) {
          nmu_pid++;
        }

        // TOF information for the first track
        if ((*itTrk)->isTofTrackValid()) {
          SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
          SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
          int nTofInfo_1 = 0;
          //@NOTE: The last value in tofTrkCol is the average value.
          for (; iter_tof != tofTrkCol.end(); iter_tof++) {
            TofHitStatus *status = new TofHitStatus;
            status->setStatus((*iter_tof)->status());
            m_fgtof_1 = 1;
            m_counter_1 = status->is_counter();
            m_isbarrel_1 = status->is_barrel();
            m_layertof_1 = status->layer();
            m_iscluster_1 = status->is_cluster();
            m_tof_1 = (*iter_tof)->tof();
            m_texe_1 = (*iter_tof)->texpElectron();
            m_texmu_1 = (*iter_tof)->texpMuon();
            m_texpi_1 = (*iter_tof)->texpPion();
            m_texk_1 = (*iter_tof)->texpKaon();
            m_texp_1 = (*iter_tof)->texpProton();
            m_dte_1 = m_tof_1 - m_texe_1;
            m_dtmu_1 = m_tof_1 - m_texmu_1;
            m_dtpi_1 = m_tof_1 - m_texpi_1;
            m_dtk_1 = m_tof_1 - m_texk_1;
            m_dtp_1 = m_tof_1 - m_texp_1;
            //   std::cout<<"There are " <<nTofInfo <<" Tof info"<< std::endl;
            //   why n is almost equal 7 ;is it require to modify?
            nTofInfo_1 += 1;
          }
          m_nTofInfo_1 = nTofInfo_1;
        }
        // EMC information for the first track
        m_evp_1 = evp;
        m_ene_1 = eemc;

        // dedx information for the first track
        if ((*itTrk)->isMdcDedxValid()) {
          RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();
          m_chi_e_1 = dedxTrk->chiE();
          m_chi_mu_1 = dedxTrk->chiMu();
          m_chi_pi_1 = dedxTrk->chiPi();
          m_chi_k_1 = dedxTrk->chiK();
          m_chi_p_1 = dedxTrk->chiP();
        } else {
          m_chi_e_1 = -100.;
          m_chi_mu_1 = -100.;
          m_chi_pi_1 = -100.;
          m_chi_k_1 = -100.;
          m_chi_p_1 = -100.;
        }

        // MUC information for the first track
        if ((*itTrk)->isMucTrackValid()) {
          RecMucTrack *mucTrk = (*itTrk)->mucTrack();
          m_maxhitsinlay_1 = mucTrk->maxHitsInLayer();
          m_numhits_1 = mucTrk->numHits();
          m_numlayers_1 = mucTrk->numLayers();
          m_depth_1 = mucTrk->depth();
          m_mucchi2_1 = mucTrk->chi2();
        } else {
          m_numhits_1 = -10;
          m_maxhitsinlay_1 = -10;
          m_numlayers_1 = -10;
          m_depth_1 = -10;
          m_mucchi2_1 = 300;
        }
      } // the end of the information of the first track
      else {
        // Mdc information for the second track
        RecMdcKalTrack *mdckaltrk = (*itTrk)->mdcKalTrack();
        for (int j = 0; j < 5; j++) {
          RecMdcKalTrack::setPidType(pidtype[j]);
          HepLorentzVector ptrk;
          ptrk.setPx(mdckaltrk->px());
          ptrk.setPy(mdckaltrk->py());
          ptrk.setPz(mdckaltrk->pz());
          double p3 = ptrk.mag();
          ptrk.setE(sqrt(p3 * p3 + xmass[j] * xmass[j]));
          P4_2[j] = ptrk;
        }

        m_charge_2 = mdckaltrk->charge();
        theta_2 = mdcTrk->theta();
        phi_2 = mdcTrk->phi();

        m_pmdc_2 = mdcTrk->p();
        m_ptmdc_2 = mdcTrk->pxy();
        for (int j = 0; j < 5; j++) {
          m_pkal_2[j] = P4_2[j].rho();
          m_ptkal_2[j] = P4_2[j].perp();
        }

        m_theta_2 = thetamdc;
        m_p3_2[0] = mdcTrk->px();
        m_p3_2[1] = mdcTrk->py();
        m_p3_2[2] = mdcTrk->pz();
        P3_2.set(mdcTrk->px(), mdcTrk->py(), mdcTrk->pz());

        // PID information for the second track
        ParticleID *pid = ParticleID::instance();
        pid->init();
        pid->setMethod(pid->methodProbability());
        pid->setChiMinCut(4);
        pid->setRecTrack(*itTrk);
        pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() |
                       pid->useEmc());
        pid->identify(pid->onlyProton() | pid->onlyKaon() | pid->onlyPion() |
                      pid->onlyElectron() | pid->onlyMuon());
        pid->calculate();
        if (pid->IsPidInfoValid()) {
          m_probe_2 = pid->probElectron();
          m_probmu_2 = pid->probMuon();
          m_probpi_2 = pid->probPion();
          m_probk_2 = pid->probKaon();
          m_probp_2 = pid->probProton();
        } else {
          m_probe_2 = -1.0;
          m_probmu_2 = -1.0;
          m_probpi_2 = -1.0;
          m_probk_2 = -1.0;
          m_probp_2 = -1.0;
        }

        double tmp_probe_2 = m_probe_2;
        double tmp_probmu_2 = m_probmu_2;
        double tmp_probpi_2 = m_probpi_2;

        if (tmp_probpi_2 > 0.001 && tmp_probpi_2 > tmp_probe_2 &&
            tmp_probpi_2 > tmp_probmu_2) {
          npi_pid++;
        }
        if (tmp_probe_2 > 0.001 && tmp_probe_2 > tmp_probpi_2 &&
            tmp_probe_2 > tmp_probmu_2) {
          nele_pid++;
        }
        if (tmp_probmu_2 > 0.001 && tmp_probmu_2 > tmp_probpi_2 &&
            tmp_probmu_2 > tmp_probe_2) {
          nmu_pid++;
        }

        // TOF information for the second track
        if ((*itTrk)->isTofTrackValid()) {
          SmartRefVector<RecTofTrack> tofTrkCol = (*itTrk)->tofTrack();
          SmartRefVector<RecTofTrack>::iterator iter_tof = tofTrkCol.begin();
          int nTofInfo_2 = 0;
          //@NOTE: The last value in tofTrkCol is the average value.
          for (; iter_tof != tofTrkCol.end(); iter_tof++) {
            TofHitStatus *status = new TofHitStatus;
            status->setStatus((*iter_tof)->status());
            m_fgtof_2 = 1;
            m_counter_2 = status->is_counter();
            m_isbarrel_2 = status->is_barrel();
            m_layertof_2 = status->layer();
            m_iscluster_2 = status->is_cluster();
            m_tof_2 = (*iter_tof)->tof();
            m_texe_2 = (*iter_tof)->texpElectron();
            m_texmu_2 = (*iter_tof)->texpMuon();
            m_texpi_2 = (*iter_tof)->texpPion();
            m_texk_2 = (*iter_tof)->texpKaon();
            m_texp_2 = (*iter_tof)->texpProton();
            m_dte_2 = m_tof_2 - m_texe_2;
            m_dtmu_2 = m_tof_2 - m_texmu_2;
            m_dtpi_2 = m_tof_2 - m_texpi_2;
            m_dtk_2 = m_tof_2 - m_texk_2;
            m_dtp_2 = m_tof_2 - m_texp_2;
            nTofInfo_2 += 1;
          }
          m_nTofInfo_2 = nTofInfo_2;
        }
        // EMC information for the second track
        m_evp_2 = evp;
        m_ene_2 = eemc;

        // dedx information for the second track
        if ((*itTrk)->isMdcDedxValid()) {
          RecMdcDedx *dedxTrk = (*itTrk)->mdcDedx();
          m_chi_e_2 = dedxTrk->chiE();
          m_chi_mu_2 = dedxTrk->chiMu();
          m_chi_pi_2 = dedxTrk->chiPi();
          m_chi_k_2 = dedxTrk->chiK();
          m_chi_p_2 = dedxTrk->chiP();
        } else {
          m_chi_e_2 = -100.;
          m_chi_mu_2 = -100.;
          m_chi_pi_2 = -100.;
          m_chi_k_2 = -100.;
          m_chi_p_2 = -100.;
        }

        // MUC information for the second track
        if ((*itTrk)->isMucTrackValid()) {
          RecMucTrack *mucTrk = (*itTrk)->mucTrack();
          m_maxhitsinlay_2 = mucTrk->maxHitsInLayer();
          m_numhits_2 = mucTrk->numHits();
          m_numlayers_2 = mucTrk->numLayers();
          m_depth_2 = mucTrk->depth();
          m_mucchi2_2 = mucTrk->chi2();
        } else {
          m_numhits_2 = -10;
          m_maxhitsinlay_2 = -10;
          m_numlayers_2 = -10;
          m_depth_2 = -10;
          m_mucchi2_2 = 300;
        }
      } // the end of the information of the second track
    }   // If EmcShowerValid
  }     // the end of the information of charged tracks

  int nchrgtrk = itrkWithValidEmcShower.size();

  m_nTrk_EMC = nchrgtrk;
  if (nchrgtrk != 2) {
    return StatusCode::SUCCESS;
  }
  m_cutflow->Fill(4);
  Ncut[4]++;
 
 
  if (itrkWithValidEmcShower[0] == itrkWithValidEmcShower[1]) {
    throw std::runtime_error("This is not supposed to happen");
  }

  if (!(m_charge_1 * m_charge_2 < 0)) {
    throw std::runtime_error("This is not supposed to happen");
  }


  // 2024.04.25,remove two-pi or two-lep event
  if (npi_pid != 1 || (nele_pid + nmu_pid) != 1) {
    return StatusCode::SUCCESS;
  }
  m_cutflow->Fill(5);
  Ncut[5]++;

  //-------------- record angle and invarm ------------------
  Hep3Vector tmp1_p[4];
  Hep3Vector tmp2_p[4];
  Hep3Vector P3_4type[4];
  Vp4 invar;
  for (int i = 0; i < nGood; i++) {
    EvtRecTrackIterator itTrka = evtRecTrkCol->begin() + iGood[i];
    RecMdcTrack *mdcTrka = (*itTrka)->mdcTrack();
    double p_tmp = mdcTrka->p();
    double pt_tmp = mdcTrka->pxy();
    Hep3Vector P3_4type[4];
    invar.clear();
    RecMdcKalTrack *mdcKalTrka = (*itTrka)->mdcKalTrack();
    if (mdcKalTrka != nullptr) {
      for (int m = 0; m < 4; m++) {
        RecMdcKalTrack::setPidType(pidtype[m]);
        Hep3Vector tmp_trk;
        tmp_trk.setX(mdcKalTrka->px());
        tmp_trk.setY(mdcKalTrka->py());
        tmp_trk.setZ(mdcKalTrka->pz());
        P3_4type[m] = tmp_trk;
        double charge_e = sqrt(tmp_trk * tmp_trk + xmass[m] * xmass[m]);
        HepLorentzVector ctrk(tmp_trk, charge_e);
        invar.push_back(ctrk);
      }
    } else {
      for (int p = 0; p < 4; p++) {
        Hep3Vector null_trk(-10., 10., -10.);
        P3_4type[p] = null_trk;
        HepLorentzVector ftrk(null_trk, -10.);
        invar.push_back(ftrk);
      }
    }

    if (i == 0) {
      for (int r = 0; r < 4; r++) {
        tmp1_p[r] = P3_4type[r];
      }
      m_px_pos_e = P3_4type[0].x();
      m_py_pos_e = P3_4type[0].y();
      m_pz_pos_e = P3_4type[0].z();
      m_px_pos_mu = P3_4type[1].x();
      m_py_pos_mu = P3_4type[1].y();
      m_pz_pos_mu = P3_4type[1].z();
      m_px_pos_pi = P3_4type[2].x();
      m_py_pos_pi = P3_4type[2].y();
      m_pz_pos_pi = P3_4type[2].z();
      m_px_pos_k = P3_4type[3].x();
      m_py_pos_k = P3_4type[3].y();
      m_pz_pos_k = P3_4type[3].z();
      m_p_pos = p_tmp;
      m_pt_pos = pt_tmp;
    } else if (i == 1) {
      for (int s = 0; s < 4; s++) {
        tmp1_p[s] = P3_4type[s];
      }
      m_px_neg_e = P3_4type[0].x();
      m_py_neg_e = P3_4type[0].y();
      m_pz_neg_e = P3_4type[0].z();
      m_px_neg_mu = P3_4type[1].x();
      m_py_neg_mu = P3_4type[1].y();
      m_pz_neg_mu = P3_4type[1].z();
      m_px_neg_pi = P3_4type[2].x();
      m_py_neg_pi = P3_4type[2].y();
      m_pz_neg_pi = P3_4type[2].z();
      m_px_neg_k = P3_4type[3].x();
      m_py_neg_k = P3_4type[3].y();
      m_pz_neg_k = P3_4type[3].z();
      m_p_neg = p_tmp;
      m_pt_neg = pt_tmp;
    }

    for (int j = 0; j < nGam; j++) {
      EvtRecTrackIterator jtTrka = evtRecTrkCol->begin() + iGam[j];
      if (!(*jtTrka)->isEmcShowerValid()) {
        throw std::runtime_error("This should not happen. The selected good "
                                 "photons are already required to have valid "
                                 "EmcShower");
      }
      RecEmcShower *emcTrka = (*jtTrka)->emcShower();
      double an_eraw = emcTrka->energy();
      double an_phi = emcTrka->phi();
      double an_the = emcTrka->theta();
      Hep3Vector gam_p(an_eraw * sin(an_the) * cos(an_phi),
                       an_eraw * sin(an_the) * sin(an_phi),
                       an_eraw * cos(an_the));
      HepLorentzVector gtrk(gam_p, an_eraw);
      m_eraw[j] = an_eraw;
      m_phi[j] = an_phi;
      m_the[j] = an_the;

      double angd_e = P3_4type[0].angle(gam_p);
      double thed_e = P3_4type[0].theta() - gam_p.theta();
      double phid_e = P3_4type[0].deltaPhi(gam_p);
      double angd_mu = P3_4type[1].angle(gam_p);
      double thed_mu = P3_4type[1].theta() - gam_p.theta();
      double phid_mu = P3_4type[1].deltaPhi(gam_p);
      double angd_pi = P3_4type[2].angle(gam_p);
      double thed_pi = P3_4type[2].theta() - gam_p.theta();
      double phid_pi = P3_4type[2].deltaPhi(gam_p);
      double angd_k = P3_4type[2].angle(gam_p);
      double thed_k = P3_4type[2].theta() - gam_p.theta();
      double phid_k = P3_4type[2].deltaPhi(gam_p);
      thed_e = fmod(thed_e + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
               CLHEP::pi;
      phid_e = fmod(phid_e + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
               CLHEP::pi;
      thed_mu = fmod(thed_mu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
                CLHEP::pi;
      phid_mu = fmod(phid_mu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
                CLHEP::pi;
      thed_pi = fmod(thed_pi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
                CLHEP::pi;
      phid_pi = fmod(phid_pi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
                CLHEP::pi;
      thed_k = fmod(thed_k + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
               CLHEP::pi;
      phid_k = fmod(phid_k + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
               CLHEP::pi;
      double dang_e = angd_e * 180 / (CLHEP::pi);
      double dthe_e = thed_e * 180 / (CLHEP::pi);
      double dphi_e = phid_e * 180 / (CLHEP::pi);
      double dang_mu = angd_mu * 180 / (CLHEP::pi);
      double dthe_mu = thed_mu * 180 / (CLHEP::pi);
      double dphi_mu = phid_mu * 180 / (CLHEP::pi);
      double dang_pi = angd_pi * 180 / (CLHEP::pi);
      double dthe_pi = thed_pi * 180 / (CLHEP::pi);
      double dphi_pi = phid_pi * 180 / (CLHEP::pi);
      double dang_k = angd_k * 180 / (CLHEP::pi);
      double dthe_k = thed_k * 180 / (CLHEP::pi);
      double dphi_k = phid_k * 180 / (CLHEP::pi);

      if (i == 0) {
        m_dang_pos_e[j] = dang_e;
        m_dthe_pos_e[j] = dthe_e;
        m_dphi_pos_e[j] = dphi_e;
        m_dang_pos_mu[j] = dang_mu;
        m_dthe_pos_mu[j] = dthe_mu;
        m_dphi_pos_mu[j] = dphi_mu;
        m_dang_pos_pi[j] = dang_pi;
        m_dthe_pos_pi[j] = dthe_pi;
        m_dphi_pos_pi[j] = dphi_pi;
        m_dang_pos_k[j] = dang_k;
        m_dthe_pos_k[j] = dthe_k;
        m_dphi_pos_k[j] = dphi_k;
        for (int l = 0; l < 4; l++) {
          if (l == 0) {
            m_m_egam_pos[j] = (invar[l] + gtrk).m();
          }
          if (l == 1) {
            m_m_mugam_pos[j] = (invar[l] + gtrk).m();
          }
          if (l == 2) {
            m_m_pigam_pos[j] = (invar[l] + gtrk).m();
          }
          if (l == 3) {
            m_m_kgam_pos[j] = (invar[l] + gtrk).m();
          }
        }
      } else if (i == 1) {
        m_dang_neg_e[j] = dang_e;
        m_dthe_neg_e[j] = dthe_e;
        m_dphi_neg_e[j] = dphi_e;
        m_dang_neg_mu[j] = dang_mu;
        m_dthe_neg_mu[j] = dthe_mu;
        m_dphi_neg_mu[j] = dphi_mu;
        m_dang_neg_pi[j] = dang_pi;
        m_dthe_neg_pi[j] = dthe_pi;
        m_dphi_neg_pi[j] = dphi_pi;
        m_dang_neg_k[j] = dang_k;
        m_dthe_neg_k[j] = dthe_k;
        m_dphi_neg_k[j] = dphi_k;
        for (int l = 0; l < 4; l++) {
          if (l == 0) {
            m_m_egam_neg[j] = (invar[l] + gtrk).m();
          }
          if (l == 1) {
            m_m_mugam_neg[j] = (invar[l] + gtrk).m();
          }
          if (l == 2) {
            m_m_pigam_neg[j] = (invar[l] + gtrk).m();
          }
          if (l == 3) {
            m_m_kgam_neg[j] = (invar[l] + gtrk).m();
          }
        }
      }
    }
  }
  double angc_ee = tmp1_p[0].angle(tmp2_p[0]);
  double thec_ee = tmp1_p[0].theta() - tmp2_p[0].theta();
  double phic_ee = tmp1_p[0].deltaPhi(tmp2_p[0]);
  double theadd_ee = tmp1_p[0].theta() + tmp2_p[0].theta();
  double cos_cal_ee = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_emu = tmp1_p[0].angle(tmp2_p[1]);
  double thec_emu = tmp1_p[0].theta() - tmp2_p[1].theta();
  double phic_emu = tmp1_p[0].deltaPhi(tmp2_p[1]);
  double theadd_emu = tmp1_p[0].theta() + tmp2_p[1].theta();
  double cos_cal_emu = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_epi = tmp1_p[0].angle(tmp2_p[2]);
  double thec_epi = tmp1_p[0].theta() - tmp2_p[2].theta();
  double phic_epi = tmp1_p[0].deltaPhi(tmp2_p[2]);
  double theadd_epi = tmp1_p[0].theta() + tmp2_p[2].theta();
  double cos_cal_epi = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_ek = tmp1_p[0].angle(tmp2_p[3]);
  double thec_ek = tmp1_p[0].theta() - tmp2_p[3].theta();
  double phic_ek = tmp1_p[0].deltaPhi(tmp2_p[3]);
  double theadd_ek = tmp1_p[0].theta() + tmp2_p[3].theta();
  double cos_cal_ek = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_mue = tmp1_p[1].angle(tmp2_p[0]);
  double thec_mue = tmp1_p[1].theta() - tmp2_p[0].theta();
  double phic_mue = tmp1_p[1].deltaPhi(tmp2_p[0]);
  double theadd_mue = tmp1_p[1].theta() + tmp2_p[0].theta();
  double cos_cal_mue = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_mumu = tmp1_p[1].angle(tmp2_p[1]);
  double thec_mumu = tmp1_p[1].theta() - tmp2_p[1].theta();
  double phic_mumu = tmp1_p[1].deltaPhi(tmp2_p[1]);
  double theadd_mumu = tmp1_p[1].theta() + tmp2_p[1].theta();
  double cos_cal_mumu = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_mupi = tmp1_p[1].angle(tmp2_p[2]);
  double thec_mupi = tmp1_p[1].theta() - tmp2_p[2].theta();
  double phic_mupi = tmp1_p[1].deltaPhi(tmp2_p[2]);
  double theadd_mupi = tmp1_p[1].theta() + tmp2_p[2].theta();
  double cos_cal_mupi = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_muk = tmp1_p[1].angle(tmp2_p[3]);
  double thec_muk = tmp1_p[1].theta() - tmp2_p[3].theta();
  double phic_muk = tmp1_p[1].deltaPhi(tmp2_p[3]);
  double theadd_muk = tmp1_p[1].theta() + tmp2_p[3].theta();
  double cos_cal_muk = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_pie = tmp1_p[2].angle(tmp2_p[0]);
  double thec_pie = tmp1_p[2].theta() - tmp2_p[0].theta();
  double phic_pie = tmp1_p[2].deltaPhi(tmp2_p[0]);
  double theadd_pie = tmp1_p[2].theta() + tmp2_p[0].theta();
  double cos_cal_pie = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_pimu = tmp1_p[2].angle(tmp2_p[1]);
  double thec_pimu = tmp1_p[2].theta() - tmp2_p[1].theta();
  double phic_pimu = tmp1_p[2].deltaPhi(tmp2_p[1]);
  double theadd_pimu = tmp1_p[2].theta() + tmp2_p[1].theta();
  double cos_cal_pimu = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_pipi = tmp1_p[2].angle(tmp2_p[2]);
  double thec_pipi = tmp1_p[2].theta() - tmp2_p[2].theta();
  double phic_pipi = tmp1_p[2].deltaPhi(tmp2_p[2]);
  double theadd_pipi = tmp1_p[2].theta() + tmp2_p[2].theta();
  double cos_cal_pipi = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_pik = tmp1_p[2].angle(tmp2_p[3]);
  double thec_pik = tmp1_p[2].theta() - tmp2_p[3].theta();
  double phic_pik = tmp1_p[2].deltaPhi(tmp2_p[3]);
  double theadd_pik = tmp1_p[2].theta() + tmp2_p[3].theta();
  double cos_cal_pik = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_ke = tmp1_p[3].angle(tmp2_p[0]);
  double thec_ke = tmp1_p[3].theta() - tmp2_p[0].theta();
  double phic_ke = tmp1_p[3].deltaPhi(tmp2_p[0]);
  double theadd_ke = tmp1_p[3].theta() + tmp2_p[0].theta();
  double cos_cal_ke = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_kmu = tmp1_p[3].angle(tmp2_p[1]);
  double thec_kmu = tmp1_p[3].theta() - tmp2_p[1].theta();
  double phic_kmu = tmp1_p[3].deltaPhi(tmp2_p[1]);
  double theadd_kmu = tmp1_p[3].theta() + tmp2_p[1].theta();
  double cos_cal_kmu = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_kpi = tmp1_p[3].angle(tmp2_p[2]);
  double thec_kpi = tmp1_p[3].theta() - tmp2_p[2].theta();
  double phic_kpi = tmp1_p[3].deltaPhi(tmp2_p[2]);
  double theadd_kpi = tmp1_p[3].theta() + tmp2_p[2].theta();
  double cos_cal_kpi = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  double angc_kk = tmp1_p[3].angle(tmp2_p[3]);
  double thec_kk = tmp1_p[3].theta() - tmp2_p[3].theta();
  double phic_kk = tmp1_p[3].deltaPhi(tmp2_p[3]);
  double theadd_kk = tmp1_p[3].theta() + tmp2_p[3].theta();
  double cos_cal_kk = sin(thec_ee / 2.0) / (sin(theadd_ee / 2.0));
  thec_ee = fmod(thec_ee + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  phic_ee = fmod(phic_ee + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  thec_emu = fmod(thec_emu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_emu = fmod(phic_emu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_epi = fmod(thec_epi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_epi = fmod(phic_epi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_ek = fmod(thec_ek + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  phic_ek = fmod(phic_ek + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  thec_mue = fmod(thec_mue + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_mue = fmod(phic_mue + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_mumu = fmod(thec_mumu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  phic_mumu = fmod(phic_mumu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  thec_mupi = fmod(thec_mupi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  phic_mupi = fmod(phic_mupi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  thec_muk = fmod(thec_muk + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_muk = fmod(phic_muk + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_pie = fmod(thec_pie + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_pie = fmod(phic_pie + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_pimu = fmod(thec_pimu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  phic_pimu = fmod(phic_pimu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  thec_pipi = fmod(thec_pipi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  phic_pipi = fmod(phic_pipi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
              CLHEP::pi;
  thec_pik = fmod(thec_pik + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_pik = fmod(phic_pik + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_ke = fmod(thec_ke + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  phic_ke = fmod(phic_ke + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  thec_kmu = fmod(thec_kmu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_kmu = fmod(phic_kmu + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_kpi = fmod(thec_kpi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  phic_kpi = fmod(phic_kpi + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi;
  thec_kk = fmod(thec_kk + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;
  phic_kk = fmod(phic_kk + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
            CLHEP::pi;

  m_costheta_ee = cos_cal_ee;
  m_cang_ee = angc_ee * 180 / (CLHEP::pi);
  m_cthe_ee = thec_ee * 180 / (CLHEP::pi);
  m_cphi_ee = phic_ee * 180 / (CLHEP::pi);
  m_costheta_emu = cos_cal_emu;
  m_cang_emu = angc_emu * 180 / (CLHEP::pi);
  m_cthe_emu = thec_emu * 180 / (CLHEP::pi);
  m_cphi_emu = phic_emu * 180 / (CLHEP::pi);
  m_costheta_epi = cos_cal_epi;
  m_cang_epi = angc_epi * 180 / (CLHEP::pi);
  m_cthe_epi = thec_epi * 180 / (CLHEP::pi);
  m_cphi_epi = phic_epi * 180 / (CLHEP::pi);
  m_costheta_ek = cos_cal_ek;
  m_cang_ek = angc_ek * 180 / (CLHEP::pi);
  m_cthe_ek = thec_ek * 180 / (CLHEP::pi);
  m_cphi_ek = phic_ek * 180 / (CLHEP::pi);
  m_costheta_mue = cos_cal_mue;
  m_cang_mue = angc_mue * 180 / (CLHEP::pi);
  m_cthe_mue = thec_mue * 180 / (CLHEP::pi);
  m_cphi_mue = phic_mue * 180 / (CLHEP::pi);
  m_costheta_mumu = cos_cal_mumu;
  m_cang_mumu = angc_mumu * 180 / (CLHEP::pi);
  m_cthe_mumu = thec_mumu * 180 / (CLHEP::pi);
  m_cphi_mumu = phic_mumu * 180 / (CLHEP::pi);
  m_costheta_mupi = cos_cal_mupi;
  m_cang_mupi = angc_mupi * 180 / (CLHEP::pi);
  m_cthe_mupi = thec_mupi * 180 / (CLHEP::pi);
  m_cphi_mupi = phic_mupi * 180 / (CLHEP::pi);
  m_costheta_muk = cos_cal_muk;
  m_cang_muk = angc_muk * 180 / (CLHEP::pi);
  m_cthe_muk = thec_muk * 180 / (CLHEP::pi);
  m_cphi_muk = phic_muk * 180 / (CLHEP::pi);
  m_costheta_pie = cos_cal_pie;
  m_cang_pie = angc_pie * 180 / (CLHEP::pi);
  m_cthe_pie = thec_pie * 180 / (CLHEP::pi);
  m_cphi_pie = phic_pie * 180 / (CLHEP::pi);
  m_costheta_pimu = cos_cal_pimu;
  m_cang_pimu = angc_pimu * 180 / (CLHEP::pi);
  m_cthe_pimu = thec_pimu * 180 / (CLHEP::pi);
  m_cphi_pimu = phic_pimu * 180 / (CLHEP::pi);
  m_costheta_pipi = cos_cal_pipi;
  m_cang_pipi = angc_pipi * 180 / (CLHEP::pi);
  m_cthe_pipi = thec_pipi * 180 / (CLHEP::pi);
  m_cphi_pipi = phic_pipi * 180 / (CLHEP::pi);
  m_costheta_pik = cos_cal_pik;
  m_cang_pik = angc_pik * 180 / (CLHEP::pi);
  m_cthe_pik = thec_pik * 180 / (CLHEP::pi);
  m_cphi_pik = phic_pik * 180 / (CLHEP::pi);
  m_costheta_ke = cos_cal_ke;
  m_cang_ke = angc_ke * 180 / (CLHEP::pi);
  m_cthe_ke = thec_ke * 180 / (CLHEP::pi);
  m_cphi_ke = phic_ke * 180 / (CLHEP::pi);
  m_costheta_kmu = cos_cal_kmu;
  m_cang_kmu = angc_kmu * 180 / (CLHEP::pi);
  m_cthe_kmu = thec_kmu * 180 / (CLHEP::pi);
  m_cphi_kmu = phic_kmu * 180 / (CLHEP::pi);
  m_costheta_kpi = cos_cal_kpi;
  m_cang_kpi = angc_kpi * 180 / (CLHEP::pi);
  m_cthe_kpi = thec_kpi * 180 / (CLHEP::pi);
  m_cphi_kpi = phic_kpi * 180 / (CLHEP::pi);
  m_costheta_kk = cos_cal_kk;
  m_cang_kk = angc_kk * 180 / (CLHEP::pi);
  m_cthe_kk = thec_kk * 180 / (CLHEP::pi);
  m_cphi_kk = phic_kk * 180 / (CLHEP::pi);

  // calculate total momentum of all tracks and photons
  HepLorentzVector pkal_tot[5];
  HepLorentzVector pkal_Trks[5];
  for (int i = 0; i < 5; i++) {
    pkal_tot[i] = P4_1[i] + P4_2[i];
    pkal_Trks[i] = P4_1[i] + P4_2[i];
    for (int j = 0; j < nGam; j++) {
      pkal_tot[i] += pGam[j];
    }
  }

  for (int i = 0; i < 5; i++) {
    m_pkal_tot_px[i] = pkal_tot[i].px();
    m_pkal_tot_py[i] = pkal_tot[i].py();
    m_pkal_tot_pz[i] = pkal_tot[i].pz();
    m_pkal_tot_e[i] = pkal_tot[i].e();
    m_mTrks[i] = pkal_Trks[i].m();
    m_m2Trks[i] = pkal_Trks[i].m2();
  }

  m_ang_mdc_acol = P3_1.polarAngle(P3_2) * 180 / CLHEP::pi;
  m_ang_mdc_acop = P3_1.deltaPhi(P3_2) * 180 / CLHEP::pi;

  m_the_add = theta_1 + theta_2;
  m_phi_diff = phi_1 - phi_2;
  if (m_phi_diff < 0)
    m_phi_diff = m_phi_diff + twopi;
  m_PTEM = (P3_1 + P3_2).perp() / (m_ecms - (P3_1 + P3_2).r());

  HepLorentzVector ecms(0.011 * m_ecms, 0, 0, m_ecms);
  // boost to cms
  Hep3Vector P3_trk1[5];
  Hep3Vector P3_trk2[5];
  for (int i = 0; i < 5; i++) {
    /*		  P4_1[i].boost(-0.011,0,0);
                      P4_2[i].boost(-0.011,0,0);
                      cms_P4_1[i]=P4_1[i];
                      cms_P4_2[i]=P4_2[i];
                      P3_trk1[i] = P4_1[i].vect();
                      P3_trk2[i] = P4_2[i].vect();
                      //boost to lab system
                      P4_1[i].boost(0.011,0,0);
                      P4_2[i].boost(0.011,0,0);  */

    cms_P4_1[i] = P4_1[i];
    cms_P4_2[i] = P4_2[i];
    cms_P4_1[i].boost(-0.011, 0, 0);
    cms_P4_2[i].boost(-0.011, 0, 0);
    P3_trk1[i] = cms_P4_1[i].vect();
    P3_trk2[i] = cms_P4_2[i].vect();
  }
  HepLorentzVector p_all(0, 0, 0, m_ecms);
  HepLorentzVector P4_ee = cms_P4_1[0] + cms_P4_2[0];
  HepLorentzVector P4_emu = cms_P4_1[0] + cms_P4_2[1];
  HepLorentzVector P4_mue = cms_P4_1[1] + cms_P4_2[0];
  HepLorentzVector P4_epi = cms_P4_1[0] + cms_P4_2[2];
  HepLorentzVector P4_pie = cms_P4_1[2] + cms_P4_2[0];
  HepLorentzVector P4_ek = cms_P4_1[0] + cms_P4_2[3];
  HepLorentzVector P4_ke = cms_P4_1[3] + cms_P4_2[0];
  HepLorentzVector P4_mumu = cms_P4_1[1] + cms_P4_2[1];
  HepLorentzVector P4_mupi = cms_P4_1[1] + cms_P4_2[2];
  HepLorentzVector P4_pimu = cms_P4_1[2] + cms_P4_2[1];
  HepLorentzVector P4_pipi = cms_P4_1[2] + cms_P4_2[2];
  HepLorentzVector P4_muk = cms_P4_1[1] + cms_P4_2[3];
  HepLorentzVector P4_kmu = cms_P4_1[3] + cms_P4_2[1];
  HepLorentzVector P4_pik = cms_P4_1[2] + cms_P4_2[3];
  HepLorentzVector P4_kpi = cms_P4_1[3] + cms_P4_2[2];
  HepLorentzVector P4_kk = cms_P4_1[3] + cms_P4_2[3];
  m_miss_m2_ee = (m_ecms - P4_ee.e()) * (m_ecms - P4_ee.e()) -
                 ((p_all - P4_ee).rho()) * ((p_all - P4_ee).rho());
  m_miss_m2_emu = (m_ecms - P4_emu.e()) * (m_ecms - P4_emu.e()) -
                  ((p_all - P4_emu).rho()) * ((p_all - P4_emu).rho());
  m_miss_m2_mue = (m_ecms - P4_mue.e()) * (m_ecms - P4_mue.e()) -
                  ((p_all - P4_mue).rho()) * ((p_all - P4_mue).rho());
  m_miss_m2_epi = (m_ecms - P4_epi.e()) * (m_ecms - P4_epi.e()) -
                  ((p_all - P4_epi).rho()) * ((p_all - P4_epi).rho());
  m_miss_m2_pie = (m_ecms - P4_pie.e()) * (m_ecms - P4_pie.e()) -
                  ((p_all - P4_pie).rho()) * ((p_all - P4_pie).rho());
  m_miss_m2_ek = (m_ecms - P4_ek.e()) * (m_ecms - P4_ek.e()) -
                 ((p_all - P4_ek).rho()) * ((p_all - P4_ek).rho());
  m_miss_m2_ke = (m_ecms - P4_ke.e()) * (m_ecms - P4_ke.e()) -
                 ((p_all - P4_ke).rho()) * ((p_all - P4_ke).rho());
  m_miss_m2_mumu = (m_ecms - P4_mumu.e()) * (m_ecms - P4_mumu.e()) -
                   ((p_all - P4_mumu).rho()) * ((p_all - P4_mumu).rho());
  m_miss_m2_mupi = (m_ecms - P4_mupi.e()) * (m_ecms - P4_mupi.e()) -
                   ((p_all - P4_mupi).rho()) * ((p_all - P4_mupi).rho());
  m_miss_m2_pimu = (m_ecms - P4_pimu.e()) * (m_ecms - P4_pimu.e()) -
                   ((p_all - P4_pimu).rho()) * ((p_all - P4_pimu).rho());
  m_miss_m2_pipi = (m_ecms - P4_pipi.e()) * (m_ecms - P4_pipi.e()) -
                   ((p_all - P4_pipi).rho()) * ((p_all - P4_pipi).rho());
  m_miss_m2_muk = (m_ecms - P4_muk.e()) * (m_ecms - P4_muk.e()) -
                  ((p_all - P4_muk).rho()) * ((p_all - P4_muk).rho());
  m_miss_m2_kmu = (m_ecms - P4_kmu.e()) * (m_ecms - P4_kmu.e()) -
                  ((p_all - P4_kmu).rho()) * ((p_all - P4_kmu).rho());
  m_miss_m2_pik = (m_ecms - P4_pik.e()) * (m_ecms - P4_pik.e()) -
                  ((p_all - P4_pik).rho()) * ((p_all - P4_pik).rho());
  m_miss_m2_kpi = (m_ecms - P4_kpi.e()) * (m_ecms - P4_kpi.e()) -
                  ((p_all - P4_kpi).rho()) * ((p_all - P4_kpi).rho());
  m_miss_m2_kk = (m_ecms - P4_kk.e()) * (m_ecms - P4_kk.e()) -
                 ((p_all - P4_kk).rho()) * ((p_all - P4_kk).rho());

  m_ee_angle = P3_trk1[0].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_emu_angle = P3_trk1[0].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mue_angle = P3_trk1[1].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_epi_angle = P3_trk1[0].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pie_angle = P3_trk1[2].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_ek_angle = P3_trk1[0].angle(P3_trk2[3]) * 180 / CLHEP::pi;
  m_ke_angle = P3_trk1[3].angle(P3_trk2[0]) * 180 / CLHEP::pi;
  m_mumu_angle = P3_trk1[1].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mupi_angle = P3_trk1[1].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pimu_angle = P3_trk1[2].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pipi_angle = P3_trk1[2].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_muk_angle = P3_trk1[1].angle(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kmu_angle = P3_trk1[3].angle(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pik_angle = P3_trk1[2].angle(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kpi_angle = P3_trk1[3].angle(P3_trk2[2]) * 180 / CLHEP::pi;
  m_kk_angle = P3_trk1[3].angle(P3_trk2[3]) * 180 / CLHEP::pi;

  m_ee_acop = P3_trk1[0].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_emu_acop = P3_trk1[0].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mue_acop = P3_trk1[1].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_epi_acop = P3_trk1[0].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pie_acop = P3_trk1[2].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_ek_acop = P3_trk1[0].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;
  m_ke_acop = P3_trk1[3].deltaPhi(P3_trk2[0]) * 180 / CLHEP::pi;
  m_mumu_acop = P3_trk1[1].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_mupi_acop = P3_trk1[1].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_pimu_acop = P3_trk1[2].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pipi_acop = P3_trk1[2].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_muk_acop = P3_trk1[1].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kmu_acop = P3_trk1[3].deltaPhi(P3_trk2[1]) * 180 / CLHEP::pi;
  m_pik_acop = P3_trk1[2].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;
  m_kpi_acop = P3_trk1[3].deltaPhi(P3_trk2[2]) * 180 / CLHEP::pi;
  m_kk_acop = P3_trk1[3].deltaPhi(P3_trk2[3]) * 180 / CLHEP::pi;

  m_ee_PTEM = (P3_trk1[0] + P3_trk2[0]).perp() /
              (m_ecms - (P3_trk1[0] + P3_trk2[0]).r());
  m_emu_PTEM = (P3_trk1[0] + P3_trk2[1]).perp() /
               (m_ecms - (P3_trk1[0] + P3_trk2[1]).r());
  m_mue_PTEM = (P3_trk1[1] + P3_trk2[0]).perp() /
               (m_ecms - (P3_trk1[1] + P3_trk2[0]).r());
  m_epi_PTEM = (P3_trk1[0] + P3_trk2[2]).perp() /
               (m_ecms - (P3_trk1[0] + P3_trk2[2]).r());
  m_pie_PTEM = (P3_trk1[2] + P3_trk2[0]).perp() /
               (m_ecms - (P3_trk1[2] + P3_trk2[0]).r());
  m_ek_PTEM = (P3_trk1[0] + P3_trk2[3]).perp() /
              (m_ecms - (P3_trk1[0] + P3_trk2[3]).r());
  m_ke_PTEM = (P3_trk1[3] + P3_trk2[0]).perp() /
              (m_ecms - (P3_trk1[3] + P3_trk2[0]).r());
  m_mumu_PTEM = (P3_trk1[1] + P3_trk2[1]).perp() /
                (m_ecms - (P3_trk1[1] + P3_trk2[1]).r());
  m_mupi_PTEM = (P3_trk1[1] + P3_trk2[2]).perp() /
                (m_ecms - (P3_trk1[1] + P3_trk2[2]).r());
  m_pimu_PTEM = (P3_trk1[2] + P3_trk2[1]).perp() /
                (m_ecms - (P3_trk1[2] + P3_trk2[1]).r());
  m_pipi_PTEM = (P3_trk1[2] + P3_trk2[2]).perp() /
                (m_ecms - (P3_trk1[2] + P3_trk2[2]).r());
  m_muk_PTEM = (P3_trk1[1] + P3_trk2[3]).perp() /
               (m_ecms - (P3_trk1[1] + P3_trk2[3]).r());
  m_kmu_PTEM = (P3_trk1[3] + P3_trk2[1]).perp() /
               (m_ecms - (P3_trk1[3] + P3_trk2[1]).r());
  m_pik_PTEM = (P3_trk1[2] + P3_trk2[3]).perp() /
               (m_ecms - (P3_trk1[2] + P3_trk2[3]).r());
  m_kpi_PTEM = (P3_trk1[3] + P3_trk2[2]).perp() /
               (m_ecms - (P3_trk1[3] + P3_trk2[2]).r());
  m_kk_PTEM = (P3_trk1[3] + P3_trk2[3]).perp() /
              (m_ecms - (P3_trk1[3] + P3_trk2[3]).r());

  // selection criteria for PTEM
  /* if ( m_epi_PTEM <= 0.05 || m_pie_PTEM <= 0.05)
        return StatusCode::SUCCESS;
     std::cout << "epi_PTEM" << m_epi_PTEM << "  " << m_pie_PTEM << std::endl;
  */

  // record the hit map in Emc
  SmartDataPtr<RecEmcHitCol> recEmcHitCol(eventSvc(),
                                          EventModel::Recon::RecEmcHitCol);
  if (!recEmcHitCol) {
    log << MSG::FATAL << "Could not find emcRecHitCol" << endreq;
    return (StatusCode::FAILURE);
  }
  RecEmcHitCol::iterator iRecEmcHit;
  int nRecEmcHits = 0;
  for (iRecEmcHit = recEmcHitCol->begin(); iRecEmcHit != recEmcHitCol->end();
       iRecEmcHit++) {
    RecEmcHit recEmcHit = *(*iRecEmcHit);
    unsigned int partId = EmcID::barrel_ec(recEmcHit.getCellId());
    unsigned int theta = EmcID::theta_module(recEmcHit.getCellId());
    unsigned int phi = EmcID::phi_module(recEmcHit.getCellId());
    double energy = recEmcHit.getEnergy();
    double time = recEmcHit.getTime();
    HepPoint3D center = recEmcHit.getCenter();
    // m_recEmcHitMap->Fill(theta, phi, energy);
    m_emc_id_phi[nRecEmcHits] = phi;
    m_emc_id_theta[nRecEmcHits] = theta;
    m_emc_energy[nRecEmcHits] = energy;
    m_emc_bc[nRecEmcHits] = partId;
    m_emc_tdc[nRecEmcHits] = time;
    m_emc_x[nRecEmcHits] = center.x();
    m_emc_y[nRecEmcHits] = center.y();
    m_emc_z[nRecEmcHits] = center.z();
    double r = std::sqrt(center.x() * center.x() + center.y() * center.y());
    m_emc_pos_theta[nRecEmcHits] =
        fmod(std::atan2(r, center.z()) + CLHEP::twopi + CLHEP::twopi + pi,
             CLHEP::twopi) -
        CLHEP::pi;
    m_emc_pos_phi[nRecEmcHits] = fmod(std::atan2(center.y(), center.x()) +
                                          CLHEP::twopi + CLHEP::twopi + pi,
                                      CLHEP::twopi) -
                                 CLHEP::pi;

    nRecEmcHits++;
    if (nRecEmcHits >= 1000) {
      break;
    }
  }
  // std::cout<<"nRecEmcHits = " << nRecEmcHits << std::endl;
  m_nRecEmcHits = nRecEmcHits;

  // reconstruct pi0 and rho, in the variable, for example mrho_mn,
  // m is the index of pi0, n is the index of charged pi
  if (eventFlag == 1) { // pi0
    double delmpi0 = 999.;
    double mpi0 = 0.1349766;
    int g1 = -1;
    int g2 = -1;
    HepLorentzVector pTot;
    for (int i = 0; i < nGam - 1; i++) {
      for (int j = i + 1; j < nGam; j++) {
        HepLorentzVector p2g = pGam[i] + pGam[j];
        double m2g = p2g.m();
        double dltm2g = abs(m2g - mpi0);
        if (delmpi0 > dltm2g) {
          delmpi0 = dltm2g;
          g1 = i;
          g2 = j;
        }
      }
    }
    m_dlt_mpi0 = delmpi0;
    m_eg_11 = pGam[g1].e();
    m_eg_12 = pGam[g2].e();
    m_eg_21 = -1.0;
    m_eg_22 = -1.0;
    m_mpi0_1 = (pGam[g1] + pGam[g2]).m();
    m_mpi0_2 = -1.0;
    m_mrho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_21 = -1.0;
    m_mrho_22 = -1.0;
    m_prho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_21 = -10.0;
    m_prho_22 = -10.0;

    // 2020.03.15, calculate the helicity angle
    HepLorentzVector ppi0_1 = pGam[g1] + pGam[g2];
    Hep3Vector bv_pi0_rec = ppi0_1.boostVector();
    HepLorentzVector p4rho_11 = P4_1[2] + pGam[g1] + pGam[g2];
    HepLorentzVector p4rho_12 = P4_2[2] + pGam[g1] + pGam[g2];
    Hep3Vector bv_rho_rec_11 = p4rho_11.boostVector();
    Hep3Vector bv_rho_rec_12 = p4rho_12.boostVector();

    // 2020.03.15, calculate acop angle
    Hep3Vector P3_rho_11 = p4rho_11.vect();
    Hep3Vector P3_rho_12 = p4rho_12.vect();
    m_erho_11_acop = P3_trk2[0].deltaPhi(P3_rho_11) * 180 / CLHEP::pi;
    m_erho_12_acop = P3_trk1[0].deltaPhi(P3_rho_12) * 180 / CLHEP::pi;
    m_murho_11_acop = P3_trk2[1].deltaPhi(P3_rho_11) * 180 / CLHEP::pi;
    m_murho_12_acop = P3_trk1[1].deltaPhi(P3_rho_12) * 180 / CLHEP::pi;

    HepLorentzVector P4_erho_11 = cms_P4_2[0] + p4rho_11;
    HepLorentzVector P4_erho_12 = cms_P4_1[0] + p4rho_12;
    HepLorentzVector P4_murho_11 = cms_P4_2[1] + p4rho_11;
    HepLorentzVector P4_murho_12 = cms_P4_1[1] + p4rho_12;

    m_miss_m2_erho_11 =
        (m_ecms - P4_erho_11.e()) * (m_ecms - P4_erho_11.e()) -
        ((p_all - P4_erho_11).rho()) * ((p_all - P4_erho_11).rho());
    m_miss_m2_erho_12 =
        (m_ecms - P4_erho_12.e()) * (m_ecms - P4_erho_12.e()) -
        ((p_all - P4_erho_12).rho()) * ((p_all - P4_erho_12).rho());
    m_miss_m2_murho_11 =
        (m_ecms - P4_murho_11.e()) * (m_ecms - P4_murho_11.e()) -
        ((p_all - P4_murho_11).rho()) * ((p_all - P4_murho_11).rho());
    m_miss_m2_murho_12 =
        (m_ecms - P4_murho_12.e()) * (m_ecms - P4_murho_12.e()) -
        ((p_all - P4_murho_12).rho()) * ((p_all - P4_murho_12).rho());

    m_erho_11_angle = P3_trk2[0].angle(P3_rho_11) * 180 / CLHEP::pi;
    m_erho_12_angle = P3_trk1[0].angle(P3_rho_12) * 180 / CLHEP::pi;
    m_murho_11_angle = P3_trk2[1].angle(P3_rho_11) * 180 / CLHEP::pi;
    m_murho_12_angle = P3_trk1[1].angle(P3_rho_12) * 180 / CLHEP::pi;

    // 2020.03.15, calculate theta_1
    P4_1[2].boost(-bv_rho_rec_11);
    P4_2[2].boost(-bv_rho_rec_12);

    m_theta_pi_11 = P4_1[2].theta();
    m_theta_pi_12 = P4_2[2].theta();
    m_theta_pi_21 = -10.;
    m_theta_pi_22 = -10.;

    m_ppi_cms_11 = P4_1[2].rho();
    m_ppi_cms_12 = P4_2[2].rho();
    m_ppi_cms_21 = -10.;
    m_ppi_cms_22 = -10.;

    // 2020.03.15, calculate theta_2
    pGam[g1].boost(-bv_pi0_rec);
    pGam[g2].boost(-bv_pi0_rec);
    ppi0_1.boost(-bv_rho_rec_11);
    if (Ncut[6] % 2 == 0) {
      m_theta_m_11 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_11 = pGam[g2].vect().angle(ppi0_1.vect());
    }

    // 2020.03.15, boost pi0 to the second rho
    ppi0_1.boost(bv_rho_rec_11);
    ppi0_1.boost(-bv_rho_rec_12);

    if (Ncut[6] % 2 == 0) {
      m_theta_m_12 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_12 = pGam[g2].vect().angle(ppi0_1.vect());
    }
    m_theta_m_21 = -10.;
    m_theta_m_22 = -10.;
  } // end of reconstruct pi0

  if (eventFlag == 2) { // rho
    double delmpi0 = 999.;
    double mpi0 = 0.1349766;
    int g1 = -1;
    int g2 = -1;
    int g3 = -1;
    int g4 = -1;
    HepLorentzVector pTot;
    for (int i = 0; i < nGam - 3; i++) {
      for (int j = i + 1; j < nGam - 2; j++) {
        for (int k = j + 1; k < nGam - 1; k++) {
          for (int m = k + 1; m < nGam; m++) {
            HepLorentzVector p2g_1 = pGam[i] + pGam[j];
            HepLorentzVector p2g_2 = pGam[k] + pGam[m];
            double m2g_1 = p2g_1.m();
            double m2g_2 = p2g_2.m();
            double dltm2g = sqrt((m2g_1 - mpi0) * (m2g_1 - mpi0) +
                                 (m2g_2 - mpi0) * (m2g_2 - mpi0));
            if (delmpi0 > dltm2g) {
              delmpi0 = dltm2g;
              g1 = i;
              g2 = j;
              g3 = k;
              g4 = m;
            }
            p2g_1 = pGam[i] + pGam[k];
            p2g_2 = pGam[j] + pGam[m];
            m2g_1 = p2g_1.m();
            m2g_2 = p2g_2.m();
            dltm2g = sqrt((m2g_1 - mpi0) * (m2g_1 - mpi0) +
                          (m2g_2 - mpi0) * (m2g_2 - mpi0));
            if (delmpi0 > dltm2g) {
              delmpi0 = dltm2g;
              g1 = i;
              g2 = k;
              g3 = j;
              g4 = m;
            }
            p2g_1 = pGam[i] + pGam[m];
            p2g_2 = pGam[j] + pGam[k];
            m2g_1 = p2g_1.m();
            m2g_2 = p2g_2.m();
            dltm2g = sqrt((m2g_1 - mpi0) * (m2g_1 - mpi0) +
                          (m2g_2 - mpi0) * (m2g_2 - mpi0));
            if (delmpi0 > dltm2g) {
              delmpi0 = dltm2g;
              g1 = i;
              g2 = m;
              g3 = j;
              g4 = k;
            }
          }
        }
      }
    }
    m_dlt_mpi0 = delmpi0;
    m_eg_11 = pGam[g1].e();
    m_eg_12 = pGam[g2].e();
    m_eg_21 = pGam[g3].e();
    m_eg_22 = pGam[g4].e();
    m_mpi0_1 = (pGam[g1] + pGam[g2]).m();
    m_mpi0_2 = (pGam[g3] + pGam[g4]).m();
    m_mrho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).m();
    m_mrho_21 = (P4_1[2] + pGam[g3] + pGam[g4]).m();
    m_mrho_22 = (P4_2[2] + pGam[g3] + pGam[g4]).m();
    m_prho_11 = (P4_1[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_12 = (P4_2[2] + pGam[g1] + pGam[g2]).rho();
    m_prho_21 = (P4_1[2] + pGam[g3] + pGam[g4]).rho();
    m_prho_22 = (P4_2[2] + pGam[g3] + pGam[g4]).rho();

    // 2020.03.15, calculate the helicity angle
    // the first pi0
    HepLorentzVector ppi0_1 = pGam[g1] + pGam[g2];
    Hep3Vector bv_pi0_1_rec = ppi0_1.boostVector();
    HepLorentzVector p4rho_11 = P4_1[2] + pGam[g1] + pGam[g2];
    HepLorentzVector p4rho_12 = P4_2[2] + pGam[g1] + pGam[g2];
    Hep3Vector bv_rho_rec_11 = p4rho_11.boostVector();
    Hep3Vector bv_rho_rec_12 = p4rho_12.boostVector();

    // 2020.03.15, calculate acop angle
    Hep3Vector P3_rho_11 = p4rho_11.vect();
    Hep3Vector P3_rho_12 = p4rho_12.vect();
    // 2020.03.15, calculate theta_1
    P4_1[2].boost(-bv_rho_rec_11);
    P4_2[2].boost(-bv_rho_rec_12);

    m_theta_pi_11 = P4_1[2].theta();
    m_theta_pi_12 = P4_2[2].theta();

    m_ppi_cms_11 = P4_1[2].rho();
    m_ppi_cms_12 = P4_2[2].rho();

    // 2020.03.15, boost back to lab
    P4_1[2].boost(bv_rho_rec_11);
    P4_2[2].boost(bv_rho_rec_12);

    // 2020.03.15, the second pi0
    HepLorentzVector ppi0_2 = pGam[g3] + pGam[g4];
    Hep3Vector bv_pi0_2_rec = ppi0_2.boostVector();
    HepLorentzVector p4rho_21 = P4_1[2] + pGam[g3] + pGam[g4];
    HepLorentzVector p4rho_22 = P4_2[2] + pGam[g3] + pGam[g4];
    Hep3Vector bv_rho_rec_21 = p4rho_21.boostVector();
    Hep3Vector bv_rho_rec_22 = p4rho_22.boostVector();

    // 2020.03.15, calculate acop angle
    Hep3Vector P3_rho_21 = p4rho_21.vect();
    Hep3Vector P3_rho_22 = p4rho_22.vect();
    m_rhorho_11_acop = P3_rho_11.deltaPhi(P3_rho_22) * 180 / CLHEP::pi;
    m_rhorho_12_acop = P3_rho_12.deltaPhi(P3_rho_21) * 180 / CLHEP::pi;

    HepLorentzVector P4_rhorho_11 = p4rho_11 + p4rho_22;
    HepLorentzVector P4_rhorho_12 = p4rho_12 + p4rho_21;

    m_miss_m2_rhorho_11 =
        (m_ecms - P4_rhorho_11.e()) * (m_ecms - P4_rhorho_11.e()) -
        ((p_all - P4_rhorho_11).rho()) * ((p_all - P4_rhorho_11).rho());
    m_miss_m2_rhorho_12 =
        (m_ecms - P4_rhorho_12.e()) * (m_ecms - P4_rhorho_12.e()) -
        ((p_all - P4_rhorho_12).rho()) * ((p_all - P4_rhorho_12).rho());

    m_rhorho_11_angle = P3_rho_11.angle(P3_rho_22) * 180 / CLHEP::pi;
    m_rhorho_12_angle = P3_rho_12.angle(P3_rho_21) * 180 / CLHEP::pi;

    // 2020.03.15, calculate theta_1
    P4_1[2].boost(-bv_rho_rec_21);
    P4_2[2].boost(-bv_rho_rec_22);

    m_theta_pi_21 = P4_1[2].theta();
    m_theta_pi_22 = P4_2[2].theta();

    m_ppi_cms_21 = P4_1[2].rho();
    m_ppi_cms_22 = P4_2[2].rho();

    // 2020.03.15, calculate theta_2
    // the first pi0
    pGam[g1].boost(-bv_pi0_1_rec);
    pGam[g2].boost(-bv_pi0_1_rec);
    ppi0_1.boost(-bv_rho_rec_11);
    if (Ncut[6] % 2 == 0) {
      m_theta_m_11 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_11 = pGam[g2].vect().angle(ppi0_1.vect());
    }

    // 2020.03.15, boost pi0 to the second rho
    ppi0_1.boost(bv_rho_rec_11);
    ppi0_1.boost(-bv_rho_rec_12);

    if (Ncut[6] % 2 == 0) {
      m_theta_m_12 = pGam[g1].vect().angle(ppi0_1.vect());
    } else {
      m_theta_m_12 = pGam[g2].vect().angle(ppi0_1.vect());
    }

    // the second pi0
    pGam[g3].boost(-bv_pi0_2_rec);
    pGam[g4].boost(-bv_pi0_2_rec);
    ppi0_2.boost(-bv_rho_rec_21);
    if (Ncut[6] % 2 == 0) {
      m_theta_m_21 = pGam[g3].vect().angle(ppi0_2.vect());
    } else {
      m_theta_m_21 = pGam[g4].vect().angle(ppi0_2.vect());
    }

    // 2020.03.15, boost pi0 to the second rho
    ppi0_2.boost(bv_rho_rec_21);
    ppi0_2.boost(-bv_rho_rec_22);

    if (Ncut[6] % 2 == 0) {
      m_theta_m_22 = pGam[g3].vect().angle(ppi0_2.vect());
    } else {
      m_theta_m_22 = pGam[g4].vect().angle(ppi0_2.vect());
    }
  }
  if (eventFlag != 1 && eventFlag != 2) {
    m_eg_11 = -10.;
    m_eg_12 = -10.;
    m_eg_21 = -10.;
    m_eg_22 = -10.;
    m_mpi0_1 = -10.;
    m_mpi0_2 = -10.;
    m_mrho_11 = -10.;
    m_mrho_12 = -10.;
    m_mrho_21 = -10.;
    m_mrho_22 = -10.;
    m_prho_11 = -10.;
    m_prho_12 = -10.;
    m_prho_21 = -10.;
    m_prho_22 = -10.;
    m_theta_pi_11 = -10.;
    m_theta_pi_12 = -10.;
    m_theta_pi_21 = -10.;
    m_theta_pi_22 = -10.;
    m_theta_m_11 = -10.;
    m_theta_m_12 = -10.;
    m_theta_m_21 = -10.;
    m_theta_m_22 = -10.;
  }

  m_tuple3->write();
  return StatusCode::SUCCESS;
}

// Thanks for Zhou Xingyu's help 2020.03.28
bool tautauAlg::mctruth() {
  SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),
                                                   "/Event/MC/McParticleCol");
  if (!mcParticleCol) {
    std::cout << "Could not retrieve McParticelCol" << std::endl;
    return (StatusCode::FAILURE);
  }

  /*
      Event::McParticleCol::iterator iter_mc=mcParticleCol->begin();
      m_idxmc=0;
      for (;iter_mc!=mcParticleCol->end();iter_mc++){
          //if((*iter_mc)->primaryParticle()) continue;
          if(!(*iter_mc)->decayFromGenerator()) continue;
          //std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdgid="<<(*iter_mc)->particleProperty()<<"\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
          m_pdgid[m_idxmc]=(*iter_mc)->particleProperty();
          m_motheridx[m_idxmc]=((*iter_mc)->mother()).trackIndex();
          m_idxmc++;
      }

  */
  /*
      Event::McParticleCol::iterator iter_mc=mcParticleCol->begin();
      int pdgid=(*iter_mc)->particleProperty();
      unsigned int idxmc;
      unsigned int motheridx;
      m_idxmc=0;
      if(pdgid==90022||pdgid==80022){
          // The first "iter_mc++" in the head of the following loop is added to
     remove the unnecessary virtual photons, so as to unify, for example, the
     following three types of processes: 1) e+ e- --> A B, 2) e+ e- --> 80022
     --> A B, and 3) e+ e- --> 90022 --> A B.
          for(iter_mc++;iter_mc!=mcParticleCol->end();iter_mc++){
              if(!(*iter_mc)->decayFromGenerator()) continue;
              //
     std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdgid="<<(*iter_mc)->particleProperty()<<"\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
              pdgid=(*iter_mc)->particleProperty();
              idxmc=(*iter_mc)->trackIndex();
              motheridx=((*iter_mc)->mother()).trackIndex();
              m_pdgid[m_idxmc]=pdgid;
              if(idxmc==motheridx||motheridx==0) m_motheridx[m_idxmc]=idxmc-1;
              else m_motheridx[m_idxmc]=motheridx-1;
              m_idxmc++;
          }
      }
      else{
          for(;iter_mc!=mcParticleCol->end();iter_mc++){
              if(!(*iter_mc)->decayFromGenerator()) continue;
              //
     std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdgid="<<(*iter_mc)->particleProperty()<<"\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
              m_pdgid[m_idxmc]=(*iter_mc)->particleProperty();
              m_motheridx[m_idxmc]=((*iter_mc)->mother()).trackIndex();
              m_idxmc++;
          }
      }
  */

  // -------------------------------DDbar-----------------------------------

  Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
  const int incpdgid0 = 9030443; // 9030443 is the PDG code of Y(4260).
  m_idxmc = 0;
  m_rec_idxmc = 0;
  int ngamma = 0;
  bool incPdcy(false);
  int rootIndex(-1);
  for (; iter_mc != mcParticleCol->end(); iter_mc++) {
    // std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdgid="<<(*iter_mc)->particleProperty()<<"\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
    // if((*iter_mc)->primaryParticle()) continue;
    if (!(*iter_mc)->decayFromGenerator())
      continue;
    if (incPdcy == false && (*iter_mc)->particleProperty() == -22) {
      m_pdgid[m_idxmc] = 222222222; // 222222222 is the default PDG code of ISR
                                    // photons in TopoAna.
      m_rec_pdgid[m_idxmc] = 22222222;
      m_motheridx[m_idxmc] = -1;
      m_rec_motheridx[m_idxmc] = -1;
      m_idxmc++;
      m_rec_idxmc++;
      ngamma++;
      continue;
    } else if ((*iter_mc)->particleProperty() == incpdgid0) {
      incPdcy = true;
      rootIndex = (*iter_mc)->trackIndex();
    } else if (incPdcy == true && (*iter_mc)->particleProperty() == 11 &&
               ((*iter_mc)->mother()).trackIndex() ==
                   (*iter_mc)->trackIndex()) {
      continue; // This statement is used to remove the redundant electron
                // following Y(4260).
    }
    if (!incPdcy)
      continue;
    m_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
    m_rec_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
    if ((*iter_mc)->particleProperty() == incpdgid0) {
      m_motheridx[m_idxmc] = -1;
      m_rec_motheridx[m_idxmc] = -1;
    } else if (((*iter_mc)->mother()).particleProperty() == incpdgid0) {
      m_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex + ngamma;
      m_rec_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex + ngamma;
    } else {
      m_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
      m_rec_motheridx[m_idxmc] =
          ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
    }
    m_idxmc++;
    m_rec_idxmc++;
  }

  // ------------------------------(hadrons,BesTwogamma)----------------------------
  if (incPdcy == false) {
    iter_mc = mcParticleCol->begin();
    const int incpdgid1 = 91; // 91 is the PDG code of cluster
    const int incpdgid2 = 92; // 92 is the PDG code of string
    m_idxmc = 0;
    m_rec_idxmc = 0;
    ngamma = 0;
    incPdcy = false;
    int rootIndex(-1);
    bool findTFRE = false; // TFRE is short for is The First Redundant Electron.
    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      // std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdgid="<<(*iter_mc)->particleProperty()<<"\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
      // if((*iter_mc)->primaryParticle()) continue;
      if (!(*iter_mc)->decayFromGenerator())
        continue;
      if (incPdcy == false && (*iter_mc)->particleProperty() == 22) {
        m_pdgid[m_idxmc] = 222222222; // 222222222 is the default PDG code of
                                      // ISR photons in TopoAna.
        m_rec_pdgid[m_idxmc] = 22222222;
        m_motheridx[m_idxmc] = -1;
        m_rec_motheridx[m_idxmc] = -1;
        m_idxmc++;
        m_rec_idxmc++;
        ngamma++;
        continue;
      } else if ((*iter_mc)->particleProperty() == incpdgid1 ||
                 (*iter_mc)->particleProperty() == incpdgid2) {
        incPdcy = true;
        rootIndex = (*iter_mc)->trackIndex();
      } else if (incPdcy == true && (*iter_mc)->particleProperty() == 11 &&
                 ((*iter_mc)->mother()).trackIndex() ==
                     (*iter_mc)->trackIndex()) {
        if (findTFRE == false) {
          m_pdgid[m_idxmc] = 11;
          m_rec_pdgid[m_idxmc] = 11;
          findTFRE = true;
        } else {
          m_pdgid[m_idxmc] = -11;
          m_rec_pdgid[m_idxmc] = -11;
        }
        m_motheridx[m_idxmc] = -1;
        m_rec_motheridx[m_idxmc] = -1;
        m_idxmc++;
        m_rec_idxmc++;
        continue;
      }
      if (!incPdcy || (*iter_mc)->particleProperty() == incpdgid1 ||
          (*iter_mc)->particleProperty() == incpdgid2)
        continue;
      m_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      m_rec_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      if (((*iter_mc)->mother()).particleProperty() == incpdgid1 ||
          ((*iter_mc)->mother()).particleProperty() == incpdgid2) {
        m_motheridx[m_idxmc] = -1;
        m_rec_motheridx[m_idxmc] = -1;
      } else {
        m_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
        m_rec_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - rootIndex - 1 + ngamma;
      }
      m_idxmc++;
      m_rec_idxmc++;
    }
  }

  // --------------------------QED (bhabha, digamma, dimu, ditau,
  // ISR)------------------
  if (incPdcy == false) {
    iter_mc = mcParticleCol->begin();
    int pdgid = (*iter_mc)->particleProperty();
    m_idxmc = 0;
    m_rec_idxmc = 0;
    int idxOfTRE = -1;
    int nOfTRQG = 0;
    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      // if(!(*iter_mc)->decayFromGenerator()) continue;
      // std::cout<<"idxmc="<<m_idxmc<<"\t"<<"pdgid="<<(*iter_mc)->particleProperty()<<"\t"<<"motheridx="<<((*iter_mc)->mother()).trackIndex()<<std::endl;
      if (pdgid == 23 && (*iter_mc)->particleProperty() == 11 &&
          ((*iter_mc)->mother()).trackIndex() == m_idxmc) {
        idxOfTRE = m_idxmc;
        continue;
      }
      if (abs((*iter_mc)->particleProperty()) == 1 ||
          abs((*iter_mc)->particleProperty()) == 2 ||
          abs((*iter_mc)->particleProperty()) == 3 ||
          abs((*iter_mc)->particleProperty()) == 4 ||
          abs((*iter_mc)->particleProperty()) == 21) {
        nOfTRQG++;
        continue;
      }
      m_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      m_rec_pdgid[m_idxmc] = (*iter_mc)->particleProperty();
      if (idxOfTRE != -1 && ((*iter_mc)->mother()).trackIndex() > idxOfTRE) {
        m_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - 1 - nOfTRQG;
        m_rec_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - 1 - nOfTRQG;
      } else {
        m_motheridx[m_idxmc] = ((*iter_mc)->mother()).trackIndex() - nOfTRQG;
        m_rec_motheridx[m_idxmc] =
            ((*iter_mc)->mother()).trackIndex() - nOfTRQG;
      }
      m_idxmc++;
      m_rec_idxmc++;
    }
  }
  if (!mcParticleCol) {
    return false;
  } else {
    Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
    for (; iter_mc != mcParticleCol->end(); iter_mc++) {
      if (!(*iter_mc)->decayFromGenerator())
        continue;
      HepLorentzVector p_mc = (*iter_mc)->initialFourMomentum();
      int self_id = (*iter_mc)->particleProperty();
      int mother_id = ((*iter_mc)->mother()).particleProperty();
      int grandmother_id = (((*iter_mc)->mother()).mother()).particleProperty();
      if (self_id == 11 && mother_id == 15) { // tau^- -> e^- nu anti-nu
        setV4(m_info_tru, p_mc, 1);
      }
      if (self_id == -11 && mother_id == -15) { // tau^+ -> e^+ nu anti-nu
        setV4(m_info_tru, p_mc, 2);
      }
      if (self_id == -13 && mother_id == -15) { // tau^+ -> mu^+ nu anti-nu
        setV4(m_info_tru, p_mc, 3);
      }
      if (self_id == 13 && mother_id == 15) { // tau^- -> mu^- nu anti-nu
        setV4(m_info_tru, p_mc, 4);
      }
      if (self_id == 15) { // tau^-
        setV4(m_info_tru, p_mc, 5);
      }
      if (self_id == -15) { // tau^+
        setV4(m_info_tru, p_mc, 6);
      }
    }
  }

  return true;
}

bool tautauAlg::correctLeptonMomentumWithFSRPhoton(
    HepLorentzVector p4lepp, std::vector<HepLorentzVector> &p4FSRCor) {
  p4FSRCor.clear();

  SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),
                                        "/Event/EvtRec/EvtRecEvent");
  if (!evtRecEvent) {
    std::cout << "Could not find EvtRecEvent" << std::endl;
    return false;
  }
  SmartDataPtr<EvtRecTrackCol> evtRecTrackCol(eventSvc(),
                                              "/Event/EvtRec/EvtRecTrackCol");
  if (!evtRecTrackCol) {
    std::cout << "Could not find EvtRecTrackCol" << std::endl;
    return false;
  }

  HepLorentzVector p4FSRLepp(0., 0., 0., 0.);
  for (int i = evtRecEvent->totalCharged(); i < evtRecEvent->totalTracks();
       i++) {
    EvtRecTrackIterator itTrk = evtRecTrackCol->begin() + i;
    if (!(*itTrk)->isEmcShowerValid())
      continue;
    if (!isGoodShower(*itTrk))
      continue;
    RecEmcShower *emcTrk = (*itTrk)->emcShower();
    double eraw = emcTrk->energy();
    double phi = emcTrk->phi();
    double the = emcTrk->theta();
    HepLorentzVector p4shower(eraw * sin(the) * cos(phi),
                              eraw * sin(the) * sin(phi), eraw * cos(the),
                              eraw);
    double dthelepp = 200.;
    dthelepp = p4lepp.vect().theta(p4shower.vect());
    dthelepp =
        fabs(fmod(dthelepp + CLHEP::twopi + CLHEP::twopi + pi, CLHEP::twopi) -
             CLHEP::pi);
    if (dthelepp < CLHEP::pi * 5. / 180.) {
      p4FSRLepp += p4shower;
      m_eleppFSRArray[m_nleppFSR] = p4shower.e();
      m_dtheleppFSRArray[m_nleppFSR] = dthelepp;
      m_nleppFSR++;
    }
  }
  p4FSRCor.push_back(p4FSRLepp);
  return true;
}

//***************************Output******************************
StatusCode tautauAlg::finalize() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in finalize()" << endmsg;
  cout << "Ntot=                    " << Ncut[0] << endl;
  cout << "nGoodTrack==2            " << Ncut[1] << endl;
  cout << "nGoodGam>=" << m_gamNumCut << "              " << Ncut[2] << endl;
  cout << "nGoodPi0<=" << m_pi0NumCut << "              " << Ncut[3] << endl;
  cout << "nGoodTrack_Emc==2        " << Ncut[4] << endl;
  cout << "npi_pid==1&&nLep_pid==1  " << Ncut[5] << endl;
  // cout << "differet_id=        " << Ncut[5] << endl;
  // cout << "opposite_charge=    " << Ncut[6] << endl;
  return StatusCode::SUCCESS;
}
