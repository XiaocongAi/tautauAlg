#ifndef Physics_Analysis_tautauAlg_H
#define Physics_Analysis_tautauAlg_H
//#include <string>

#include "GaudiKernel/ITHistSvc.h"
#include "tautauAlg/util/MyConst.h"
#include <vector>
#include "TH1.h"
#include "TH1F.h"
using namespace std;

class tautauAlg : public Algorithm {
public:
  tautauAlg(const std::string &name, ISvcLocator *pSvcLocator);
  StatusCode initialize();
  StatusCode execute();
  StatusCode finalize();
  bool mctruth();
  bool
  correctLeptonMomentumWithFSRPhoton(HepLorentzVector p4lepp,
                                     std::vector<HepLorentzVector> &p4FSRCor);
  void setPxPyPzE(HepLorentzVector &p4, double px, double py, double pz,
                  double e) const;

private:
  // event Num
  int m_nEvtDisp;
  int m_nEvt;

  bool m_writeGenOnly;
  bool m_mctruth;
  bool m_appFSRCorrection;
  double m_ecms;
  int m_testMC;

  NTuple::Tuple *m_tuple3;

  NTuple::Item<int> m_rec_run;
  NTuple::Item<int> m_rec_event;
  NTuple::Item<int> m_rec_nchrg;
  NTuple::Item<int> m_rec_nneu;
  NTuple::Array<double> m_gamcon;
  NTuple::Array<double> m_eta_gamee;
  NTuple::Array<double> m_pi0_gamee;
  NTuple::Array<double> m_elp_info;
  NTuple::Array<double> m_mum_info;
  NTuple::Array<double> m_elp_prob;
  NTuple::Array<double> m_mum_prob;
  NTuple::Item<double> m_esum;
  NTuple::Item<double> m_egam_max;
  NTuple::Item<int> m_id_gam_max;

  NTuple::Item<int> m_rec_idxmc;
  NTuple::Array<int> m_rec_pdgid;
  NTuple::Array<int> m_rec_motheridx;

  NTuple::Item<int> m_nsig_shower;
  NTuple::Matrix<double> m_shower_info;

  NTuple::Item<int> m_nleppFSR;
  NTuple::Array<double> m_eleppFSRArray;
  NTuple::Array<double> m_dtheleppFSRArray;

  NTuple::Array<double> m_gam_fromelp;
  NTuple::Array<double> m_elp_corfsr;
  NTuple::Array<double> m_info_elp_mum;
  NTuple::Matrix<double> m_info_emu;

  NTuple::Item<int> m_charge_1;
  NTuple::Item<int> m_charge_2;
  NTuple::Item<double> m_pmdc_1;
  NTuple::Item<double> m_pmdc_2;
  NTuple::Item<double> m_ptmdc_1;
  NTuple::Item<double> m_ptmdc_2;

  NTuple::Array<double> m_pkal_1;
  NTuple::Array<double> m_pkal_2;
  NTuple::Array<double> m_ptkal_1;
  NTuple::Array<double> m_ptkal_2;
  NTuple::Array<double> m_p3_1;
  NTuple::Array<double> m_p3_2;

  NTuple::Item<double> m_probe_1;
  NTuple::Item<double> m_probe_2;
  NTuple::Item<double> m_probmu_1;
  NTuple::Item<double> m_probmu_2;
  NTuple::Item<double> m_probpi_1;
  NTuple::Item<double> m_probpi_2;
  NTuple::Item<double> m_probk_1;
  NTuple::Item<double> m_probk_2;
  NTuple::Item<double> m_probp_1;
  NTuple::Item<double> m_probp_2;

  NTuple::Item<int> m_nTofInfo_1;
  NTuple::Item<int> m_nTofInfo_2;
  NTuple::Item<int> m_fgtof_1;
  NTuple::Item<double> m_counter_1;
  NTuple::Item<double> m_isbarrel_1;
  NTuple::Item<double> m_layertof_1;
  NTuple::Item<double> m_iscluster_1;
  NTuple::Item<double> m_tof_1;
  NTuple::Item<double> m_texe_1;
  NTuple::Item<double> m_texmu_1;
  NTuple::Item<double> m_texpi_1;
  NTuple::Item<double> m_texk_1;
  NTuple::Item<double> m_texp_1;
  NTuple::Item<double> m_dte_1;
  NTuple::Item<double> m_dtmu_1;
  NTuple::Item<double> m_dtpi_1;
  NTuple::Item<double> m_dtk_1;
  NTuple::Item<double> m_dtp_1;
  NTuple::Item<double> m_evp_1;
  NTuple::Item<double> m_ene_1;
  NTuple::Item<int> m_nTrk_EMC;

  NTuple::Item<int> m_fgtof_2;
  NTuple::Item<double> m_counter_2;
  NTuple::Item<double> m_isbarrel_2;
  NTuple::Item<double> m_layertof_2;
  NTuple::Item<double> m_iscluster_2;
  NTuple::Item<double> m_tof_2;
  NTuple::Item<double> m_texe_2;
  NTuple::Item<double> m_texmu_2;
  NTuple::Item<double> m_texpi_2;
  NTuple::Item<double> m_texk_2;
  NTuple::Item<double> m_texp_2;
  NTuple::Item<double> m_dte_2;
  NTuple::Item<double> m_dtmu_2;
  NTuple::Item<double> m_dtpi_2;
  NTuple::Item<double> m_dtk_2;
  NTuple::Item<double> m_dtp_2;
  NTuple::Item<double> m_evp_2;
  NTuple::Item<double> m_ene_2;

  NTuple::Item<double> m_theta_1;
  NTuple::Item<double> m_chi_e_1;
  NTuple::Item<double> m_chi_mu_1;
  NTuple::Item<double> m_chi_pi_1;
  NTuple::Item<double> m_chi_k_1;
  NTuple::Item<double> m_chi_p_1;

  NTuple::Item<double> m_theta_2;
  NTuple::Item<double> m_chi_e_2;
  NTuple::Item<double> m_chi_mu_2;
  NTuple::Item<double> m_chi_pi_2;
  NTuple::Item<double> m_chi_k_2;
  NTuple::Item<double> m_chi_p_2;

  NTuple::Item<long> m_maxhitsinlay_1;
  NTuple::Item<long> m_numhits_1;
  NTuple::Item<long> m_numlayers_1;
  NTuple::Item<double> m_depth_1;
  NTuple::Item<double> m_mucchi2_1;
  NTuple::Item<long> m_maxhitsinlay_2;
  NTuple::Item<long> m_numhits_2;
  NTuple::Item<long> m_numlayers_2;
  NTuple::Item<double> m_depth_2;
  NTuple::Item<double> m_mucchi2_2;

  NTuple::Array<float> m_pkal_tot_px;
  NTuple::Array<float> m_pkal_tot_py;
  NTuple::Array<float> m_pkal_tot_pz;
  NTuple::Array<float> m_pkal_tot_e;
  NTuple::Array<float> m_mTrks;
  NTuple::Array<float> m_m2Trks;

  NTuple::Item<double> m_ang_mdc_acol;
  NTuple::Item<double> m_ang_mdc_acop;
  NTuple::Item<double> m_the_add;
  NTuple::Item<double> m_phi_diff;
  NTuple::Item<double> m_PTEM;
  NTuple::Item<double> m_ee_acop;
  NTuple::Item<double> m_emu_acop;
  NTuple::Item<double> m_mue_acop;
  NTuple::Item<double> m_epi_acop;
  NTuple::Item<double> m_pie_acop;
  NTuple::Item<double> m_ek_acop;
  NTuple::Item<double> m_ke_acop;
  NTuple::Item<double> m_mumu_acop;
  NTuple::Item<double> m_mupi_acop;
  NTuple::Item<double> m_pimu_acop;
  NTuple::Item<double> m_pipi_acop;
  NTuple::Item<double> m_muk_acop;
  NTuple::Item<double> m_kmu_acop;
  NTuple::Item<double> m_pik_acop;
  NTuple::Item<double> m_kpi_acop;
  NTuple::Item<double> m_kk_acop;
  NTuple::Item<double> m_erho_11_acop;
  NTuple::Item<double> m_erho_12_acop;
  NTuple::Item<double> m_murho_11_acop;
  NTuple::Item<double> m_murho_12_acop;
  NTuple::Item<double> m_rhorho_11_acop;
  NTuple::Item<double> m_rhorho_12_acop;

  NTuple::Item<double> m_ee_angle;
  NTuple::Item<double> m_emu_angle;
  NTuple::Item<double> m_mue_angle;
  NTuple::Item<double> m_epi_angle;
  NTuple::Item<double> m_pie_angle;
  NTuple::Item<double> m_ek_angle;
  NTuple::Item<double> m_ke_angle;
  NTuple::Item<double> m_mumu_angle;
  NTuple::Item<double> m_mupi_angle;
  NTuple::Item<double> m_pimu_angle;
  NTuple::Item<double> m_pipi_angle;
  NTuple::Item<double> m_muk_angle;
  NTuple::Item<double> m_kmu_angle;
  NTuple::Item<double> m_pik_angle;
  NTuple::Item<double> m_kpi_angle;
  NTuple::Item<double> m_kk_angle;
  NTuple::Item<double> m_erho_11_angle;
  NTuple::Item<double> m_erho_12_angle;
  NTuple::Item<double> m_murho_11_angle;
  NTuple::Item<double> m_murho_12_angle;
  NTuple::Item<double> m_rhorho_11_angle;
  NTuple::Item<double> m_rhorho_12_angle;

  NTuple::Item<double> m_miss_m2_ee;
  NTuple::Item<double> m_miss_m2_emu;
  NTuple::Item<double> m_miss_m2_mue;
  NTuple::Item<double> m_miss_m2_epi;
  NTuple::Item<double> m_miss_m2_pie;
  NTuple::Item<double> m_miss_m2_ek;
  NTuple::Item<double> m_miss_m2_ke;
  NTuple::Item<double> m_miss_m2_mumu;
  NTuple::Item<double> m_miss_m2_mupi;
  NTuple::Item<double> m_miss_m2_pimu;
  NTuple::Item<double> m_miss_m2_pipi;
  NTuple::Item<double> m_miss_m2_muk;
  NTuple::Item<double> m_miss_m2_kmu;
  NTuple::Item<double> m_miss_m2_pik;
  NTuple::Item<double> m_miss_m2_kpi;
  NTuple::Item<double> m_miss_m2_kk;
  NTuple::Item<double> m_miss_m2_erho_11;
  NTuple::Item<double> m_miss_m2_erho_12;
  NTuple::Item<double> m_miss_m2_murho_11;
  NTuple::Item<double> m_miss_m2_murho_12;
  NTuple::Item<double> m_miss_m2_rhorho_11;
  NTuple::Item<double> m_miss_m2_rhorho_12;

  NTuple::Item<double> m_ee_PTEM;
  NTuple::Item<double> m_emu_PTEM;
  NTuple::Item<double> m_mue_PTEM;
  NTuple::Item<double> m_epi_PTEM;
  NTuple::Item<double> m_pie_PTEM;
  NTuple::Item<double> m_ek_PTEM;
  NTuple::Item<double> m_ke_PTEM;
  NTuple::Item<double> m_mumu_PTEM;
  NTuple::Item<double> m_mupi_PTEM;
  NTuple::Item<double> m_pimu_PTEM;
  NTuple::Item<double> m_pipi_PTEM;
  NTuple::Item<double> m_muk_PTEM;
  NTuple::Item<double> m_kmu_PTEM;
  NTuple::Item<double> m_pik_PTEM;
  NTuple::Item<double> m_kpi_PTEM;
  NTuple::Item<double> m_kk_PTEM;
  NTuple::Item<double> m_dlt_mpi0;
  NTuple::Item<double> m_eg_11;
  NTuple::Item<double> m_eg_12;
  NTuple::Item<double> m_eg_21;
  NTuple::Item<double> m_eg_22;
  NTuple::Item<double> m_mpi0_1;
  NTuple::Item<double> m_mpi0_2;
  NTuple::Item<double> m_mrho_11;
  NTuple::Item<double> m_mrho_12;
  NTuple::Item<double> m_mrho_21;
  NTuple::Item<double> m_mrho_22;
  NTuple::Item<double> m_prho_11;
  NTuple::Item<double> m_prho_12;
  NTuple::Item<double> m_prho_21;
  NTuple::Item<double> m_prho_22;
  NTuple::Item<double> m_theta_pi_11;
  NTuple::Item<double> m_theta_pi_12;
  NTuple::Item<double> m_theta_pi_21;
  NTuple::Item<double> m_theta_pi_22;
  NTuple::Item<double> m_ppi_cms_11;
  NTuple::Item<double> m_ppi_cms_12;
  NTuple::Item<double> m_ppi_cms_21;
  NTuple::Item<double> m_ppi_cms_22;
  NTuple::Item<double> m_theta_m_11;
  NTuple::Item<double> m_theta_m_12;
  NTuple::Item<double> m_theta_m_21;
  NTuple::Item<double> m_theta_m_22;

  // 2020.03.27 for topology
  NTuple::Item<int> m_nTrack;
  NTuple::Array<int> ma_trackID;
  NTuple::Array<int> ma_trackIndex;
  NTuple::Array<int> ma_motherID;
  NTuple::Array<int> ma_motherIndex;
  NTuple::Array<int> ma_fromGenerator;
  NTuple::Array<int> ma_primaryParticle;
  NTuple::Array<int> ma_px;
  NTuple::Array<int> ma_py;
  NTuple::Array<int> ma_pz;
  NTuple::Array<int> ma_e;

  NTuple::Tuple *m_tuple2;
  NTuple::Item<int> m_run;
  NTuple::Item<int> m_event;
  NTuple::Item<int> m_idxmc;
  NTuple::Array<int> m_pdgid;
  NTuple::Array<int> m_motheridx;
  NTuple::Matrix<double> m_info_tru;

  NTuple::Tuple *m_tuple1;
  NTuple::Item<int> m_within_kine_bounds;
  NTuple::Item<double> m_pi_px;
  NTuple::Item<double> m_pi_py;
  NTuple::Item<double> m_pi_pz;
  NTuple::Item<double> m_pi_e;

  NTuple::Item<double> m_tau_px;
  NTuple::Item<double> m_tau_py;
  NTuple::Item<double> m_tau_pz;
  NTuple::Item<double> m_tau_e;

  NTuple::Item<double> m_gam_px;
  NTuple::Item<double> m_gam_py;
  NTuple::Item<double> m_gam_pz;
  NTuple::Item<double> m_gam_e;

  NTuple::Item<double> m_nu_px;
  NTuple::Item<double> m_nu_py;
  NTuple::Item<double> m_nu_pz;
  NTuple::Item<double> m_nu_e;

  NTuple::Item<double> m_x;
  NTuple::Item<double> m_y;
  NTuple::Item<double> m_t;
  NTuple::Item<double> m_s;
  NTuple::Item<double> m_E;
  NTuple::Item<double> m_z;

  NTuple::Item<double> m_f_IB;
  NTuple::Item<double> m_f_VV;
  NTuple::Item<double> m_f_AV;
  NTuple::Item<double> m_f_AA;
  NTuple::Item<double> m_f_IB_V;
  NTuple::Item<double> m_f_IB_A;
  NTuple::Item<double> m_Gam_IB;
  NTuple::Item<double> m_Gam_SD;
  NTuple::Item<double> m_Gam_int;
  NTuple::Item<int> m_ngamfromTau;
  NTuple::Item<int> m_nefrommu;
  NTuple::Item<int> m_nefromtau;
  NTuple::Array<int> m_gam_motherIndex;
  NTuple::Item<int> m_nGam;

  NTuple::Item<int> m_npi0;
  NTuple::Array<double> m_m_pi0before;
  NTuple::Array<double> m_m_pi0fit;
  NTuple::Array<double> m_chi_pi0fit;
  NTuple::Array<int> m_i1gam_pi0;
  NTuple::Array<int> m_i2gam_pi0;

  NTuple::Array<double> m_eraw;
  NTuple::Array<double> m_phi;
  NTuple::Array<double> m_the;

  NTuple::Array<double> m_dang_pos_e;
  NTuple::Array<double> m_dthe_pos_e;
  NTuple::Array<double> m_dphi_pos_e;
  NTuple::Array<double> m_dang_neg_e;
  NTuple::Array<double> m_dthe_neg_e;
  NTuple::Array<double> m_dphi_neg_e;
  NTuple::Array<double> m_dang_pos_mu;
  NTuple::Array<double> m_dthe_pos_mu;
  NTuple::Array<double> m_dphi_pos_mu;
  NTuple::Array<double> m_dang_neg_mu;
  NTuple::Array<double> m_dthe_neg_mu;
  NTuple::Array<double> m_dphi_neg_mu;
  NTuple::Array<double> m_dang_pos_pi;
  NTuple::Array<double> m_dthe_pos_pi;
  NTuple::Array<double> m_dphi_pos_pi;
  NTuple::Array<double> m_dang_neg_pi;
  NTuple::Array<double> m_dthe_neg_pi;
  NTuple::Array<double> m_dphi_neg_pi;
  NTuple::Array<double> m_dang_pos_k;
  NTuple::Array<double> m_dthe_pos_k;
  NTuple::Array<double> m_dphi_pos_k;
  NTuple::Array<double> m_dang_neg_k;
  NTuple::Array<double> m_dthe_neg_k;
  NTuple::Array<double> m_dphi_neg_k;

  NTuple::Item<double> m_px_pos_e;
  NTuple::Item<double> m_py_pos_e;
  NTuple::Item<double> m_pz_pos_e;
  NTuple::Item<double> m_px_pos_mu;
  NTuple::Item<double> m_py_pos_mu;
  NTuple::Item<double> m_pz_pos_mu;
  NTuple::Item<double> m_px_pos_pi;
  NTuple::Item<double> m_py_pos_pi;
  NTuple::Item<double> m_pz_pos_pi;
  NTuple::Item<double> m_px_pos_k;
  NTuple::Item<double> m_py_pos_k;
  NTuple::Item<double> m_pz_pos_k;
  NTuple::Item<double> m_p_pos;
  NTuple::Item<double> m_pt_pos;

  NTuple::Item<double> m_px_neg_e;
  NTuple::Item<double> m_py_neg_e;
  NTuple::Item<double> m_pz_neg_e;
  NTuple::Item<double> m_px_neg_mu;
  NTuple::Item<double> m_py_neg_mu;
  NTuple::Item<double> m_pz_neg_mu;
  NTuple::Item<double> m_px_neg_pi;
  NTuple::Item<double> m_py_neg_pi;
  NTuple::Item<double> m_pz_neg_pi;
  NTuple::Item<double> m_px_neg_k;
  NTuple::Item<double> m_py_neg_k;
  NTuple::Item<double> m_pz_neg_k;
  NTuple::Item<double> m_p_neg;
  NTuple::Item<double> m_pt_neg;

  NTuple::Item<double> m_cang_ee;
  NTuple::Item<double> m_cthe_ee;
  NTuple::Item<double> m_cphi_ee;
  NTuple::Item<double> m_costheta_ee;
  NTuple::Item<double> m_cang_emu;
  NTuple::Item<double> m_cthe_emu;
  NTuple::Item<double> m_cphi_emu;
  NTuple::Item<double> m_costheta_emu;
  NTuple::Item<double> m_cang_epi;
  NTuple::Item<double> m_cthe_epi;
  NTuple::Item<double> m_cphi_epi;
  NTuple::Item<double> m_costheta_epi;
  NTuple::Item<double> m_cang_ek;
  NTuple::Item<double> m_cthe_ek;
  NTuple::Item<double> m_cphi_ek;
  NTuple::Item<double> m_costheta_ek;
  NTuple::Item<double> m_cang_mue;
  NTuple::Item<double> m_cthe_mue;
  NTuple::Item<double> m_cphi_mue;
  NTuple::Item<double> m_costheta_mue;
  NTuple::Item<double> m_cang_mumu;
  NTuple::Item<double> m_cthe_mumu;
  NTuple::Item<double> m_cphi_mumu;
  NTuple::Item<double> m_costheta_mumu;
  NTuple::Item<double> m_cang_mupi;
  NTuple::Item<double> m_cthe_mupi;
  NTuple::Item<double> m_cphi_mupi;
  NTuple::Item<double> m_costheta_mupi;
  NTuple::Item<double> m_cang_muk;
  NTuple::Item<double> m_cthe_muk;
  NTuple::Item<double> m_cphi_muk;
  NTuple::Item<double> m_costheta_muk;
  NTuple::Item<double> m_cang_pie;
  NTuple::Item<double> m_cthe_pie;
  NTuple::Item<double> m_cphi_pie;
  NTuple::Item<double> m_costheta_pie;
  NTuple::Item<double> m_cang_pimu;
  NTuple::Item<double> m_cthe_pimu;
  NTuple::Item<double> m_cphi_pimu;
  NTuple::Item<double> m_costheta_pimu;
  NTuple::Item<double> m_cang_pipi;
  NTuple::Item<double> m_cthe_pipi;
  NTuple::Item<double> m_cphi_pipi;
  NTuple::Item<double> m_costheta_pipi;
  NTuple::Item<double> m_cang_pik;
  NTuple::Item<double> m_cthe_pik;
  NTuple::Item<double> m_cphi_pik;
  NTuple::Item<double> m_costheta_pik;
  NTuple::Item<double> m_cang_ke;
  NTuple::Item<double> m_cthe_ke;
  NTuple::Item<double> m_cphi_ke;
  NTuple::Item<double> m_costheta_ke;
  NTuple::Item<double> m_cang_kmu;
  NTuple::Item<double> m_cthe_kmu;
  NTuple::Item<double> m_cphi_kmu;
  NTuple::Item<double> m_costheta_kmu;
  NTuple::Item<double> m_cang_kpi;
  NTuple::Item<double> m_cthe_kpi;
  NTuple::Item<double> m_cphi_kpi;
  NTuple::Item<double> m_costheta_kpi;
  NTuple::Item<double> m_cang_kk;
  NTuple::Item<double> m_cthe_kk;
  NTuple::Item<double> m_cphi_kk;
  NTuple::Item<double> m_costheta_kk;

  NTuple::Array<double> m_m_egam_pos;
  NTuple::Array<double> m_m_mugam_pos;
  NTuple::Array<double> m_m_pigam_pos;
  NTuple::Array<double> m_m_kgam_pos;
  NTuple::Array<double> m_m_egam_neg;
  NTuple::Array<double> m_m_mugam_neg;
  NTuple::Array<double> m_m_pigam_neg;
  NTuple::Array<double> m_m_kgam_neg;

  ITHistSvc *m_histSvc;

  TH1F *m_cutflow;

  LocalPhotonSelector photonSelector;
};

inline void tautauAlg::setPxPyPzE(HepLorentzVector &p4, double px, double py,
                                  double pz, double e) const {
  p4.setPx(px);
  p4.setPy(py);
  p4.setPz(pz);
  p4.setE(e);
}
#endif
