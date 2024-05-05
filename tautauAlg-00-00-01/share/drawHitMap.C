

void drawHitMap(){

TFile* file = TFile::Open("tau_decay_pi_pi0_nu.root", "read");


TH1F* h_nRecEmcHits = new TH1F("nRecEmcHits_hist", "", 500, 0, 500);

//barrel
TH2F* id_map_barrel = new TH2F("id_barrel", "", 43, 0, 43, 119, 0, 119);
//endcap east
TH2F* id_map_east = new TH2F("id_east", "", 5, 0, 5, 95, 0, 95);
//endcap west
TH2F* id_map_west = new TH2F("id_west", "", 5, 0, 5, 95, 0, 95);

TH2F* pos_map = new TH2F("pos_map", "", 40, 0, M_PI, 40, -M_PI, M_PI);


TTreeReader*  treeReader = new TTreeReader("rec", file); 
TTreeReaderValue<int>* nGam = new TTreeReaderValue<int>(*treeReader, "numGam");
TTreeReaderValue<int>* nPi0 = new TTreeReaderValue<int>(*treeReader, "numpi0");
TTreeReaderValue<int>* nRecEmcHits = new TTreeReaderValue<int>(*treeReader, "nRecEmcHits");
TTreeReaderArray<int>* id_theta = new TTreeReaderArray<int>(*treeReader, "emcHit_id_theta");
TTreeReaderArray<int>* id_phi = new TTreeReaderArray<int>(*treeReader, "emcHit_id_phi");
TTreeReaderArray<int>* bc = new TTreeReaderArray<int>(*treeReader, "emcHit_bc");
TTreeReaderArray<int>* tdc = new TTreeReaderArray<int>(*treeReader, "emcHit_tdc");
TTreeReaderArray<double>* energy = new TTreeReaderArray<double>(*treeReader, "emcHit_energy");
TTreeReaderArray<double>* pos_theta = new TTreeReaderArray<double>(*treeReader, "emcHit_pos_theta");
TTreeReaderArray<double>* pos_phi = new TTreeReaderArray<double>(*treeReader, "emcHit_pos_phi");

int ievent =0;
while (treeReader->Next()) {
   h_nRecEmcHits->Fill(*(*nRecEmcHits)); 
   if(ievent==5){
     std::cout<<"nGam = " << *(*nGam) <<", npi0 = " << *(*nPi0) << std::endl; 
     for(int ih=0; ih<*(*nRecEmcHits); ++ih){
       if((*bc)[ih]==1)
       id_map_barrel->Fill((*id_theta)[ih], (*id_phi)[ih], (*energy)[ih]);
       if((*bc)[ih]==0)
       id_map_east->Fill((*id_theta)[ih], (*id_phi)[ih], (*energy)[ih]);
       if((*bc)[ih]==2)
       id_map_west->Fill((*id_theta)[ih], (*id_phi)[ih], (*energy)[ih]);

       pos_map->Fill((*pos_theta)[ih], (*pos_phi)[ih], (*energy)[ih]);
     }
   }
   ievent++; 
}

std::cout<<"ievent = " << ievent << std::endl;


TCanvas* c1 = new TCanvas("c1", "", 600, 500);
h_nRecEmcHits->Draw();

TCanvas* c2 = new TCanvas("c2", "", 600, 500);
id_map_west->Draw("colz");

TCanvas* c3 = new TCanvas("c3", "", 600, 500);
id_map_barrel->Draw("colz");

TCanvas* c4 = new TCanvas("c4", "", 600, 500);
id_map_east->Draw("colz");

TCanvas* c5 = new TCanvas("c5", "", 600, 500);
pos_map->Draw("colz");

}
