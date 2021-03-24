#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TPave.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLatex.h"

const double L = 220; // cm


 

///////////////////////////

const double Vinf = 15.90; // velocity [cm/s]
//const double Vinf = 17.00;

//double VINF[20]= {Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf};

// TOF1
// Alexander : double VINF[20]= {16.014,15.964,15.912,15.906,15.870,15.872,15.889,15.724,15.924,15.878,15.934,15.929,15.889,15.839,15.821,15.958,15.939,15.919,15.879,16.014};
//double VINF[20]= { 15.8956, 15.902, 15.8996, 15.899, 15.8961, 15.905, 15.9036, 15.9009, 15.8992, 15.9002, 15.9019, 15.8968, 15.9094, 15.901, 15.8986, 15.9021, 15.8986, 15.899, 15.8953, 15.8938  }; 
//double VINF[20]= { 15.891, 15.904, 15.899, 15.897, 15.892, 15.910, 15.907, 15.903, 15.898, 15.901, 15.904, 15.894, 15.919, 15.901, 15.898, 15.904, 15.897, 15.898, 15.891, 15.887  }; 
//double VINF[20]= {15.904,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf};
double VINF[20]= {16.014,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf,Vinf};


// TOF2
//double VINF[20]= {16.053,15.896,15.908,15.902,15.908,15.829,15.815,15.798,15.880,15.915,15.907,15.911,15.860,15.926,15.896,15.912,15.858,15.968,15.940,15.989};

// TOF3
//double VINF[20]= {16.010,15.917,15.944,15.890,15.865,15.897,15.889,15.889,15.886,15.864,15.867,15.827,15.871,15.828,15.904,15.960,15.976,15.966,16.023,15.871};

// TOF4
//double VINF[20]= {16.079,16.034,15.859,15.749,16.070,15.985,15.959,15.949,15.845,15.817,15.888,15.863,15.932,15.809,15.922,15.965,15.820,15.999,16.016,16.079};



///////////////////////////

//double DX[20]= {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

const double Dx = 3.5/2;
//double DX[20]= {-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx};

// TOF1
// Alexander --> double DX[20]= {-1.091,2.903,-1.356,2.045,-2.117,1.302,-2.404,0.624,-1.570,1.623,-2.274,1.932,-1.187,2.086,-1.671,1.473,-1.302,1.068,-1.904,1.678};
//double DX[20]= {-1.63831, 1.89605, -2.15735, 2.30059, -1.90448, 1.93373, -1.29697, 2.35799, -1.52673, 1.40584, -1.89098, 1.82633, -2.78917, 1.09952, -2.21297, 1.35722, -1.60793, 1.89879, -0.781469, 2.23056};
//double DX[20]= { -1.935, 1.098, -2.732, 2.759, -1.683, 2.567, -0.257, 3.852, -1.367, 1.177, -1.755, 1.699, -4.374, 0.039, -2.827, 1.113, -1.916, 2.516, 0.100, 2.635  };
//double DX[20]= {-1.862,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx};
double DX[20]= {-1.091,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx,-Dx,Dx};

// TOF2
//double DX[20]= {-1.296, 2.740,-1.581,1.946,-2.815,1.418,-1.733,1.516,-1.371,1.971,-2.028,1.632,-0.674,1.112,-1.985,1.980,-1.081,1.306,-1.760,0.637};

// TOF3
//double DX[20]= {-1.369,2.804,-1.739,1.807,-2.814,1.383,-1.969,1.880,-1.250,1.546,-2.458,1.860,-1.209,2.092,-1.147,1.949,-0.960,0.699,-1.913,0.764};

// TOF4
//double DX[20]= {-1.310,2.545,-1.338,1.852,-3.178,1.850,-1.758,1.258,-1.018,1.909,-2.514,1.859,-1.627,2.468,-2.174,1.301,-1.251,1.107,-2.163,1.389};




////////////////////////////////////////////////////////////////

double GetV( double x = L, int ibar = 0 );
double GetL1(double t1,double t2, int ibar = 0 );
double GetL2(double t2, int ibar = 0 );
void DrawBars();



void DrawTrackOnly(TString sname="")
{
    gStyle->SetOptStat(1111);
    
    gStyle->SetOptFit(11111);
    /////////////////// in ///////////////////////////////////

    //sname="../datafiles/Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF1_Data_7_21_2020_Binary.root";
    //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF1_Data_7_21_2020_Binary_1track.root";
    //sname="./TOF1_default_v15p90_and_pos_selectedtracks.root";
    //sname= "./TOF1_default_v15p90_and_pos_selectedtracks_cal2.root";
    sname = "./TOF1_cal_bar0.root";


    //sname="../datafiles/Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF2_Data_7_31_2020_Binary.root";
    //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF2_Data_7_31_2020_Binary_1track.root";
    
    //sname="../datafiles/Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF3_Data_8_6_2020_Binary.root";
    //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF3_Data_8_6_2020_Binary_1track.root";
    
    //sname="../datafiles/Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF4_Data_8_12_2020_Binary.root";
    //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF4_Data_8_12_2020_Binary_1track.root";
    
    //this has 36 channels!
    //sname="../datafiles/Run_bar220cm_110cm_36ch_2sides_V111_bot_TOF5_Data_9_30_2020_Binary.root";

    //sname="../datafiles/Run_bar220cm_110cm_40ch_2sides_V111_top_TOF6_Data_10_2_2020_Binary.root";

    
    ///////////////// Open input file /////////////////////////////////////
    
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(sname.Data());
    if (!f) f = new TFile(sname.Data());
    if( f->IsOpen() == 0 ) { printf("\nERROR: No file %s\n\n", sname.Data() ); exit(0); }
    printf("\nOpen file: %s\n", sname.Data() );
    
    
    /////////////// input tree ///////////////////////////////////////////
    
    f->cd();
    TTree *tree = 0;
    tree = (TTree*)f->Get("tree");

   

    const int NCH=40;
    
    //Declaration of leaves types
    Int_t           iev;
    Int_t           nch;
    Int_t           ich[NCH];
    Double_t        A[NCH];
    Double_t        AtrMin;
    Double_t        tbin[NCH];
    Double_t        t6[NCH];
    Double_t        t8[NCH];
    Double_t        t10[NCH];
    Double_t        t12[NCH];
    Double_t        t20[NCH];
    Double_t        t30[NCH];
    Double_t        tLE10[NCH];
    Double_t        tLE20[NCH];
    Double_t        tLE30[NCH];
    Double_t        tToT10[NCH];
    Double_t        tToT20[NCH];
    Double_t        tToT30[NCH];
    Double_t        tr6;
    Double_t        tr8;
    Double_t        tr10;
    Double_t        tr12;
    Double_t        tr20;
    Double_t        tr30;
    Double_t        dtr6;
    Double_t        dtr8;
    Double_t        dtr10;
    Double_t        dtr12;
    Double_t        dtr20;
    Double_t        dtr30;
    Double_t        BL[NCH];
    Double_t        RMS[NCH];
    Double_t        rt[NCH];
    Double_t        Max[NCH];
    Double_t        Min[NCH];
    
    // Set branch addresses.
    tree->SetBranchAddress("iev",&iev);
    tree->SetBranchAddress("nch",&nch);
    tree->SetBranchAddress("ich",ich);
    tree->SetBranchAddress("A",A);
    tree->SetBranchAddress("AtrMin",&AtrMin);
    tree->SetBranchAddress("tbin",tbin);
    tree->SetBranchAddress("t6",t6);
    tree->SetBranchAddress("t8",t8);
    tree->SetBranchAddress("t10",t10);
    tree->SetBranchAddress("t12",t12);
    tree->SetBranchAddress("t20",t20);
    tree->SetBranchAddress("t30",t30);
    tree->SetBranchAddress("tLE10",tLE10);
    tree->SetBranchAddress("tLE20",tLE20);
    tree->SetBranchAddress("tLE30",tLE30);
    tree->SetBranchAddress("tToT10",tToT10);
    tree->SetBranchAddress("tToT20",tToT20);
    tree->SetBranchAddress("tToT30",tToT30);
    tree->SetBranchAddress("tr6",&tr6);
    tree->SetBranchAddress("tr8",&tr8);
    tree->SetBranchAddress("tr10",&tr10);
    tree->SetBranchAddress("tr12",&tr12);
    tree->SetBranchAddress("tr20",&tr20);
    tree->SetBranchAddress("tr30",&tr30);
    tree->SetBranchAddress("dtr6",&dtr6);
    tree->SetBranchAddress("dtr8",&dtr8);
    tree->SetBranchAddress("dtr10",&dtr10);
    tree->SetBranchAddress("dtr12",&dtr12);
    tree->SetBranchAddress("dtr20",&dtr20);
    tree->SetBranchAddress("dtr30",&dtr30);
    tree->SetBranchAddress("BL",BL);
    tree->SetBranchAddress("RMS",RMS);
    tree->SetBranchAddress("rt",rt);
    tree->SetBranchAddress("Max",Max);
    tree->SetBranchAddress("Min",Min);
    
    // tree->SetBranchStatus("*",0);  // disable all branches
    // TTreePlayer->SetBranchStatus("branchname",1);  // activate branchname
    
    TFile newfile("small.root", "recreate");
    auto newtree = tree->CloneTree(0);
  
    
    //////////////////// Make histograms /////////////////////////////////////////////
    
    char ctmp[100];
    
    //////////////////////////////////
    
    double AMAX=1.2;
    int ANB=100;
    
    TH1 *hA[100];
    for( int i=0; i<NCH; i++ ) {
        sprintf(ctmp,"hA%i",i);
        hA[i] = (TH1*)gDirectory->Get(ctmp);
        if(hA[i]!=0) delete hA[i];
        hA[i] = new TH1F(ctmp,"",ANB,0.,AMAX);
        hA[i]->SetFillColor(kYellow);
        sprintf(ctmp,"Amplitude-%i [V]",i);
        hA[i]->SetXTitle(ctmp);
    }
    
    
    TH1 *hDTL[100];
    for( int i=0; i<NCH/2; i++ ) {
        sprintf(ctmp,"hDTL%i",i);
        hDTL[i] = (TH1*)gDirectory->Get(ctmp);
        if(hDTL[i]!=0) delete hDTL[i];
        hDTL[i] = new TH1F(ctmp,"",100, -10,10 );
        hDTL[i]->SetFillColor(kYellow);
        sprintf(ctmp,"dt [ns]");
        hDTL[i]->SetXTitle(ctmp);
    }

    
    TH1 *hDT[100];
    for( int i=0; i<NCH/2-1; i++ ) {
        sprintf(ctmp,"hDTL_%i_%i",i,i+1);
        hDT[i] = (TH1*)gDirectory->Get(ctmp);
        if(hDT[i]!=0) delete hDT[i];
        hDT[i] = new TH1F(ctmp,"",100, -5,5 );
        hDT[i]->SetFillColor(kYellow);
        sprintf(ctmp,"dt [ns]");
        hDT[i]->SetXTitle(ctmp);
    }


    TH2 *hX_vs_Xtr[100];
    for( int i=0; i<NCH/2; i++ ) {
        sprintf(ctmp,"hX_vs_Xtr_%i",i);
        //std::cout << ctmp << std::endl;
        hX_vs_Xtr[i] = (TH2*)gDirectory->Get(ctmp);
        if(hX_vs_Xtr[i]!=0) delete hX_vs_Xtr[i];
        hX_vs_Xtr[i] = new TH2F(ctmp,"",230,-115,115, 20,-10, 10 );
        hX_vs_Xtr[i]->SetMarkerStyle(21);
        hX_vs_Xtr[i]->SetMarkerSize(0.6);
        sprintf(ctmp,"x_{hit} [cm]");
        hX_vs_Xtr[i]->SetXTitle(ctmp);
        sprintf(ctmp,"X_{hit} - X_{track} [cm]");
        hX_vs_Xtr[i]->SetYTitle(ctmp);
    }

     

    TH1 *hXres[100];
    for( int i=0; i<NCH/2; i++ ) {
        sprintf(ctmp,"hXres_%i",i);
        hXres[i] = (TH1*)gDirectory->Get(ctmp);
        if(hXres[i]!=0) delete hXres[i];
        hXres[i] = new TH1F(ctmp,"Xtrack - X hit residuals",200, -10,10 );
        hXres[i]->SetLineWidth(3);
        sprintf(ctmp,"Xtrack - X hit [cm]");
        hXres[i]->SetXTitle(ctmp);
    }
    

    
    sprintf(ctmp,"hXY");
    TH2 *hXY = (TH2*)gDirectory->Get(ctmp);
    if(hXY!=0) delete hXY;
    hXY = new TH2F(ctmp,"Y vs X upside-down",2*7*5, -7,7, NCH/2,0,NCH/2 );
    hXY->SetFillColor(kRed);
    hXY->SetXTitle("X [ns]");
    hXY->SetYTitle("bar number");
    
    
    sprintf(ctmp,"hXYl");
    TH2 *hXYl = (TH2*)gDirectory->Get(ctmp);
    if(hXYl!=0) delete hXYl;
    hXYl = new TH2F(ctmp,"Y vs X upside-down", 120,-115,115, NCH/2+1,-1.5,NCH/2-0.5 );
    hXYl->GetYaxis()->SetNdivisions(21);
    hXYl->GetXaxis()->SetLabelSize(0.03);
    hXYl->GetYaxis()->SetLabelSize(0.025);
    hXYl->GetXaxis()->SetTickLength(0.02);
    //hXYl->GetYaxis()->SetBinLabel(1,"x");
    hXYl->GetYaxis()->SetAxisColor(kWhite);
    hXYl->SetFillColor(kRed);
    hXYl->SetXTitle("x [cm]");
    hXYl->SetYTitle("bar number");
    //
    sprintf(ctmp,"hXYl_all");
    TH2 *hXYl_all = (TH2*)gDirectory->Get(ctmp);
    if(hXYl_all!=0) delete hXYl_all;
    hXYl_all = new TH2F(ctmp,"Y vs X upside-down", 120,-115,115, NCH/2+1,-1.5,NCH/2-0.5 );
    hXYl_all->GetYaxis()->SetNdivisions(21);
    hXYl_all->GetXaxis()->SetLabelSize(0.03);
    hXYl_all->GetYaxis()->SetLabelSize(0.025);
    hXYl_all->GetXaxis()->SetTickLength(0.02);
    //hXYl_all->GetYaxis()->SetBinLabel(1,"x");
    hXYl_all->GetYaxis()->SetAxisColor(kWhite);
    hXYl_all->SetFillColor(kRed);
    hXYl_all->SetXTitle("x [cm]");
    hXYl_all->SetYTitle("bar number");
    //
    sprintf(ctmp,"hXYl_inv");
    TH2 *hXYl_inv = (TH2*)gDirectory->Get(ctmp);
    if(hXYl_inv!=0) delete hXYl_inv;
    hXYl_inv = new TH2F(ctmp,"Y vs X", 120,-115,115, NCH/2+1,-1.5,NCH/2-0.5 );
    hXYl_inv->GetYaxis()->SetNdivisions(21);
    hXYl_inv->GetXaxis()->SetLabelSize(0.03);
    hXYl_inv->GetYaxis()->SetLabelSize(0.);
    hXYl_inv->GetYaxis()->SetLabelColor(kWhite);
    hXYl_inv->GetXaxis()->SetTickLength(0.02);
    //hXYl_inv->GetYaxis()->SetBinLabel(1,"x");
    hXYl_inv->GetYaxis()->SetAxisColor(kWhite);
    hXYl_inv->SetFillColor(kRed);
    hXYl_inv->SetXTitle("x [cm]");
    hXYl_inv->SetYTitle("bar number");
    
    
    sprintf(ctmp,"hChi2");
    TH1 *hChi2 = (TH1*)gDirectory->Get(ctmp);
    if(hChi2!=0) delete hChi2;
    hChi2 = new TH1F(ctmp,"", 500,0,10 );
    hChi2->GetYaxis()->SetNdivisions(20);
    hChi2->SetFillColor(kGreen);
    hChi2->SetXTitle("#chi^{2}");
    
    
    TF1* fun_tr = new TF1("fun_tr","pol1",-110,110);
    fun_tr->SetLineColor(kBlue);
    
    
    /*sprintf(ctmp,"hX_Xtr");
    TH2 *hX_Xtr = (TH2*)gDirectory->Get(ctmp);
    if(hX_Xtr!=0) delete hX_Xtr;
    hX_Xtr = new TH2F(ctmp,"X vs Xtr", 120,-115,115, 20,-10, 10);
    //hX_Xtr->GetYaxis()->SetNdivisions(21);
    hX_Xtr->GetXaxis()->SetLabelSize(0.03);
    hX_Xtr->GetYaxis()->SetLabelSize(0.03);
    //hX_Xtr->GetYaxis()->SetLabelSize(0.);
    //hX_Xtr->GetYaxis()->SetLabelColor(kWhite);
    hX_Xtr->GetXaxis()->SetTickLength(0.02);
    hX_Xtr->GetYaxis()->SetTickLength(0.02);
    //hXYl_inv->GetYaxis()->SetBinLabel(1,"x");
    //hX_Xtr->GetYaxis()->SetAxisColor(kWhite);
    //hX_Xtr->SetFillColor(kRed);
    hX_Xtr->SetXTitle("x_{hit} [cm]");
    hX_Xtr->SetYTitle("X_{hit} - X_{track}");
    */

    sprintf(ctmp,"hX_Xtr_sel");
    TH2 *hX_Xtr_sel = (TH2*)gDirectory->Get(ctmp);
    if(hX_Xtr_sel!=0) delete hX_Xtr_sel;
    hX_Xtr_sel = new TH2F(ctmp,"X vs Xtr (chi2<4)", 120,-115,115, 20,-10, 10);
    hX_Xtr_sel->GetXaxis()->SetLabelSize(0.03);
    hX_Xtr_sel->GetYaxis()->SetLabelSize(0.03);
    hX_Xtr_sel->GetXaxis()->SetTickLength(0.02);
    hX_Xtr_sel->GetYaxis()->SetTickLength(0.02);
    hX_Xtr_sel->SetXTitle("x_{hit} [cm]");
    hX_Xtr_sel->SetYTitle("X_{hit} - X_{track}");
    



    TH1F *hX   = new TH1F("hX", " X hit   ", 120, -115, 115);
    TH1F *hXtr = new TH1F("hXtr", " Xtrack   ", 120, -115, 115);
    //TH1F *hXres = new TH1F("hXres", " Xtrack - X hit residuals   ", 200, -10, 10);

    TH1F *hXres_sel = new TH1F("hXres_sel", " Xtrack - X hit residuals  (chi2<4) ", 200, -10, 10);  
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TCanvas *c1 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c1");
    if(c1!=0) delete c1;
    //c1 = new TCanvas("c1","c1",0,0,800,520);
    c1 = new TCanvas("c1","c1",0,0, 600,600);
    c1->cd();
    
    gPad->Clear();
    gPad->SetFillColor(0);
    
    //gPad->SetFrameFillColor(0);
    //gPad->SetFillColor(0);
    
    
    //gPad->SetBorderMode(0);
    //gPad->SetBorderSize(0);
    gPad->SetFrameBorderMode(0);
    gPad->SetFrameBorderSize(0);
    //gPad->SetFrameLineWidth(2);
    gPad->SetFrameLineWidth(0);
    
    gPad->SetLeftMargin(0.08);
    gPad->SetRightMargin(0.02);
    gPad->SetTopMargin(0.07);
    gPad->SetBottomMargin(0.08);
    
    //gPad->SetGridy();
    
    
    
    /////////////////////////////////////////////////
    ////////////////////////////////  Loop over events
    
    int NSEL = 0;
    
    Long64_t nentries = tree->GetEntries();
    //nentries=2000;
    printf("Number of events in file %i\n",(int)nentries);

    //double a_bar_v1sum[20]; // array storing the sum of the estimated v1 for each bar
    //double a_bar_v2sum[20]; // array storing the sum of the estimated v2 for each bar --> to do the average at the end
    //int  a_bar_evts[20];  // to count how many hits we have in that bar (the selection can kill some hits so nhits != nentries)

    double a_par0[20];
    double a_par1[20];
    double a_bar_number[20];  

    double a_gaus_mu[20];
    double a_gaus_sigma[20];


    /*
    for (int ibar =0; ibar < 20; ibar++){
        a_bar_v1sum[ibar] = 0;   
        a_bar_v2sum[ibar] = 0;  
        a_bar_evts[ibar]  = 0;   
        a_bar_number[ibar] = ibar;
    }*/


    for ( Long64_t IE=0; IE<nentries; IE++ ) {
        
        
        if( tree->GetEntry(IE) == 0 ) break;
        
        printf("---------- %i ----------\n",(int)IE);
        
        int NGOOD = 0;
        double X[NCH];
        
        for( int i=0; i<NCH; i++ ) {
            
            if(A[i]>0.01) hA[i]->Fill( A[i] );
            
            if(A[i]>0.01)
            if(i<10) hDTL[i]->Fill( (t8[i]-t8[i+10])/2 );
            
            
            double cutA = 0.5; // the amplitude cut in Volts
           
            if( i%2==0
               && A[i]>cutA && A[i+1]>cutA
               ) {
                
                int ibar = ich[i]/2;
                
                hXY->Fill( (t8[i]-t8[i+1])/2, ibar );
                
                double dxb = DX[ibar];
                
                double l1_b1 = GetL1(t8[i],t8[i+1],ibar) - dxb;
                
                //
                //double v1 = GetV(l1_b1,ibar);
                //double v2 = GetV(L-l1_b1,ibar);
                //std::cout << " bar  " << ibar << "  v1 " << v1 << "  v2 " << v2 <<  std::endl;    

                //a_bar_v1sum[ibar]= a_bar_v1sum[ibar] + v1;
                //a_bar_v2sum[ibar]= a_bar_v2sum[ibar] + v2;
                
                //a_bar_evts[ibar] = a_bar_evts[ibar]+1;
                
                //

                l1_b1 -= (L+2*2.5)/2;
                //X[i] = l1_b1;
                X[NGOOD] = l1_b1;
                
                hXYl    ->Fill( l1_b1, ibar );
                hXYl_all->Fill( l1_b1, ibar );
                hXYl_inv->Fill( l1_b1, 19-ibar );
                //std::cout << "l1_b1  "  << l1_b1 << std::endl;
                //std::cout << " NGOOD " << NGOOD<< "X[NGOOD]  "  << X[NGOOD] << std::endl;
                NGOOD++;
            }
            
        } // end of loop over channels
        
        
        double xtr[20];
        for( int i=0; i<20; i++ ) xtr[i] = 999;
        
        if(NGOOD==20) {
    
            newtree->Fill();

            //hXY->Draw("box");
            //hXYl->Draw("box");
            hXYl_inv->Draw("box");
            
            DrawBars();
           
            hXYl_inv->Fit("fun_tr","Q");
            
            double chi2 = fun_tr->GetChisquare() / fun_tr->GetNDF();
            //printf("%.2f\n",chi2);
            
            sprintf(ctmp,"sel %i, ev %i",NSEL,(int)IE);
            //sprintf(ctmp,"sel %i",NSEL);
            hXYl->SetTitle(ctmp);
            
            hChi2->Fill( chi2 );
            
            // line assumed
            double a = fun_tr->GetParameter(1);
            double b = fun_tr->GetParameter(0);
            for(int itr=0; itr<20; itr++) {
                xtr[itr] = (itr-b)/a;
                //printf("%2i %.2f\n",itr,xtr[itr]);
                
            }
            
            
            
            if( chi2 < 100 ) {
                
                c1->Update(); gSystem->ProcessEvents(); gSystem->Exec("sleep 1");
              
                // make pdf file
                //sprintf(ctmp,"ev_%i.pdf",(int)IE); c1->Print(ctmp);
                //sprintf(ctmp,"ev_%i_%i.pdf",NSEL,(int)IE); c1->Print(ctmp);
                
                printf("%i  %.2f\n",(int)IE,chi2);
                
            
                for(int itr=0; itr<20; itr++) {
                
 
                hX_vs_Xtr[itr]->Fill( X[19-itr], X[19-itr]-xtr[itr]);
                

                hX->Fill( X[itr]);
                hXtr->Fill(xtr[itr]);
                hXres[itr]->Fill(X[19-itr]-xtr[itr]);

                if( chi2 < 4 ) {    
                    hX_Xtr_sel->Fill( X[19-itr], X[19-itr]-xtr[itr]);
                    hXres_sel->Fill(X[19-itr]-xtr[itr]);

                }


                }   


                NSEL++;
            }
            
        } // end of if NGOOD==20
        
  
        if( IE != nentries-1 ) {
            hXY->Reset();
            hXYl->Reset();
            hXYl_inv->Reset();
        }
        
    } // end of loop over IE
    



    //double a_bar_v1mean[20] ;
    //double a_bar_v2mean[20] ;

    
    // printout the mean value for the v1 and v2 for each bar
    //for (int ibar=0; ibar < 20; ibar++){
       
        //a_bar_v1mean[ibar] = a_bar_v1sum[ibar]/a_bar_evts[ibar];
        //a_bar_v2mean[ibar] = a_bar_v2sum[ibar]/a_bar_evts[ibar];


        //std::cout << " bar " << ibar << ":  v1_mean " << a_bar_v1mean[ibar] << ";  v2_mean " << a_bar_v2mean[ibar] << std::endl;
    }//   

    
    //////////////////////////////////
  

    TCanvas *c2 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c2");
    if(c2!=0) delete c2;
    //c1 = new TCanvas("c1","c1",0,0,800,520);
    c2 = new TCanvas("c2","c2",0,0, 600,600);
    c2->cd();


    hXYl_all->Draw("box");
    
    TCanvas *c3 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c3");
    if(c3!=0) delete c3;
    //c1 = new TCanvas("c1","c1",0,0,800,520);
    c3 = new TCanvas("c3","c3",0,0, 600,600);
    c3->cd();
    
    hChi2->Draw("");
    

    TCanvas *c4 = (TCanvas*)gROOT->GetListOfCanvases()->FindObject("c4");
    if(c4!=0) delete c4;
    //c1 = new TCanvas("c1","c1",0,0,800,520);
    c4 = new TCanvas("c4","c4",0,0, 600,600);
    c4->cd();

    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);
    
    for (int ibar=0; ibar<20; ibar++){
    //hX_Xtr->SetMarkerStyle(21);
    //hX_Xtr->SetMarkerSize(0.6);
    
        std::cout << " bar " << ibar << std::endl; 
 
        hX_vs_Xtr[ibar]->Draw();

        TF1* fun_pol1 = new TF1("fun_pol1","pol1",-110,110);
        fun_pol1->SetLineColor(kRed);

        //hX_vs_Xtr[ibar]->Fit("pol1"); 
        
        hX_vs_Xtr[ibar]->Fit("fun_pol1");     

        //retrieve the values and store them in array to have an easy printout at the end to use in the macro
        a_par0[ibar] = fun_pol1->GetParameter(0); 
        a_par1[ibar] = fun_pol1->GetParameter(1);

        std::cout<< " bar : " << ibar << "  par0:  "<< fun_pol1->GetParameter(0) << "  par 1:  " << fun_pol1->GetParameter(1) << std::endl; 

    }

    // execute the printout

    std::cout << " par 0 [X pos]: [ " << std::endl;
    std::cout << " { " ; 
    std::cout << std::setprecision(3) << std::fixed;    

    for (int ibar=0; ibar<20; ibar++) { 
        std::cout << DX[ibar] - a_par0[ibar] << ", ";
    }
    std::cout << " }; " << std::endl;

    std::cout << " par 1 [Vinf]: "<< std::endl; 
    std::cout << " { " ;
    std::cout << std::setprecision(3) << std::fixed;

    for (int ibar=0; ibar<20; ibar++) { 
        std::cout << VINF[ibar] - a_par1[ibar] << ", ";
    }
    std::cout << " }; " << std::endl;


    TCanvas *c5 = new TCanvas("c5", "c5 ", 600, 600);
    c5->cd();
    hX->SetLineColor(4);
    hX->Draw();
    hXtr->SetLineColor(kGreen-6);
    hXtr->Draw("same");
    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->AddEntry("hX", "X hit", "l");
    leg->AddEntry("hXtr", "X track", "l");
    leg->Draw("same");

    TCanvas *c6 = new TCanvas("c6", "c6 ", 600, 600);
    c6->cd();
    gStyle->SetOptFit(11111);
    
    for (int ibar =0; ibar<20; ibar ++){
    hXres[ibar]->SetLineColor(1);
    hXres[ibar]->SetLineWidth(3);

    hXres[ibar]->Draw();

    TF1* fun_gaus = new TF1("fun_gaus","gaus",-10,10);
    fun_gaus->SetLineColor(kRed);
    hXres[ibar]->Fit("fun_gaus");     

    a_gaus_mu[ibar]    = fun_gaus->GetParameter(1); 
    a_gaus_sigma[ibar] = fun_gaus->GetParameter(2);
    }


    TGraph *g_gaus_mean = new TGraph(20, a_bar_number, a_gaus_mu);
    TGraph *g_gaus_sigma = new TGraph(20, a_bar_number, a_gaus_sigma);
    
    g_gaus_mean->SetMarkerStyle(21);  
    g_gaus_mean->SetMarkerColor(kBlue);  
    g_gaus_mean->SetMarkerSize(0.6);

    g_gaus_sigma->SetMarkerStyle(21);  
    g_gaus_sigma->SetMarkerColor(kOrange);  
    g_gaus_sigma->SetMarkerSize(0.6);     

    g_gaus_mean->Draw("AP");
    g_gaus_sigma->Draw("P same");

    g_gaus_mean->SetName("g_gaus_mean");
    g_gaus_sigma->SetName("g_gaus_sigma");

    g_gaus_mean->Write();
    g_gaus_sigma->Write();
    

    TCanvas *c7 = new TCanvas("c7", "c7 ", 600, 600);
    c7->cd();
    hXres_sel->SetLineColor(1);
    hXres_sel->SetLineWidth(3);

    hXres_sel->Draw();

    TCanvas *c8 = new TCanvas("c8", "c8 ", 600, 600);
    c8->cd();
    hX_Xtr_sel->Draw();



    //////////////////////////////////
    
    printf("NSEL = %i\n",NSEL);
    
    //c1->Print("c1.pdf"); printf(".! open c1.pdf\n");
    newfile.Write();


    return;
}






void DrawBars()
{
    double x1,x2,y1,y2;
    
    //for(int i=0;i<20;i++) {
    for(int i=19;i>-1;i--) {
        x1 = -L/2;
        x1 += DX[i]/2;
        x2 = x1+L;
        y1 = i-0.5;
        y2 = y1+1;
        
        TPave *pave = new TPave( x1,y1, x2,y2, 1,"br");
        pave->SetFillStyle(0);
        pave->SetLineColor(kBlack);
        pave->Draw();
    }
    
    char ctmp[100];
    //for(int i=0;i<20;i++) {
    for(int i=19;i>-1;i--) {
        sprintf(ctmp,"%i",19-i);
        TLatex* tex = new TLatex( -120, i-0.2, ctmp);
        tex->SetTextFont(42);
        tex->SetTextSize(0.025);
        tex->Draw();
    }
    
    return;
}


double GetV( double x, int ibar)
{
    //return 16;
    
    const double v0 = 14.54;
    //double vinf = ::Vinf;
    double vinf = ::VINF[ibar];
    
    //printf("%i %f\n",ibar,vinf);
    
    double v = vinf + ( v0 - vinf ) * exp( - x / 55.51 );
    
    return v;
}


double GetL1(double t1,double t2, int ibar )
{
    double v1 = GetV(L,ibar);
    double v2 = v1;
    double l1 = ( (t1-t2)*v2 + L ) / 2;
    //return l1;
    
    v1 = GetV(l1,ibar);
    v2 = GetV(L-l1,ibar);
    l1 = v1/(v1+v2) * ( (t1-t2)*v2 + L );

    //std::cout << " bar  " << ibar << "  v1 " << v1 << "  v2 " << v2 <<  std::endl;

    return l1;
}


/*
double GetL1(double t1, int ibar = 0 )
{
    double v1 = GetV(L,ibar);
    t1 -= 82;
    double l1 = t1*v1;
    
    printf("l1=%f\n",l1);
    return l1;
    
    //v1 = GetV(l1,ibar);
    //v2 = GetV(L-l1,ibar);
    //l1 = v1/(v1+v2) * ( (t1-t2)*v2 + L );
    
    return l1;
}
*/


double GetL2(double t2, int ibar )
{
    double v2 = GetV(L,ibar);
    t2 -= 70;
    double l2 = L-t2*v2;
    
    printf("l2=%f\n",l2);
    return l2;
    
    //v1 = GetV(l1,ibar);
    //v2 = GetV(L-l1,ibar);
    //l1 = v1/(v1+v2) * ( (t1-t2)*v2 + L );
    
    return l2;
}



int IsSingleTrack(double X[])
{
    return 1;
}


void ChangeRange(TF1* fun)
{
    fun->SetRange(-50,-40);
    return;
}

