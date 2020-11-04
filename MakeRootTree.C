#include "TString.h"
#include "TFile.h"
//#include "TString.h"
//#include "TCanvas.h"
//#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
//#include "TMath.h"
//#include "TGraphErrors.h"
//#include "TStyle.h"
//#include "TRint.h"
#include "TROOT.h"
#include "TSystem.h"
//#include "TLine.h"
//#include "Math/Integrator.h"
//#include "TSpline.h"
//#include "TProfile.h"
//#include "TArrow.h"
//#include "TLatex.h"
//#include "TDatime.h"
//#include "TFile.h"
//#include "TLegend.h"


#include "read_data_WC.h"
//#include "/Users/alex/work/DRS4/C_scripts/read_data_WC.h"

int NWRITTEN = 0;
int NTOT = 0;


const int DATA_SAMPLES = 1024;

double baseline[N_CHN_MAX];
double Maximum[N_CHN_MAX];
double t6[N_CHN_MAX],t8[N_CHN_MAX],t10[N_CHN_MAX],t12[N_CHN_MAX], t20[N_CHN_MAX], t30[N_CHN_MAX];
double tbin[N_CHN_MAX];

double tr6,tr8,tr10,tr12,tr20,tr30;
double dtr6,dtr8,dtr10,dtr12,dtr20,dtr30;
double Amp[N_CHN_MAX];
double rt[N_CHN_MAX];

const double c = 299792458; // m/s


//////////////////////////////////
inline void GetBaseLine(int ich, double tl,double tr, double tbin, double&, double& );
inline void GetMaximum(int ich, int& im, double& max);
inline void GetMinimum(int ich, int& im, double& min);
void GetTimeCFD(int ich,double thr,double bl, double tbin, double&t,double&A);
void GetTimeLE(int ich,double thr, double, double, double&t, double& );
inline double GetRiseTime(int ich,double tbin);
void ResetSamplingFrequency(int ich, int ISFS);


void MakeRootTree(TString sname="")
{
    //gStyle->SetOptStat(11111111);
    //gStyle->SetOptStat(0);
    //gStyle->SetOptTitle(0);
    
    ///////////////////////////////////////////////////////
    
    NWRITTEN = 0;
    NTOT = 0;
    
    
    ////// Open WaveCatcher *.bin file ////////////////////
    
    int FEXTERN = 0;
    if( sname == "" ) FEXTERN = 0; else FEXTERN = 1;
    
    if( !FEXTERN ) {
        
        //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF1_Data_7_21_2020_Binary.bin";
        
        //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF2_Data_7_31_2020_Binary.bin";
        
        //sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF3_Data_8_6_2020_Binary.bin";
        
        sname="Run_bar220cm_110cm_40ch_2sides_vert_111V_TOF4_Data_8_12_2020_Binary.bin";
        
    }
    
    
    OpenDataFileWC( sname );
    
   
    
    //////// Open output file ///////////////////////
    
    int idot = sname.Index(".");
    TString fname = sname(0,idot);
    fname += ".root";
    
    if( fname == sname ) { printf("WError: the same name\n"); exit(0);}
    
    TFile *f=(TFile*)gROOT->GetListOfFiles()->FindObject(fname.Data());
    if(f==0) f = new TFile(fname.Data(),"RECREATE");
    if( !f->IsOpen() ) {printf("No file %s\n",fname.Data()); exit(0);}
    printf("Output root file %s\n",fname.Data());
    printf("\n");
    
    
    /////////////////////////////////////////////////////////////
    
    int nch = N_CHN;
    int ICH[N_CHN_MAX]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
                        21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39};
    double BL[N_CHN_MAX],RMS[N_CHN_MAX],Max[N_CHN_MAX],Min[N_CHN_MAX];
    double AtrMin=0;
    double tLE10[N_CHN_MAX],tLE20[N_CHN_MAX],tLE30[N_CHN_MAX];
    double tToT10[N_CHN_MAX],tToT20[N_CHN_MAX],tToT30[N_CHN_MAX];
    //double tToTLE10[N_CHN_MAX];
    int iev;
    
  
    
    TTree *tree = new TTree("tree","uToF Events");
    
    tree->Branch("iev" ,&iev ,"iev/I");
    tree->Branch("nch" ,&nch ,"nch/I");
    tree->Branch("ich" ,&ICH ,"ich[nch]/I");
    
    tree->Branch("A" ,&Amp ,"A[nch]/D");
    tree->Branch("AtrMin" ,&AtrMin ,"AtrMin/D");
    
    tree->Branch("tbin" ,&tbin ,"tbin[nch]/D");
    
    tree->Branch("t6" ,&t6 ,"t6[nch]/D");
    tree->Branch("t8" ,&t8 ,"t8[nch]/D");
    tree->Branch("t10",&t10,"t10[nch]/D");
    tree->Branch("t12",&t12,"t12[nch]/D");
    tree->Branch("t20",&t20,"t20[nch]/D");
    tree->Branch("t30",&t30,"t30[nch]/D");
    
    tree->Branch("tLE10" ,&tLE10 ,"tLE10[nch]/D");
    tree->Branch("tLE20" ,&tLE20 ,"tLE20[nch]/D");
    tree->Branch("tLE30" ,&tLE30 ,"tLE30[nch]/D");
    
    tree->Branch("tToT10" ,&tToT10 ,"tToT10[nch]/D");
    tree->Branch("tToT20" ,&tToT20 ,"tToT20[nch]/D");
    tree->Branch("tToT30" ,&tToT30 ,"tToT30[nch]/D");
    
    //tree->Branch("tToTLE10" ,&tToTLE10 ,"tToTLE10[nch]/D");
    
    tree->Branch("tr6" ,&tr6 ,"tr6/D");
    tree->Branch("tr8" ,&tr8 ,"tr8/D");
    tree->Branch("tr10",&tr10,"tr10/D");
    tree->Branch("tr12",&tr12,"tr12/D");
    tree->Branch("tr20",&tr20,"tr20/D");
    tree->Branch("tr30",&tr30,"tr30/D");
    
    tree->Branch("dtr6" ,&dtr6 ,"dtr6/D");
    tree->Branch("dtr8" ,&dtr8 ,"dtr8/D");
    tree->Branch("dtr10",&dtr10,"dtr10/D");
    tree->Branch("dtr12",&dtr12,"dtr12/D");
    tree->Branch("dtr20",&dtr20,"dtr20/D");
    tree->Branch("dtr30",&dtr30,"dtr30/D");
    
    tree->Branch("BL" ,&BL ,"BL[nch]/D"); // baseline
    tree->Branch("RMS",&RMS,"RMS[nch]/D"); // rms of noise
    
    tree->Branch("rt" ,&rt ,"rt[nch]/D"); // rise time
    
    tree->Branch("Max",&Max,"Max[nch]/D");
    tree->Branch("Min",&Min,"Min[nch]/D");
    
    
    
    /////////////////////////////////////////////////////////////
    ///////////// loop over events //////////////////////////////
    
    for( int i=0; ; i++ ) {
    //for(int i=0;i<1;i++) {
        
        if( !ReadDataEventWC() ) break;
        NTOT++;
        iev = i;
        
        
        ////////////////////////////////////////////////
        
        //printf("----- %i %f -----\n",i,EpochTimeWC);
        if( i%1000==0 )
            printf("%5i  Date: %4i %2i %2i   Time: %2i %2i %2i\n",i, YearWC, MonthWC,DayWC, HourWC,MinuteWC,SecondWC);
    
        
        ////////////////////////////////////////////////
        
        double t6_[N_CHN_MAX],t8_[N_CHN_MAX],t10_[N_CHN_MAX],t12_[N_CHN_MAX], t20_[N_CHN_MAX], t30_[N_CHN_MAX];
        
        int im;
        double dtmp;
        
        for(int ich=0; ich<N_CHN; ich++) {
            
            tbin[ich] = TimeBinWC;
           
            
            /////////////////////////////////////////////////////////////
            
            //baseline[ich] = 0;
            GetBaseLine(ich, 40.,70., tbin[ich], baseline[ich],RMS[ich]);
            
            BL[ich]=baseline[ich];
            
            //printf("%2i %f\n",ich, baseline[ich]);
            
            GetMinimum( ich, im, Min[ich]);
            GetMaximum( ich, im, Maximum[ich]);
            Max[ich]=Maximum[ich];
            
            //printf("old: %i %f\n",im, Maximum[ich]);
           
            GetTimeCFD( ich, 0.06, baseline[ich], tbin[ich], t6 [ich],dtmp );
            GetTimeCFD( ich, 0.08, baseline[ich], tbin[ich], t8 [ich],dtmp );
            GetTimeCFD( ich, 0.10, baseline[ich], tbin[ich], t10[ich],dtmp );
            GetTimeCFD( ich, 0.12, baseline[ich], tbin[ich], t12[ich],dtmp );
            GetTimeCFD( ich, 0.20, baseline[ich], tbin[ich], t20[ich],dtmp );
            GetTimeCFD( ich, 0.30, baseline[ich], tbin[ich], t30[ich],dtmp );
            
            //if( t6[ich] != t6_[ich] ) {printf("%4i %i %.2f %.2f\n\n",i,ich,t6[ich],t6_[ich]);}
            //exit(1);
            
            GetTimeLE( ich, 0.01, baseline[ich], tbin[ich], tLE10[ich],tToT10[ich] );
            GetTimeLE( ich, 0.02, baseline[ich], tbin[ich], tLE20[ich],tToT20[ich] );
            GetTimeLE( ich, 0.03, baseline[ich], tbin[ich], tLE30[ich],tToT30[ich] );
            
            rt[ich] = GetRiseTime(ich,tbin[ich]);
            
            //printf("%i %f\n",ich,baseline[ich]);
            //printf("%i %f\n",ich,Maximum[ich]);
            //printf("%i %f\n",ich,t[ich]);
            
            //Amp[ich] = Maximum[ich];
            Amp[ich] = Maximum[ich] - baseline[ich];
            
            //printf("%i %f %f %f\n",im, Maximum[ich],baseline[ich],Amp[ich]);
            
        } // end of loop over channels
        
        
        //for(int ich=0; ich<N_CHN; ich++) printf("%6.2f ",t8[ich]); cout<<endl; exit(1);
        
        ////////////////////////////////////////////////
        
        tree->Fill();
        NWRITTEN++;
        
        //printf("--------------\n");
        
    } // end of loop over events
    
    
    //////////// close file /////////////////////////////////
    
    printf("\n");
    CloseDataFileWC(); // close input file
    
    printf("Statistics:  NTOT = %i,  NWRITTEN = %i,  %.1f %%\n",NTOT,NWRITTEN, 100.*NWRITTEN/NTOT );
    printf("\n");
    
    f->cd();
    f->Write("",TObject::kOverwrite);
    f->Close();  // close output file
    
    ///////////////////////////////////////////
    
    //if( !FEXTERN ) exit(1);
    return;
}








//////////// calculate time with dCFD /////////////

void GetTimeCFD(int ich,double thr,double bl, double tbin, double&t,double&A)
{
    int iMax=-999;
    double max=-999;
    
    GetMaximum(ich, iMax, max);
    //printf("new: %i %f\n",iMax, max);
    
    if(iMax>DATA_SAMPLES) {printf("ERROR1 in GetTimeCFD\n"); exit(0);}
    if(iMax<0           ) {printf("ERROR2 in GetTimeCFD\n"); exit(0);}
    
    double ythr = bl + thr*(max-bl);
    t=0;
    double y1=99999,y2=0;
    
    for(int i=iMax;i>0;i--) {
        //y2 = -(hit.DataSamplesRaw[i]-BaselineSAM);
        y2 = waveform[ich][i];
        if( y2 < ythr ) {
            double x1 = (i+1)*tbin;
            double x2 =     i*tbin;
            double a = (y1-y2)/(x1-x2);
            double b = y1-a*x1;
            t = (ythr-b)/a;
            //printf("new: x1=%.2f, x2=%.2f,  y1=%.4f, %.4f, y2=%.4f\n",x1,x2, y1,ythr,y2);
            //printf("new: ich=%i, x2=%.2f, x1=%.2f,  y2=%.4f, %.4f, y1=%.4f\n",ich, x2,x1, y2,ythr,y1);
            break;
        }
        y1=y2;
    }
    
    //t += hit.Cell0TimeStamp;
    A = ythr;
    
    return;
}



void GetTimeLE(int ich,double thr, double bl, double tbin, double&t , double& ToT )
{
    ToT=0;
    
    bl = 0;
    double ythr = bl+thr;
    t=0;
    double y1=0,y2=0,y3=0;
    
    int i=0;
    
    for(i=5; i<DATA_SAMPLES; i++ ) {
        //y2 = -(hit.DataSamplesRaw[i]-BaselineSAM);
        y2 = waveform[ich][i];
        y3 = waveform[ich][i+1];
        if( y2 > ythr && y3 > y2 ) {
            double x1 = (i-1)*tbin;
            double x2 =     i*tbin;
            double a = (y1-y2)/(x1-x2);
            double b = y1-a*x1;
            t = (ythr-b)/a;
            //printf("y1=%f, y2=%f\n",y1,y2);
            break;
        }
        y1=y2;
    }
    
    for( i=i+2; i<DATA_SAMPLES; i++ ) {
        //y2 = -(hit.DataSamplesRaw[i]-BaselineSAM);
        y2 = waveform[ich][i];
        
        if( y2 < ythr ) {
            double x1 = (i-1)*tbin;
            double x2 =     i*tbin;
            double a = (y1-y2)/(x1-x2);
            double b = y1-a*x1;
            double t2 = (ythr-b)/a;
            ToT = t2-t;
            //printf("y1=%f, y2=%f\n",y1,y2);
            break;
        }
        y1=y2;
    }
    
    
    //printf("t = %5.1f,  ToT = %4.1f\n",t,ToT);
    
    return;
}



void GetMaximum(int ich, int& im, double& max)
{
    //double X[DATA_SAMPLES];
    //for(int i=0;i<DATA_SAMPLES;i++) X[i] = i*TimeBinWC;
    
    im = -99; max = -1e99;
    for(int i=0;i<DATA_SAMPLES;i++) {
        
        //if(X[i]<70||X[i]>150) continue;
        
        double y = waveform[ich][i];
        if(y>max){max=y;im=i;}
    }
    return;
}


void GetMinimum(int ich, int& im, double& min)
{
    im = -99; min = 1e99;
    for(int i=0;i<DATA_SAMPLES;i++) {
        double y = waveform[ich][i];
        if(y<min){min=y;im=i;}
    }
    return;
}


void GetBaseLine(int ich, double tl,double tr, double tbin, double& bl, double& rms)
{
    double X[DATA_SAMPLES];
    for(int i=0;i<DATA_SAMPLES;i++) X[i] = i*tbin;
    
    unsigned int n = 0;
    double sum = 0;
    
    for(int i=0;i<DATA_SAMPLES;i++) {
        if(X[i]<tl) continue;
        if(X[i]>tr) break;
        //printf("%f \n", X[i]);
        //cout<<X[i]<<endl;
        sum += waveform[ich][i];
        n++;
    }
    
    //cout<<tbin<<endl;
    //exit(1);
    
    bl = sum/n;
    
    /////// RMS
    
    sum=0; n=0;
    
    for(int i=0; i<DATA_SAMPLES; i++ ) {
        if(X[i]<tl) continue;
        if(X[i]>tr) break;
        sum += pow( bl-waveform[ich][i], 2 );
        n++;
    }
    
    rms = sum/n;
    rms = sqrt(rms);
    
    return;
}



double GetRiseTime(int ich,double tbin)
{
    double yB=0.1, yT=0.9;
    double tB=0, tT=0;
    
    double vMax = -99999, vMin = 0;
    int iMax = -99999;
    int iFirst = 6;
    
    GetMaximum(ich, iMax, vMax);
    
    double yMin = vMin + yB*(vMax-vMin);
    double yMax = vMin + yT*(vMax-vMin);
    //printf("yMin=%.3f V, yMax=%.3f V\n",yMin,yMax);
    
    
    ////////////////////////////////////////////////
    double X[DATA_SAMPLES];
    for(int i=0;i<DATA_SAMPLES;i++) X[i] = i*tbin;
    
    
    ///////////////////////////////////////////////////
    ///////////// new, right -> left ///////////
    int i=0;
    int iTop = 0;

    for( i = iMax; i>iFirst; i-- ) {
      if( waveform[ich][i] < yMax && waveform[ich][i+1] >= yMax ) {
        double a = (waveform[ich][i]-waveform[ich][i+1]) / (X[i]-X[i+1]);
        double b = waveform[ich][i] - a*X[i];
        tT = (yMax - b) / a;
        //printf("Rise Time(1): i=%i, tT = %.2f ns\n",i,tT);
        iTop = i;
        //break;
      }
    }

    for( i = iTop; i>iFirst; i-- ) {
      //printf("",);
      if( waveform[ich][i] < yMin ) {
        double a = (waveform[ich][i]-waveform[ich][i+1]) / (X[i]-X[i+1]);
        double b =  waveform[ich][i] - a*X[i];
        tB = (yMin - b) / a;
        //printf("Rise Time(2): i=%i, tB = %.2f ns\n",i,tB);
        break;
      }
    }
    
    double dt = tT - tB;
    
    //exit(0);
    return dt;
}




void ResetSamplingFrequency(int ich,int ISFS)
{
    double wftmp[N_CHN_MAX][DATA_SAMPLES];
    
    for(int i=0; i<DATA_SAMPLES; i++) {
        wftmp[ich][i] = waveform[ich][i];
        waveform[ich][i] = 0;
    }
    
    /////////////////////////////
    
    for(int i=0; i<DATA_SAMPLES/ISFS; i++) waveform[ich][i] = wftmp[ich][ISFS*i];
    
    return;
}
