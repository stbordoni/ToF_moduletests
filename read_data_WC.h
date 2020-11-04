#ifndef read_WC_h
#define read_WC_h

#include <stdio.h>
#include <stdlib.h>
//#include <iostream>


//#include "/Users/alex/work/DRS4/C_scripts/read_data_arr.h"
long int N_EVENT = 0;
int N_TRIGGER = 0;
int N_BOARD = 0;
int N_CHN = 0;
//const int N_CHN_MAX = 18;
const int N_CHN_MAX = 48;
double waveform[N_CHN_MAX][1024], timeAr[N_CHN_MAX][1024];
double fr_sampl = 0;



/*
long int N_EVENT = 0;
int N_CHN = 0;
const int N_CHN_MAX = 100;
double waveform[N_CHN_MAX][1024], timeAr[N_CHN_MAX][1024];
double fr_sampl = 0;
*/

FILE *DataFileWC = 0, *DataOutFileWC = 0;

char SoftVersWC[10]; // SOFTWARE VERSION
double GainWC = 0;
int N_CHN_MAX_WC = 0; // max number of channels in the WC module
double TimeBinWC = 0; // ns
int N_EVENTWC = 0;
double EpochTimeWC = 0;
unsigned long int TDCperiodWC;
unsigned int YearWC=0, MonthWC=0,DayWC=0,HourWC=0,MinuteWC=0,SecondWC=0,MillisecondWC=0;

int ichannelWC[N_CHN_MAX];
int ieventWC[N_CHN_MAX];
int icellWC[N_CHN_MAX];

float baselineWC[N_CHN_MAX];
float amplitudeWC[N_CHN_MAX];
float chargeWC[N_CHN_MAX];
float risetimeWC[N_CHN_MAX];
float falltimeWC[N_CHN_MAX];
float ratecounterWC[N_CHN_MAX];

int CHANNEL_ACC[N_CHN_MAX];

////////////////////////////////

void OpenDataFileWC(const char* cname);
void OpenDataOutFileWC(const char* cname);
bool ReadDataEventWC();
void WriteDataEventWC();
void CloseDataFileWC();
void CloseDataOutFileWC();

////////////////////////////////


void read_data_WC()
{
  //OpenDataFileWC("Run_LG_bar50cm_25cm_55V_h5GeV_Data_5_5_2017_Binary.bin");
  OpenDataFileWC("in.bin");

  // set channels to write
  for(int ich=0; ich<N_CHN_MAX; ich++ ) CHANNEL_ACC[ich] = 0;
  CHANNEL_ACC[0] = 1;
  CHANNEL_ACC[1] = 1;
  CHANNEL_ACC[2] = 1;
  CHANNEL_ACC[3] = 1;

  OpenDataOutFileWC("out.bin");

  //while( ReadDataEventWC() );
  for(int i=0;i<10;i++) {
    ReadDataEventWC();
    //printf("\n");
    WriteDataEventWC();
    //printf("--------------\n");
  }

  CloseDataFileWC();
  CloseDataOutFileWC();

  return;
}



void OpenDataOutFileWC(const char* cname)
{
  //////////// open binary file ///////////////////////

  if( DataOutFileWC != 0 ) fclose(DataOutFileWC);

  DataOutFileWC = fopen(cname, "w");
  if ( DataOutFileWC == 0 ) {
    printf("Cannot write data file %s\n", cname );
    exit(0);
  } else {
    printf("------- Writing to file %s\n", cname );
  }
  printf("\n");


  ////////// write ASCII header //////////////////////////

  fprintf(DataOutFileWC,"=== DATA FILE SAVED WITH SOFTWARE VERSION: %s ===\n",SoftVersWC);

  fprintf(DataOutFileWC,"=== WAVECATCHER SYSTEM OF TYPE 1 WITH 18 CHANNELS AND GAIN: %.1f ===\n",GainWC);

  //fprintf(DataOutFileWC,"=== Rate coincidence masks (4 bits) for ch 0 to 16: 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 === ");
  fprintf(DataOutFileWC,"=== Rate coincidence masks (4 bits) for ch 0 to 15: 1 2 4 8 1 2 4 8 1 2 4 8 1 2 4 8 === ");

  int N_CHN_ACC = 0;
  for(int ich=0; ich<N_CHN; ich++ ) if( CHANNEL_ACC[ich] == 1 ) N_CHN_ACC++;
  printf("Writing: N_CHN_ACC = %i\n",N_CHN_ACC);

  // !!
  fprintf(DataOutFileWC,"Posttrig in ns for SamBlock 0 to 8: 140 140 140 140 140 140 140 140 0 ===\n");

  fprintf(DataOutFileWC,"=== DATA SAMPLES [1024] in Volts == NB OF CHANNELS ACQUIRED: %i == Sampling Period: %.1f ps == INL Correction: 1 ===\n",N_CHN_ACC,1000*TimeBinWC);

  ///////////////////////////////

  printf("%i channels to write: ", N_CHN_ACC);
  for(int ich=0; ich<N_CHN; ich++ ) printf("%i ",CHANNEL_ACC[ich]);
  printf("\n");

  ///////////////////////////////
  //exit(1);

  return;
}





/////////////////////////////////////////////////

void OpenDataFileWC(const char* cname)
{
  printf("\n");

  for(int ich=0; ich<N_CHN_MAX; ich++ ) CHANNEL_ACC[ich] = 1;
  
  /////////////////////

  N_EVENT = 0;

  //////////// open binary file ///////////////////////

  if( DataFileWC != 0 ) CloseDataFileWC();

  DataFileWC = fopen(cname, "r");
  if ( DataFileWC == 0 ) {
    printf("Cannot find data file %s\n", cname );
    exit(0);
  } else {
    printf("------- Reading WAVECATCHER file %s\n", cname );
  }
  printf("\n");

  ////////// read ASCII header //////////////////////////

  //char *line = NULL;
  //size_t len = 0;
  //ssize_t iread=0;
  //iread = getline(&line, &len, DataFileWC);
  //printf("%s\n",line);
  //free(line);

  int items=0;
  char ctmp[100];

  for(int i=0;i<7;i++) items=fscanf(DataFileWC,"%s",ctmp);
  items=fscanf(DataFileWC,"%s",SoftVersWC);

  for(int i=0;i<12;i++) items=fscanf(DataFileWC,"%s",ctmp);
  items=fscanf(DataFileWC,"%lf",&GainWC);

  for(int i=0;i<11;i++) items=fscanf(DataFileWC,"%s",ctmp);
  //printf("%s\n",ctmp);
  items=fscanf(DataFileWC,"%i",&N_CHN_MAX_WC);
  //printf("%i\n",N_CHN_MAX_WC);
  N_CHN_MAX_WC++;
  //printf("%i\n",N_CHN_MAX_WC);
  //exit(1);

    
  for(int i=0;i<N_CHN_MAX_WC;i++) items=fscanf(DataFileWC,"%s",ctmp);

  for(int i=0;i<10-1;i++) items=fscanf(DataFileWC,"%s",ctmp);
  int SamBlock=0;
  fscanf(DataFileWC,"%i",&SamBlock);
  //printf("SamBlock = %i\n",SamBlock); //exit(1); // SamBlock
  
    
  //for(int i=0;i<N_CHN_MAX_WC/2;i++) items=fscanf(DataFileWC,"%s",ctmp);
  for(int i=0;i<SamBlock+1;i++) items=fscanf(DataFileWC,"%s",ctmp);
  //printf("%s\n",ctmp); //exit(1);

  for(int i=0;i<13;i++) items=fscanf(DataFileWC,"%s",ctmp);
  //printf("%s\n",ctmp); //exit(1);  // ACQUIRED:

  items=fscanf(DataFileWC,"%i",&N_CHN);
  //printf("%i\n",N_CHN); //exit(1);


  for(int i=0;i<3;i++) items=fscanf(DataFileWC,"%s",ctmp);
  //printf("%s\n",ctmp);
  items=fscanf(DataFileWC,"%lf",&TimeBinWC);

  for(int i=0;i<6;i++) items=fscanf(DataFileWC,"%s",ctmp);
  //printf("%s\n",ctmp);

  int nitem = 0;
  nitem = fread(ctmp, 1, 1, DataFileWC); // remove 1 byte

  printf("SOFTWARE VERSION = %s\n",SoftVersWC);
  printf("GAIN = %.1f\n",GainWC);
  printf("Max Numb of CHANNELS = %i\n",N_CHN_MAX_WC);
  printf("CHANNELS ACQUIRED = %i\n",N_CHN);
  //if(N_CHN>8) N_CHN=8; // !!!
  //N_CHN -= 2; // !!!
  //printf("CHANNELS ACQUIRED (by hand!!!) = %i\n",N_CHN);

  TimeBinWC /= 1000; // ps -> ns
  printf("Sampling bin = %.5f ns\n",TimeBinWC);
  printf("Time window = %.1f ns\n",TimeBinWC*1024);
  fr_sampl = 1/TimeBinWC;
  printf("Sampling frequency = %.3f GHz\n",fr_sampl);
  printf("\n");

  //exit(0);

  //////////////////////////////////////
  
  if( fr_sampl==0 || isnan(fr_sampl) ) {
    printf("Something is wrong with the sampling frequency!!!\n");
    exit(0);
  }
  
  ///////// set time bins //////////////////////////////////

  for( int ich=0; ich<N_CHN; ich++ ) {
    for(int icell=0; icell<1024; icell++ ) {
      timeAr[ich][icell] = icell*TimeBinWC;
    }
  }

  //////////////////////////////////////////////////////////
  printf("\n");
  //exit(1);

  return;
}




void WriteDataEventWC()
{
  int nitem=0;

  nitem = fwrite(&N_EVENTWC,  sizeof(int),   1, DataOutFileWC);

  nitem = fwrite(&EpochTimeWC  , sizeof(double)      , 1, DataOutFileWC);
  nitem = fwrite(&YearWC       , sizeof(unsigned int), 1, DataOutFileWC);
  nitem = fwrite(&MonthWC      , sizeof(unsigned int), 1, DataOutFileWC);
  nitem = fwrite(&DayWC        , sizeof(unsigned int), 1, DataOutFileWC);
  nitem = fwrite(&HourWC       , sizeof(unsigned int), 1, DataOutFileWC);
  nitem = fwrite(&MinuteWC     , sizeof(unsigned int), 1, DataOutFileWC);
  nitem = fwrite(&SecondWC     , sizeof(unsigned int), 1, DataOutFileWC);
  nitem = fwrite(&MillisecondWC, sizeof(unsigned int), 1, DataOutFileWC);

  ///////////////////////

  if( SoftVersWC[4] != '7' ) {

    ///printf("Writing...");

    nitem = fwrite(&TDCperiodWC, sizeof(unsigned long int), 1, DataOutFileWC);
    //printf("%li\n",TDCperiodWC);

    unsigned int iich=0;
    nitem = fwrite(&iich, sizeof(unsigned int), 1, DataOutFileWC);
    //printf("%i %i\n",(int)ich,(int)sizeof(unsigned int));
    //printf("\n");
  }

  //int ichannel=0;
  //int ievent =0;
  //int icell =0;

  for(int ich=0; ich<N_CHN; ich++) {
  //for(int ich=0; ich<N_CHN_MAX; ich++) {

    if( CHANNEL_ACC[ich] == 0 ) continue;

    nitem = fwrite(&ichannelWC[ich], sizeof(int), 1, DataOutFileWC);
    nitem = fwrite(&ieventWC[ich],   sizeof(int), 1, DataOutFileWC);
    nitem = fwrite(&icellWC[ich],    sizeof(int), 1, DataOutFileWC);
    //printf("ichannel=%i, ievent=%i, icell=%i\n",ichannelWC[ich],ieventWC[ich],icellWC[ich]);

    nitem = fwrite(&baselineWC[ich]   , sizeof(float), 1, DataOutFileWC);
    nitem = fwrite(&amplitudeWC[ich]  , sizeof(float), 1, DataOutFileWC);
    nitem = fwrite(&chargeWC[ich]     , sizeof(float), 1, DataOutFileWC);
    nitem = fwrite(&risetimeWC[ich]   , sizeof(float), 1, DataOutFileWC);
    nitem = fwrite(&falltimeWC[ich]   , sizeof(float), 1, DataOutFileWC);
    nitem = fwrite(&ratecounterWC[ich], sizeof(float), 1, DataOutFileWC);



    signed short itmp=0;
    double coef = 2.5 / (4096 * 10); // 2^12=4096
    for(int icell=0; icell<1024; icell++) {
      //nitem = fwrite(&itmp, sizeof(signed short), 1, DataOutFileWC);
      //printf(" %i %i\n ",icell,itmp);
      //waveform[ich][icell] = coef*itmp;
      //printf(" %4i %9f\n ",icell,waveform[ich][icell]);
      itmp = waveform[ich][icell] / coef;
      nitem = fwrite(&itmp, sizeof(signed short), 1, DataOutFileWC);
    }

    //exit(0);
  }

  return;
}



bool ReadDataEventWC()
{
  //printf("ReadDataEventWC()\n");

  int nitem = 0;

  nitem = fread(&N_EVENTWC, sizeof(int), 1, DataFileWC);
  if( nitem == 0 ) {
    printf("ReadDataEventWC() ERROR: nitem = %i,  sizeof(th)=%li\n",nitem, sizeof(N_EVENTWC) );
    //exit(0);
    return false;
  }
  //printf("N_EVENTWC = %i\n",N_EVENTWC);


  nitem = fread(&EpochTimeWC  , sizeof(double)      , 1, DataFileWC);
  nitem = fread(&YearWC       , sizeof(unsigned int), 1, DataFileWC);
  nitem = fread(&MonthWC      , sizeof(unsigned int), 1, DataFileWC);
  nitem = fread(&DayWC        , sizeof(unsigned int), 1, DataFileWC);
  nitem = fread(&HourWC       , sizeof(unsigned int), 1, DataFileWC);
  nitem = fread(&MinuteWC     , sizeof(unsigned int), 1, DataFileWC);
  nitem = fread(&SecondWC     , sizeof(unsigned int), 1, DataFileWC);
  nitem = fread(&MillisecondWC, sizeof(unsigned int), 1, DataFileWC);

  /*
  printf("EpochTimeWC = %f\n",EpochTimeWC);
  printf("YearWC = %i\n",YearWC);
  printf("MonthWC = %i\n",MonthWC);
  printf("DayWC = %i\n",DayWC);
  printf("HourWC = %i\n",HourWC);
  printf("MinuteWC = %i\n",MinuteWC);
  printf("SecondWC = %i\n",SecondWC);
  printf("MillisecondWC = %i\n",MillisecondWC);
  //printf("\n");
  //exit(1);
  */

  if( SoftVersWC[4] != '7' ) {
    nitem = fread(&TDCperiodWC, sizeof(unsigned long int), 1, DataFileWC);
    //printf("%li\n",TDCperiodWC);

    unsigned int ich;
    nitem = fread(&ich, sizeof(unsigned int), 1, DataFileWC);
    //printf("%i %i\n",(int)ich,(int)sizeof(unsigned int));
    //printf("\n");
  }



  ///////////////////////

  //int ichannel=0;
  //int ievent =0;
  //int icell =0;

  for(int ich=0; ich<N_CHN; ich++) {
    nitem = fread(&ichannelWC[ich], sizeof(int), 1, DataFileWC);
    nitem = fread(&ieventWC[ich],   sizeof(int), 1, DataFileWC);
    nitem = fread(&icellWC[ich],    sizeof(int), 1, DataFileWC);
    //printf("ichannel=%i, ievent=%i, icell=%i\n",ichannelWC[ich],ieventWC[ich],icellWC[ich]);

    nitem = fread(&baselineWC[ich]   , sizeof(float), 1, DataFileWC);
    nitem = fread(&amplitudeWC[ich]  , sizeof(float), 1, DataFileWC);
    nitem = fread(&chargeWC[ich]     , sizeof(float), 1, DataFileWC);
    nitem = fread(&risetimeWC[ich]   , sizeof(float), 1, DataFileWC);
    nitem = fread(&falltimeWC[ich]   , sizeof(float), 1, DataFileWC);
    nitem = fread(&ratecounterWC[ich], sizeof(float), 1, DataFileWC);

    //printf("%2i  baseline=%f, amplitude=%f, charge=%f\n",ich,baselineWC[ich],amplitudeWC[ich],chargeWC[ich]);
    //printf(" risetime=%f, falltime=%f, ratecounter=%f\n",risetimeWC[ich],falltimeWC[ich],ratecounterWC[ich]);

    signed short itmp=0;
    double coef = 2.5 / (4096 * 10); // 2^12=4096
    for(int icell=0; icell<1024; icell++) {
      nitem = fread(&itmp, sizeof(signed short), 1, DataFileWC);
      //printf(" %i %i\n ",icell,itmp);
      waveform[ich][icell] = coef*itmp;
      //printf("%4i ",itmp);
      //printf(" %4i %9f\n ",icell,waveform[ich][icell]);

    }
    //printf("\n");

    //exit(0);
  }
  //exit(0);

  ///////////////////////

  N_EVENT++;

  return true;
}




void CloseDataFileWC()
{
  printf("CloseDataFileWC()\n\n");
  fclose(DataFileWC);
  DataFileWC = 0;

  printf("N_EVENT = %li, N_EVENTWC = %i\n",N_EVENT,N_EVENTWC);
  printf("\n");

  return;
}


void CloseDataOutFileWC()
{
  printf("CloseDataOutFileWC()\n\n");
  fclose(DataOutFileWC);
  DataOutFileWC = 0;

  //printf("N_EVENT = %li, N_EVENTWC = %i\n",N_EVENT,N_EVENTWC);
  //printf("\n");

  return;
}


#endif
