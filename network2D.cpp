//Rulkov, N., Timofeev, I. and Bazhenov, M. "Oscillations in large-scale cortical 
//networks: map-based model" (2004) Journal of Computational Neuroscience 17, 203-224, 2004.

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <iostream.h> 

//---------------Parameters initialization--------------------------------------------------- 
//--------------- Network Geometry (2-dimensional 2-layer network) -------------------------- 
#define I_CX    1     //  0 - No layer, 1 - Add Layer of Pyramidal neurons 
#define I_IN    1     //  0 - No layer, 1 - Add Layer of Interneurons
 
#define Mcx      256      //  Number of CX cells (Pyramidal neurons) in X direction
#define Min      128      //  Number of IN cells (Interneurons) in X direction
#define Mcx1     256      //  Number of CX cells in Y direction 
#define Min1     128      //  Number of IN cells in Y direction  
#define Max     ((M > Mcx) ? M : Mcx) 
#define Max1    ((M1 > Mcx1) ? M1 : Mcx1) 
 
#define SELFINGcx  0   //  0 - CX without self excitation; 1 - with ...  
 
//-------------- Define the connections between cells ----------------- 
 
#define MS_CX_CX 8   //radius of CX->CX connections in X direction
#define MS_CX_CX1 8  //radius of CX->CX connections in Y direction
#define MS_CX_CX_MAX  ((MS_CX_CX > MS_CX_CX1) ? MS_CX_CX : MS_CX_CX1)
#define N_CX_CX  (2*MS_CX_CX+1)*(2*MS_CX_CX1+1) //total number of CX connections accepted by one CX

#define MS_CX_IN 8   //radius of CX->IN connections in X direction
#define MS_CX_IN1 8  //radius of CX->IN connections in Y direction
#define MS_CX_IN_MAX  ((MS_CX_IN > MS_CX_IN1) ? MS_CX_IN : MS_CX_IN1)
#define N_CX_IN  (2*MS_CX_IN+1)*(2*MS_CX_IN1+1) //total number of CX connections accepted by one IN

#define MS_IN_CX 2   //radius of IN->CX connections in X direction
#define MS_IN_CX1 2  //radius of IN->CX connections in Y direction
#define MS_IN_CX_MAX  ((MS_IN_CX > MS_IN_CX1) ? MS_IN_CX : MS_IN_CX1) 
#define N_IN_CX  (2*MS_IN_CX+1)*(2*MS_IN_CX1+1)
 
//-------------------RS CELL (main class to describe regular spiking neuron)-----------------------
class RS{
  double xp, xpp, mu, sigma_n, beta_n;
  double sigma_dc, beta_dc;
public:
  double x, y, alpha, sigma, sigma_e, beta_e, Idc, S_CX_DEND;
  int spike;

  RS(){ 
	mu=0.0005;
	spike=0;
	alpha=3.65;
	sigma=6.00E-2;
	sigma_e=1.0;
        sigma_dc=1.0;
	beta_e=0.133;
        beta_dc=0.133;
        Idc = 0;
        S_CX_DEND = 165.0e-6;
        xpp=-1+sigma;
        xp=-1+sigma;
        x=-1+sigma;
        y= x - alpha/(1-x);
        }
 
  void init(){ 
	mu=0.0005;
	spike=0;
	alpha=3.65;
	sigma=6.00E-2;
	sigma_e=1.0;
        sigma_dc=1.0;
	beta_e=0.133;
        beta_dc=0.133;
        Idc = 0;
        S_CX_DEND = 165.0e-6;
        xpp=-1+sigma;
        xp=-1+sigma;
        x=-1+sigma;
        y= x - alpha/(1-x);
        }
  void calc(double); 
};  

void RS::calc(double I){
     beta_n = beta_e * I + beta_dc * Idc;
     sigma_n = sigma_e * I + sigma_dc * Idc;

     if (beta_n < -0.0001) beta_n = -0.0001;
     if (beta_n > 1.0) beta_n = 1.0;

     if(xp <= 0.0) { 
        x = alpha / (1.0 - xp) + y +beta_n;
	spike = 0;
     }
     else{
     if(xp <= alpha + y +beta_n && xpp <= 0.0) { 
	x = alpha + y + beta_n;
	spike = 1;
	} 
	else {
	x = -1;
	spike = 0;
	}
     }
     
       	y = y - mu* (xp +1.0) + mu * sigma + mu * sigma_n;
       	xpp = xp;
	xp = x;
	//     y=1;
}


//-------------------FS1 CELL (main class to describe fast spiking interneuron)-----------------
class FS1{ 
  double xp, xpp, mu, sigma_n, beta_n, ii, gg; 
  double sigma_dc, beta_dc, beta_ii, sigma_ii; 
public: 
  double x, y, alpha, sigma, sigma_e, beta_e, Idc, S_CX_DEND; 
  int spike; 
 
  FS1(){  
        ii=0;
	mu=0.002;
	spike=0; 
	alpha=3.8;     //3.87;
	sigma=-1.75E-2;
	sigma_e=1.0; 
	sigma_dc=1.0; 
	beta_e=0.1;
	beta_dc=0.0;
        Idc = 0; 
        S_CX_DEND = 165.0e-6; 
        xpp=-1+sigma; 
        xp=-1+sigma;
        x=-1+sigma; 
        y= x - alpha/(1-x); 
        gg = 0.5;
        beta_ii = 0.5;
        sigma_ii = 0.5;
        } 
  void init(){  
        ii=0;
	mu=0.002;
	spike=0; 
	alpha=3.8;     //3.87;
	sigma=-1.75E-2;
	sigma_e=1.0; 
	sigma_dc=1.0; 
	beta_e=0.1;
	beta_dc=0.0;
        Idc = 0; 
        S_CX_DEND = 165.0e-6; 
        xpp=-1+sigma; 
        xp=-1+sigma;
        x=-1+sigma; 
        y= x - alpha/(1-x); 
        gg = 0.5;
        beta_ii = 0.5;
        sigma_ii = 0.5;
        } 
  void calc(double);  
};
 
void FS1::calc(double I){ 

     if(spike > 0.1){  
         ii =  0.60 * ii - gg;
       }
       else{
         ii = 0.60 * ii;
       }       

     beta_n = beta_e * I + beta_dc * Idc; 
     sigma_n = sigma_e * I + sigma_dc * Idc; 
     if (beta_n < -0.0001) beta_n = -0.0001;
     if (beta_n > 1.0) beta_n = 1.0;
     beta_n = beta_n + beta_ii *ii;
     sigma_n = sigma_n + sigma_ii *ii;
     if(xp <= 0.0) {  
        x = alpha / (1.0 - xp) + y +beta_n; 
	spike = 0; 
     } 
     else{ 
     if(xp <= alpha + y +beta_n && xpp <= 0.0) {  
	x = alpha + y + beta_n; 
	spike = 1; 
	}  
	else { 
	x = -1; 
	spike = 0; 
	} 
     }

//      	y = y - mu* (xp +1.0) + mu * sigma + mu * sigma_n;
        y = -2.9;
       	xpp = xp;
	xp = x; 
} 

//------------Main class to describe first order kinetic model for AMPA map-based------------ 
//--------------------------------------synapse with depression--------------------------------
class AMPAmapD {
   static double E_AMPA;
   static double gamma;
   double d;   // depression variable
public:
   double I, d_dep, d_rec;
   AMPAmapD() {
       I=0;
       d=1;
       d_dep=0.5;
       d_rec=0.010;
  }
   void calc(double, double, int);
};
double AMPAmapD::E_AMPA = 0, AMPAmapD::gamma = 0.6;

void AMPAmapD::calc(double g_AMPA, double x_post, int spike) {
         if(spike > 0.1){
           I = gamma * I - d*g_AMPA * (x_post - E_AMPA);
           d = (1.0 - d_dep)*d;
        }
        else{
          I = gamma * I;
          d = 1.0 - (1.0 - d_rec) * (1.0 - d);
        }
 }

 
//---------Main class to describe first order kinetic model for AMPA map-based synapse------- 
class AMPAmap { 
  static double E_AMPA; 
  static double gamma; 
 
public: 
  double I; 
  AMPAmap() { 
      I=0; 
  } 
  void calc(double, double, int); 
}; 
double AMPAmap::E_AMPA = 0, AMPAmap::gamma = 0.6; 
void AMPAmap::calc(double g_AMPA, double x_post, int spike) { 
 
       if(spike > 0.1){  
         I = gamma * I - g_AMPA * (x_post - E_AMPA); 
       } 
       else{ 
         I = gamma * I; 
       } 
} 
 
//--------Main class to describe first order kinetic model for GABA-A map-based synapse---------- 
class GABAAmap { 
  static double E_GABAA; 
  static double gamma; 
 
public: 
  double I; 
  GABAAmap() { 
      I=0; 
  } 
  void calc(double, double, int); 
}; 
double GABAAmap::E_GABAA = -1.1, GABAAmap::gamma = 0.8;
void GABAAmap::calc(double g_GABA, double x_post, int spike) { 
 
       if(spike > 0.1){  
         I =  gamma * I - g_GABA * (x_post - E_GABAA); 
       } 
       else{ 
         I = gamma * I; 
       } 
} 
 
//-----Main class to describe first order kinetic model for AMPA synapse used for external---------
//--------------------------- stimulation---------------------------------------------------------- 
class Extern_ampa { 
  static double Cdur, Cmax, Deadtime, Prethresh;  
  double R, C, R0, R1; 
  double lastrelease; 
  double q, Rinf, Rtau; 
  double TR, w, wom, RRR; 
  double exptable(double z) 
    { 
    if((z > -10) && (z < 10)) return( exp(z) ); 
    else return( 0 ); 
    } 
public: 
  double g, Alpha, Beta; 
  Extern_ampa() { 
    Alpha = 0.94; 
    Beta = 0.18; 
    R = 0, C = 0, R0 = 0, R1 = 0; 
    lastrelease = -100; 
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta); 
    Rtau = 1 / ((Alpha * Cmax) + Beta); 
    TR = 1000, w=0.01, wom=0; 
  } 
  void iii(unsigned int seek) {srand(seek);} 
  void calc(double g_Extern_ampa, double x); 
}; 
double Extern_ampa::Cdur = 0.3, Extern_ampa::Cmax = 0.5, Extern_ampa::Deadtime = 1; 
double Extern_ampa::Prethresh = 0; 
void Extern_ampa::calc(double g_Extern_ampa, double x) { 
 
        q = ((x - lastrelease) - Cdur);          
        if (q > Deadtime) { 
                if ((x - lastrelease) > TR) {         
                        C = Cmax;                 
                        R0 = R; 
                        lastrelease = x; 
//                        RRR = 1.0*rand()/(RAND_MAX + 1.0); 
// 	  	  	  if(RRR < 0.000001) RRR = 0.000001; 
//                        TR = -(log(RRR))/(w+(w/2)*sin(wom*x)); 
                } 
        } else if (q < 0) {                      
        } else if (C == Cmax) {                   
                R1 = R; 
                C = 0.; 
        } 
        if (C > 0) {                             
           R = Rinf + (R0 - Rinf) * exptable (-(x - lastrelease) / Rtau); 
        } else {                               
           R = R1 * exptable (-Beta * (x - (lastrelease + Cdur))); 
        } 
       g = g_Extern_ampa * R; 
} 
 
//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++ 
 
//----------external functions---------------------------------------------- 
void fun(double); 
 
//----------external variables --------------------------------------------- 
double *g_ampa_cx_cx[Mcx][Mcx1], *g_ampa_cx_in[Min][Min1], *g_gaba_a_in_cx[Mcx][Mcx1];   
double g_ext_cx, g_ext_in; 
FILE *f1, *f11, *f13, *f98, *f99; 

int C_CXCX[2*MS_CX_CX+1][2*MS_CX_CX1+1];
int C_CXIN[2*MS_CX_IN+1][2*MS_CX_IN1+1]; 
int C_INCX[2*MS_IN_CX+1][2*MS_IN_CX1+1]; 

int k_CXCX, k_CXIN, k_INCX; 
int k_CXCXmax=0, k_CXINmax=0, k_INCXmax=0; 
 
//----------external classes (beginning of initialization)------------------ 

//---Synapses---------------------
AMPAmapD       *a_cx_cx[Mcx][Mcx1]; 
AMPAmap        *a_cx_in[Min][Min1]; 
GABAAmap       *ga_in_cx[Mcx][Mcx1]; 
 
//---Neurons-----------------------
RS            cx_cell[Mcx][Mcx1]; 
FS1           in_cell[Min][Min1];          

//---External input---------------- 
Extern_ampa  a_ext1; 
 
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
main(int argc,char **argv) 
{ 
//---------general parameters---------------------------------------------- 
  double t = 0, tmax, t3D, ttime; 
  double h = 0.02, R, av=0; 
  int i, j, i1, j1, k, ii = 0, jj=0, i_gip = 0; 
  double scale; 
  double g_AMPA_CX_CX, g_AMPA_CX_IN, g_GABA_A_IN_CX; 
  double g_Extern_ampa; 
  double TTR=0, aver=0;
  int TTRS=0;
  char namein[10], namecx[10];
//-------parameter initialization (from file)---------------------------------- 

  f1=fopen("input2D.txt","r");
  fscanf(f1, "%lf %lf %lf %lf %lf %lf %lf ",
       &tmax, &t3D, &ttime, 
       &g_AMPA_CX_CX, &g_AMPA_CX_IN, &g_GABA_A_IN_CX,   
       &g_Extern_ampa);
  printf("\n param: %lf %lf %lf AM_CX_CX=%lf AM_CX_IN=%lf GA_IN_CX=%lf %lf ",
       tmax, t3D, ttime,  
       g_AMPA_CX_CX, g_AMPA_CX_IN, g_GABA_A_IN_CX,   
       g_Extern_ampa); 
  fclose(f1); 
 

//----------classes initialization (continue)---------------------------- 
  printf("\n CLASS INIT CONT"); 
  for(i=0; i < Mcx; ++i) 
  for(j=0; j < Mcx1; ++j){  
    a_cx_cx[i][j] = new AMPAmapD[N_CX_CX]; 
    ga_in_cx[i][j] = new GABAAmap[N_IN_CX];  
    g_ampa_cx_cx[i][j] = new double[N_CX_CX]; 
    g_gaba_a_in_cx[i][j] = new double[N_IN_CX]; 
 
  } 
  for(i=0; i < Min; ++i) 
  for(j=0; j < Min1; ++j){  
    a_cx_in[i][j] = new AMPAmap[N_CX_IN]; 
    g_ampa_cx_in[i][j] = new double[N_CX_IN]; 
 
  } 
 
//--------set synaptic conductances to zero (will be updated later)---------------- 
 
if(I_CX == 1){ 
  for(i=0; i < Mcx; ++i) 
   for(j=0; j < Mcx1; ++j) 
    for(k=0; k < N_CX_CX; ++k) 
      g_ampa_cx_cx[i][j][k] = 0; 
 
  for(i=0; i < Mcx; ++i) 
   for(j=0; j < Mcx1; ++j) 
    for(k=0; k < N_IN_CX; ++k) 
      g_gaba_a_in_cx[i][j][k] = 0; 
  g_ext_cx = 0; 
} 
 
if(I_IN == 1){ 
  for(i=0; i < Min; ++i) 
   for(j=0; j < Min1; ++j) 
    for(k=0; k < N_CX_IN; ++k) 
      g_ampa_cx_in[i][j][k] = 0; 
  g_ext_in = 0; 
} 

//----------here we are changing some variables to introduce random variability--------- 

srand(6);

   for(i=0; i < Mcx; ++i)
    for(j=0; j < Mcx1; ++j){
     R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
     cx_cell[i][j].sigma =
                    cx_cell[i][j].sigma + R * 0.02;
    }
cout << "\n -------------------------------------------------"; 
 
printf("\n VARIABILITY"); 
printf("\n VARIABILITY"); 
 
//--------------open ALL files------------------------------------- 
  f11 = fopen("time_cx", "w");  
  f13 = fopen("time_in", "w"); 

  sprintf(namecx,"graf_cx%d", TTRS);
  f98 = fopen(namecx, "w");
  sprintf(namein,"graf_in%d", TTRS);
  f99 = fopen(namein, "w");

//----------Connection matrix------------------------------- 
printf("\n Begin Connect Matrix"); 

//---set connection matrix to zero--------------------------
for(i=0; i<(2*MS_CX_CX+1); i++) 
  for(j=0; j<(2*MS_CX_CX1+1); j++){ 
     k_CXCX=0; 
     C_CXCX[i][j]=0; 
  } 
 
for(i=0; i<(2*MS_CX_IN+1); i++) 
  for(j=0; j<(2*MS_CX_IN1+1); j++){ 
     k_CXIN=0; 
     C_CXIN[i][j]=0; 
  } 
 
for(i=0; i<(2*MS_IN_CX+1); i++) 
  for(j=0; j<(2*MS_IN_CX1+1); j++){ 
     k_INCX=0;
     C_INCX[i][j]=0; 
  } 

//---start updating connection matrix-------------------------
for(i1=-MS_CX_CX; i1<=MS_CX_CX; i1++) 
for(j1=-MS_CX_CX1; j1<=MS_CX_CX1; j1++){ 
  scale = sqrt((double) (i1*i1 + j1*j1)); 
  if(scale <= MS_CX_CX_MAX){ 
       C_CXCX[i1+MS_CX_CX][j1+MS_CX_CX1]=1; 
       k_CXCX=k_CXCX+1;} 
  if( (scale == 0) && (SELFINGcx == 0)  ){ 
       C_CXCX[i1+MS_CX_CX][j1+MS_CX_CX1]=0; 
       k_CXCX=k_CXCX-1;} 
  } 

for(i1=-MS_IN_CX; i1<=MS_IN_CX; i1++) 
for(j1=-MS_IN_CX1; j1<=MS_IN_CX1; j1++){ 
  scale = sqrt((double) (i1*i1 + j1*j1)); 
  if(scale <= MS_IN_CX_MAX){ 
       C_INCX[i1+MS_IN_CX][j1+MS_IN_CX1]=1; 
       k_INCX=k_INCX+1;} 
  } 

for(i1=-MS_CX_IN; i1<=MS_CX_IN; i1++) 
for(j1=-MS_CX_IN1; j1<=MS_CX_IN1; j1++){ 
  scale = sqrt((double) (i1*i1 + j1*j1)); 
  if(scale <= MS_CX_IN_MAX){ 
       C_CXIN[i1+MS_CX_IN][j1+MS_CX_IN1]=1; 
       k_CXIN=k_CXIN+1;} 
  } 

k_CXCXmax=k_CXCX;
k_CXINmax=k_CXIN;
k_INCXmax=k_INCX;

printf("\n End Connect Matrix");

//----------------CALCULATION----------------------------------------
printf("\n CALCULATION IN PROGRESS!!!: t= %lf: tmax= %lf", t,tmax);
t=0;

printf("\n start simulation"); 
while( t < tmax){
    t++;
    fun(t);

// External (DC type) stimulation to 10 pyramidal neurons (only applied once)
    for(i=0; i<10; i++){
        if(t>300) cx_cell[i][Mcx1-1].Idc=0.2;
        if(t>320) cx_cell[i][Mcx1-1].Idc=0;
           }

//------- Calculate average response--------------
aver=0;
for(i=Mcx/2-4; i<=Mcx/2+5; i++)
  for(j=Mcx1/2-4; j<=Mcx1/2+5; j++){
     aver=aver+cx_cell[i][j].x;
}


//--------scaling of all conductances-----------------------------------
if( (t>200.1) && (t<=201.5) ){
printf("\n CONDUCTANCE INITIALIZATION");
if(I_CX == 1){
  g_ampa_cx_cx[0][0][0] = g_AMPA_CX_CX / cx_cell[0][0].S_CX_DEND / k_CXCXmax;
  g_gaba_a_in_cx[0][0][0] = g_GABA_A_IN_CX / cx_cell[0][0].S_CX_DEND / k_INCXmax;
}
if(I_IN == 1){
  g_ampa_cx_in[0][0][0] = g_AMPA_CX_IN / in_cell[0][0].S_CX_DEND / k_CXINmax;
}
 printf("\n COUPLING %lf %lf %lf",g_ampa_cx_cx[0][0][0], g_gaba_a_in_cx[0][0][0], g_ampa_cx_in[0][0][0]);
}


//-----------------------------------------------------------------------
   av=0;
   if((t/(5000))*(5000) == t){
      printf("\n %lf", t/tmax);
      }    

//----Print output to files----------------

   if((t > ttime) && (((int)t/(2))*(2) == (int)t)){
     if(I_CX == 1) fprintf(f11,"%lf ",t);
     if(I_IN == 1) fprintf(f13,"%lf ",t);

     for(i = 0; i < Mcx; ++i){
         if(I_CX == 1) fprintf(f11,"%lf ", cx_cell[i][0].x);
     }
     for(i = 0; i < Min; ++i){
       if(I_IN == 1) fprintf(f13,"%lf ", in_cell[i][0].x);
     }
     fprintf(f11,"\n");
     fprintf(f13,"\n");
   }

   if((t > t3D) && (((int)t/(10))*(10) == (int)t)){

     for(i = 0; i < Mcx; ++i){
       for(j = 0; j < Mcx1; ++j){
         if(I_CX == 1) fprintf(f98,"%lf ", cx_cell[i][j].x);
       }
       fprintf(f98,"\n");
     } 

     for(i = 0; i < Min; ++i){
       for(j = 0; j < Min1; ++j){
         if(I_IN == 1) fprintf(f99,"%lf ", in_cell[i][j].x);
       }
       fprintf(f99,"\n");
     }
   }
 }

//--------------------END CALCULATION------------------------------- 
 
//-----------------close ALL files----------------------------------- 
 
  fclose(f11); 
  fclose(f13); 

  fclose(f98);
  fclose(f99);

  printf("\n");  
} 
 
 
//+++++++++++ Function to calculate the right sides for ALL ODE +++++++++++++++ 
void fun(double t){ 
int i, j, k, kmax, i1, j1, ii, jj; 
double IsynCXCX[Mcx][Mcx1], IsynCXIN[Min][Min1], IsynINCX[Mcx][Mcx1]; 
 
for(i = 0; i < Mcx; ++i) 
for(j = 0; j < Mcx1; ++j){ 
 
//-----calculating reciprocal AMPA-mediated synaptic currents between CX cells--------------
if(I_CX == 1){ 
  for(i1 = -MS_CX_CX, k = 0; i1 <= MS_CX_CX; ++i1) 
  for(j1 = -MS_CX_CX1; j1 <= MS_CX_CX1; ++j1){ 
     ii=i+i1, jj=j+j1;
     if((ii>=0) && (ii<Mcx) && (jj>=0) && (jj<Mcx1)){ 
       if(C_CXCX[i1+MS_CX_CX][j1+MS_CX_CX1] > 0){ 
          a_cx_cx[i][j][k].calc(g_ampa_cx_cx[0][0][0], cx_cell[i][j].x, cx_cell[ii][jj].spike); 

        ++k; } 
     }} 
  kmax = k; 
  IsynCXCX[i][j] = 0; 
  for(k = 0; k < kmax; ++k){ 
    if(kmax == 0) { }  //NO synapses 
    else {  
       IsynCXCX[i][j] = IsynCXCX[i][j] + a_cx_cx[i][j][k].I;
   } 
  } 
 } 

//------calculating GABA-A type current from IN to CX cells-------------------------------- 
if((I_CX == 1) && (I_IN == 1)){ 
  for(i1 = -MS_IN_CX, k = 0; i1 <= MS_IN_CX; ++i1) 
  for(j1 = -MS_IN_CX1; j1 <= MS_IN_CX1; ++j1){ 
     ii=i*Min/Mcx+i1, jj=j*Min1/Mcx1+j1;
     if((ii>=0) && (ii<Min) && (jj>=0) && (jj<Min1)){ 
              if(C_INCX[i1+MS_IN_CX][j1+MS_IN_CX1] > 0){ 
         ga_in_cx[i][j][k].calc(g_gaba_a_in_cx[0][0][0], cx_cell[i][j].x, in_cell[ii][jj].spike); 

        ++k; } 
     } 
   }
  kmax = k; 
  IsynINCX[i][j] = 0; 
  for(k = 0; k < kmax; ++k){ 
    if(kmax == 0) { }  //NO synapses 
    else {  
       IsynINCX[i][j] = IsynINCX[i][j] + ga_in_cx[i][j][k].I;
  } 
  } 
} 
   
 }
 
for(i = 0; i < Min; ++i) 
for(j = 0; j < Min1; ++j){ 

//------calculating AMPA type current from CX to IN cells-------------------------------- 
if((I_CX == 1) && (I_IN == 1)){ 
  for(i1 = -MS_CX_IN, k = 0; i1 <= MS_CX_IN; ++i1) 
  for(j1 = -MS_CX_IN1; j1 <= MS_CX_IN1; ++j1){
    ii=i*Mcx/Min+i1, jj=j*Mcx1/Min1+j1;
    if((ii>=0) && (ii<Mcx) && (jj>=0) && (jj<Mcx1)){ 
            if(C_CXIN[i1+MS_CX_IN][j1+MS_CX_IN1] > 0){  
        a_cx_in[i][j][k].calc(g_ampa_cx_in[0][0][0], in_cell[i][j].x, cx_cell[ii][jj].spike); 

        ++k; } 
    } 
    }
 
  kmax = k; 
  IsynCXIN[i][j] = 0; 
  for(k = 0; k < kmax; ++k){ 
    if(kmax == 0) { }  //NO synapses 
    else {  
       IsynCXIN[i][j] = IsynCXIN[i][j] + a_cx_in[i][j][k].I;

   } 
  } 
 }   
 }

//------calculating the state of all neurons at the next iteration based on synaptic input-----

for(i=0; i < Mcx; ++i) 
  for(j=0; j < Mcx1; ++j){ 
  if(I_CX == 1) cx_cell[i][j].calc(IsynCXCX[i][j]+IsynINCX[i][j]); 
} 
 
for(i=0; i < Min; ++i) 
  for(j=0; j < Min1; ++j){ 
     if(I_IN == 1) in_cell[i][j].calc(IsynCXIN[i][j]); 
  } 

}
 
