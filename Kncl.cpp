/*******************************************************************

Conductance based neocortical network model
Modified version of code by MB (2006)

********************************************************************/

// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "params.h"

#define NUM_THREADS 40

// Constants
#define PI 3.141592654

#define input_start_time 100000
#define input_end_time 100400

#define Check_time 1000

#define no_clusters 10
#define Mcx_cluster_size 50
#define Min_cluster_size 10

// Network Geometry

#define I_CX    1     //  0 - No layer, 1 - Add Layer
#define I_IN    1     //  0 - No layer, 1 - Add Layer
#define I_GB    0     //  0 - No GABAb from IN to CX, 1 - Yes

#define Mcx      (Mcx_cluster_size*no_clusters)   //  Number of CX cells in X direction
#define Min      (Min_cluster_size*no_clusters)    //  Number of IN cells in X direction
#define Mcx1     1     //  Number of CX cells in Y direction
#define Min1     1     //  Number of IN cells in Y direction
#define Max     Mcx
#define Max1    Mcx1


#define Num_Types 2
#define Num_Cells ((Min*Min1)+(Mcx*Mcx1))
int cell_sizes[2][2] = {{Mcx,Mcx1},
                        {Min,Min1}}; 

double Cluster_conn_matrix[no_clusters][no_clusters]; //cluster edge i->o == Cluster_conn_matrix[i][o]

// Boundary conditions
#define BOUND      0   // 0 - flow; 1 - periodic; 9 - no boundary elements
#define BOUNDcx    0   // 0 - flow; 1 - periodic; 9 - no boundary elements
#define SELFINGcx  0   // 0 - CX without self excitation; 1 - with ...

// Connections
#define MS_CX_CX 5
#define MS_CX_CX1 0
#define MS_CX_CX_MAX  ((MS_CX_CX > MS_CX_CX1) ? MS_CX_CX : MS_CX_CX1)
#define N_CX_CX  1

#define MS_CX_CX_NMDA 5
#define MS_CX_CX_NMDA1 0
#define MS_CX_CX_NMDA_MAX  ((MS_CX_CX_NMDA > MS_CX_CX_NMDA1) ? MS_CX_CX_NMDA : MS_CX_CX_NMDA1)
#define N_CX_CX_NMDA  1

#define MS_CX_IN_NMDA 1
#define MS_CX_IN_NMDA1 0
#define MS_CX_IN_NMDA_MAX  ((MS_CX_IN_NMDA > MS_CX_IN_NMDA1) ? MS_CX_IN_NMDA : MS_CX_IN_NMDA1)
#define N_CX_IN_NMDA 1

#define MS_CX_IN 1
#define MS_CX_IN1 0
#define MS_CX_IN_MAX  ((MS_CX_IN > MS_CX_IN1) ? MS_CX_IN : MS_CX_IN1)
#define N_CX_IN  1

#define MS_IN_CX 5
#define MS_IN_CX1 0
#define MS_IN_CX_MAX  ((MS_IN_CX > MS_IN_CX1) ? MS_IN_CX : MS_IN_CX1)
#define N_IN_CX  ((2*MS_IN_CX+1)+4)*Min/Mcx

#define MS_IN_IN 2
#define MS_IN_IN1 0
#define MS_IN_IN_MAX  ((MS_IN_IN > MS_IN_IN1) ? MS_IN_IN : MS_IN_IN1)
#define N_IN_IN  (2*MS_IN_IN+1)

// Number of ODE for each cell
#define N_RE 7
#define N_TC 12
#define N_GB 2

#define N_DEND   17
#define N_SOMA   9
#define N_CX     (N_DEND + N_SOMA)
#define N_IN     N_CX

#define N_EQ2  (N_CX*I_CX + N_IN_CX*N_GB*I_CX*I_IN*I_GB)*Mcx*Mcx1
#define N_EQ3  (N_IN*I_IN)*Min*Min1
#define N_EQ   (N_EQ2 + N_EQ3)    //  Total number of ODE

// Reversal potentials
#define EE_Na 50
#define EE_l -61
#define EE_K -95
#define EE_h -50
#define E_LK = -58

int **index_loc_array;
int seizure=0;
double syn_factor=1.5;

// Ratio of the CX-CX, CX-IN connection for between cluster
double Long_range_cxcx_strength=1;
double Long_range_cxin_strength=0.0;


// Declare all the current classes
#include "current_classes.h"

//+++++++++++++++++++ MAIN PROGRAM +++++++++++++++++++++++++++++++++++++++++++

//----------external functuons----------------------------------------------
void index_array(int***);
void rk(unsigned, void (unsigned, double, double*, double*), double, double,
        double*, double*, double*, double*);
void calcit(int index,unsigned neq, double x, double *y_ini, double *f_ini);
void calc_syn(int index,unsigned neq, double x, double *y_ini, double *f_ini);
void fun(unsigned, double, double*, double*);
double gaussrand();

//----------external variables ---------------------------------------------
double g_ampa_cx_cx[Mcx][Mcx1], g_nmda_cx_cx, g_nmda_cx_in;
double  g_ampa_cx_cx_MINICE, g_ampa_cx_in_MINICE;
double g_ampa_cx_in[Min][Min1], g_gaba_a_in_cx[Mcx][Mcx1];
double g_gaba_b_in_cx, g_gaba_a_in_in;
double g_ext_cx1, g_ext_in1;
FILE *finput, *f_cx, *f_cx_syn, *f_cx_syn_inh, *f_in, *f_py_rates, *g_py_py,*f_in_rates, *g_py_in, *MyFile, *g_in_py, *g_GHVA,*g_GNAP,*g_GKCA, *SynCurrent, *myconnect,*f1000,*f_gaba, *f_par, *log_file, *fsp;

int no_cx[Mcx][Mcx1][N_CX]; 
int no_in[Min][Min1][N_IN], no_gcx[Mcx][Mcx1][N_IN_CX][N_GB];

int C_CXCX[Mcx][Mcx1][Mcx][Mcx1];
int C_CXIN[Mcx][Mcx1][Min][Min1];
int C_INCX[Min][Min1][Mcx][Mcx1];
int C_ININ[Min][Min1][Min][Min1];

double Cluster_strength_CXCX[Mcx][Mcx1][Mcx][Mcx1];
double Cluster_strength_CXIN[Mcx][Mcx1][Min][Min1];

int k_CXCX[Mcx][Mcx1];
int k_CXIN[Min][Min1];
int k_INCX[Mcx][Mcx1];
int k_ININ[Min][Min1];

int k_CXCXmax=0;
int k_CXINmax=0;
int k_INCXmax=0;
int k_ININmax=0;
int firstflag=0;

double prev[Mcx][Mcx1];
double spikecount[Mcx][Mcx1];

double average_rate_py=0;

double prev_in[Min][Min1];
double spikecount_in[Min][Min1];

double average_rate_in=0;
double E[Mcx][Mcx1];

double GPY = 0;
double GIN = 0;
double GPYPY = 0;
double GPYIN = 0;

int excounter = 0;

int randseed = 0;
//----------external classes (beginning of initialization)------------------
AMPA_D2         a_cx_cx[Mcx][Mcx1];
NMDA_D1         nmda_cx_cx[Mcx][Mcx1];
NMDA_D1         nmda_cx_in[Min][Min1];
AMPA_D2         a_cx_in[Min][Min1];
Gaba_A_D2       ga_in_in[Min][Min1];
Gaba_A_D2       ga_in_cx[Mcx][Mcx1];

GB           *gb_in_cx[Mcx][Mcx1];

CX           cx_cell[Mcx][Mcx1];
CX           in_cell[Min][Min1];

Extern_ampa1  a_ext33[Mcx][Mcx1], a_ext44[Min][Min1];
double my_E_GABA[Mcx][Mcx1];
double my_E_GABA1[Mcx][Mcx1];

double p_gNals, p_gKls, p_gNald, p_gKld, p_gld, p_cldrive, p_cltau, p_extinp;
double networks_ampa[1][no_clusters+2]; //double is needed really for the first column only
//todo merge two funs below
double network_get_node_ampa(int n) { //return multiplicative ampa change if the array contains node n
  if (networks_ampa[0][0]==0) return 1;
  int i=0; while (networks_ampa[0][++i]!=0) {
    if (networks_ampa[0][i]==n) {/*printf("C:%dL:%lf ",n,networks_gkld[0][0]);*/ return networks_ampa[0][0];}
  }
  return 1;
}

double networks_gkld[1][no_clusters+2]; //double is needed really for the first column only
double network_gkld_get_node_leak(int n, double p_gKld) { //return leak (in [0][0] if the array contains node n
  if (networks_gkld[0][0]==0) return p_gKld;
  int i=0; while (networks_gkld[0][++i]!=0) {
    if (networks_gkld[0][i]==n) {/*printf("C:%dL:%lf ",n,networks_gkld[0][0]);*/ return networks_gkld[0][0];}
  }
  return p_gKld;
}

// ************************************************************************
int main(int argc,char **argv)
{
  //---------allocate place for ALL variables and functions-----------
  double y_ini[N_EQ], f_ini[N_EQ];

  //---------allocate place for TWO temporal arrays using by RK.C solver
  double y1[N_EQ], y2[N_EQ];

  //---------general parameters----------------------------------------------
  double t = 0, tmax, tmax1, ttime, TAU;
  double h = 0.04, g_Extern_ampa1, R;
  int i, j, i1, j1, k, ih, ii = 0, i_inp;
  double scale;
  double g_AMPA_CX_CX, g_NMDA_CX_CX, g_NMDA_CX_IN, g_AMPA_CX_CX_MINICE, g_AMPA_CX_IN_MINICE;
  double g_AMPA_CX_IN, g_GABA_A_IN_CX, g_GABA_B_IN_CX;
  double g_GABA_A_IN_IN, g_Extern_ampa;

  double *tmpa;

  int N_deaff = 10;
  int myvar = 0;
  double myvard=0.0;
  ttime=0;
  int no_input=0;

  index_array(&index_loc_array);

  log_file = fopen("logfile", "w");
  //----------arrays initialization----------------------------------------------
  for(i=N_EQ-1; i>=0; --i){
    y_ini[i] = 0, f_ini[i] = 0;
    y1[i] = 0, y2[i] = 0; }

  //-------parameter initialization (from file)----------------------------------
  if (argc <= 1) {
    puts("Command parameters");
    puts("-----------------------");
    puts("Input File"); }

  f_par = fopen("par", "r");
  no_input=fscanf(f_par, "%lf %lf %lf %lf %lf %lf %lf %lf",
                  &p_gKld, &p_gNald, &p_gKls, &p_gNals, &p_gld, &p_cldrive, &p_cltau, &p_extinp);
  printf("\n %d Var Parameters: p_gKls %lf, p_gNals %lf, p_gKld %lf, p_gNald %lf, p_gld %lf, p_cldrive %lf, p_cltau %lf, p_extinp %lf", 
         no_input, p_gKls, p_gNals, p_gKld, p_gNald, p_gld, p_cldrive, p_cltau, p_extinp);
  fclose(f_par);

  if (!(finput=fopen(argv[1],"r"))) {
    printf("%s doesn't exist\n",argv[1]);
    exit(0); }

  no_input=fscanf(finput, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
                  &tmax, &g_AMPA_CX_CX, &g_NMDA_CX_CX,&g_AMPA_CX_IN, &g_NMDA_CX_IN, &g_GABA_A_IN_CX,
                  &g_GABA_B_IN_CX, &g_GABA_A_IN_IN,&g_Extern_ampa, &g_Extern_ampa1, &i_inp);
  printf("\n %d Parameters: %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d",
         no_input, tmax, g_AMPA_CX_CX, g_NMDA_CX_CX,g_AMPA_CX_IN, g_NMDA_CX_IN, g_GABA_A_IN_CX,
         g_GABA_B_IN_CX, g_GABA_A_IN_IN,g_Extern_ampa, g_Extern_ampa1, i_inp);
  fclose(finput);

  string clusters; //global network connectivity, local (within cluster) will be generated below
  string network_gkld_file; //set of global nodes with different gKld value
  string network_ampa_file; //set of global nodes with different ampa strength
  double cluster_scale=1.0;
  double cluster_cutoff=0.01;
  add_int_param(randseed);
  add_string_param(clusters);
  add_double_param(cluster_scale);
  add_double_param(p_gKld);
  add_string_param(network_gkld_file);
  add_string_param(network_ampa_file);
  //assert(load_parameters(argv[1]));
  assert(cmdline_parameters(argc,argv));
  print_parameters();

  // double **Cluster_conn_matrix;
  // Cluster_conn_matrix = new double *[no_clusters];
  // for (int cli=0; cli<no_clusters; cli++)
  //   Cluster_conn_matrix[cli]=new double[no_clusters];


  if (!clusters.empty()) {
    FILE *fconn;
    fconn = fopen(clusters.c_str(),"r");
    double conn; int rlines=0;
    for (int i=0; i<no_clusters; i++)
      for (int j=0;j<no_clusters; j++){
        if (fscanf(fconn, "%lf", &conn)!=1)
         printf("Network input error\n"),exit(1);
	 
        conn=conn*cluster_scale;
        if (conn>cluster_cutoff)
         Cluster_conn_matrix[i][j]=conn;
      }//for j
    
  } else {

  FILE *fconn, *ftmp;
  fconn = fopen("conn","r");
  ftmp = fopen("rconn","w");
  int cx, cy; double conn; int rlines=0;
  while(1){
    if (fscanf(fconn, "%d %d %lf", &cx, &cy, &conn)!=3)
      break;
    else{
      rlines=rlines+1;
      Cluster_conn_matrix[cx-1][cy-1]=conn;
      fprintf(ftmp, "Connectivity %lf \n", Cluster_conn_matrix[cx-1][cy-1]);
    }
  }

  printf ("Read %d lines from conn file ", rlines);
  if (rlines==0){
    printf("ERROR!!!");
      exit(0);
  }
  fclose(ftmp);

  } //cluster param empty

  // Expected file format: <g_kld value> <number of N rsns nodes> <node 1 index> ... <node N index> 
  // Node indexing starts from 1 to stay in tune with Matlab scripts
  if (!network_gkld_file.empty()) {
    FILE *fconn;
    fconn = fopen(network_gkld_file.c_str(),"r");

    double leak;        
    int nnodes,node;
    if (fscanf(fconn, "%lf %d", &leak, &nnodes)!=2) printf("RSNs leak input error\n"),exit(1);
    networks_gkld[0][0]=leak;

    for (int i=0;i<nnodes;i++) {
      if (fscanf(fconn, "%d", &node)!=1) printf("RSNs leak input error\n"),exit(1);
      networks_gkld[0][i+1]=node;
    }

    printf("\nRSNs leakage: %lf, Nodes: ",networks_gkld[0][0]);
    i=1; while (networks_gkld[0][i]!=0) printf("%d ",(int)networks_gkld[0][i++]); printf("\n");

  } //network_gkld_file

  // Expected file format: <multiplicative ampa value> <number of N rsns nodes> <node 1 index> ... <node N index> 
  // Node indexing starts from 1 to stay in tune with Matlab scripts
  if (!network_ampa_file.empty()) {
    FILE *fconn;
    fconn = fopen(network_ampa_file.c_str(),"r");

    double ampa;        
    int nnodes,node;
    if (fscanf(fconn, "%lf %d", &ampa, &nnodes)!=2) printf("RSNs ampa input error\n"),exit(1);
    networks_ampa[0][0]=ampa;

    for (int i=0;i<nnodes;i++) {
      if (fscanf(fconn, "%d", &node)!=1) printf("RSNs ampa input error\n"),exit(1);
      networks_ampa[0][i+1]=node;
    }

    printf("\nRSNs ampas: %lf, Nodes: ",networks_ampa[0][0]);
    i=1; while (networks_ampa[0][i]!=0) printf("%d ",(int)networks_ampa[0][i++]); printf("\n");

  } //network_ampa_file

  for (i=0; i<Mcx; i++)
    for (j=0; j<Mcx1; j++)
      for (k=0; k<N_CX; k++){
        no_cx[i][j][k] =   k + (j + i*Mcx1) * N_CX;
        //printf("\ni:%d j:%d k:%d no_cx:%d",i,j,k,no_cx[i][j][k]);
      }

  for (i=0; i<Min; i++)
    for (j=0; j<Min1; j++)
      for (k=0; k<N_IN; k++){
        no_in[i][j][k] =   Mcx*Mcx1*N_CX*I_CX + k + (j + i*Min1) * N_IN;
        //printf("\ni:%d j:%d k:%d no_in:%d",i,j,k,no_in[i][j][k]);
      }

  //---variable initialization (additional for standart constructor)----------

  for(i=Mcx-1; i>=0; --i)
    if(I_CX == 1){
      for(j=Mcx1-1; j>=0; --j){
        spikecount[i][j]=0;
        prev[i][j]=-10000;
        cx_cell[i][j].init(y_ini+no_cx[i][j][0]);

        tmpa=y_ini+no_cx[i][0][0];

        a_cx_cx[i][j].iii(1+i*(j));

        if(I_IN == 1) {
          ga_in_cx[i][j].iii(1+27+i*j);
        }
	a_ext33[i][j].iii(randseed+i*1000);
      }
    }

 
  for(j=Min1-1; j>=0; --j)
    for(i=Min-1; i>=0; --i)
    {
      if(I_IN == 1){
        spikecount_in[i][j]=0;
        prev_in[i][j]=0;

        in_cell[i][j].init(y_ini+no_in[i][j][0]);
        double tmp=in_cell[i][0].v_SOMA;
        double *tmp1=y_ini+no_in[i][0][0];
	
        if(I_CX == 1) {
          a_cx_in[i][j].iii(1+17+i*j);
        }
	a_ext44[i][j].iii(randseed+i*10200);
      }
    }

  if(I_CX == 1){
    for(j=Mcx1-1; j >=0; --j)
      for(i=Mcx-1; i >=0; --i)
      {

        g_ampa_cx_cx[i][j] = 0;
        g_gaba_a_in_cx[i][j] = 0;
      }
    g_nmda_cx_cx = 0;
  }

  if(I_IN == 1){
    for(j=Min1-1; j >=0; --j)
      for(i=Min-1; i >=0; --i)


        g_ampa_cx_in[i][j] = 0;

    g_gaba_a_in_in = 0;
  }

  //----------here we are changing some variables----------------------

  srand(3000+randseed);

  for(j=Min1-1; j >=0; --j)
    for(i=Min-1; i >=0; --i){
      in_cell[i][j].rho = 50;
      in_cell[i][j].CX_DEND::G_Nap = 0.0;
      in_cell[i][j].CX_SOMA::G_Nap = 0.0;

      in_cell[i][j].I_Stim1 = 0;
      in_cell[i][j].G_h = 0.00;

      in_cell[i][j].G_kl = 0.045; // **** K leak dendrite
      in_cell[i][j].G_Nal = 0.024; // **** Na leak dendrite
      in_cell[i][j].g_kl = 0.045;  // **** K leak soma
      in_cell[i][j].g_Nal = 0.024; // **** Na leak soma

      in_cell[i][j].CX_DEND::Imax = 20;
      in_cell[i][j].CX_SOMA::Imax = 20;

      in_cell[i][j].CX_DEND::ICl::Cl_inf = 5.0;
      in_cell[i][j].CX_DEND::ICl::Cl_inf1 = 5.0; 
      in_cell[i][j].CX_DEND::ICl::Tau_inf = 0.5;
      in_cell[i][j].CX_DEND::ICl::Cl_drive = 50;
      in_cell[i][j].CX_DEND::ICl::Cl_Tauinf = 20000;

      in_cell[i][j].CX_SOMA::Na_sc = 500; 
      in_cell[i][j].CX_SOMA::K_sc = 500;

      in_cell[i][j].CX_DEND::G_l = 0.010;
      in_cell[i][j].CX_SOMA::g_l = 0.0;
      in_cell[i][j].CX_SOMA::G_Na = 3200;
      in_cell[i][j].CX_SOMA::G_Kv = 200.0;
      in_cell[i][j].CX_DEND::G_Na = 1.0;

      in_cell[i][j].CX_SOMA::G_KNa = 0.0;
      in_cell[i][j].CX_DEND::G_Nap = 0.0;
      in_cell[i][j].CX_SOMA::G_Nap = 0.0;
      in_cell[i][j].ICa::D = 0.0;

      in_cell[i][j].G_h = 0.1;
      in_cell[i][j].G_KCa = 0.0;
      in_cell[i][j].G_HVA = 0.0;
      in_cell[i][j].CX_DEND::G_Km = 0.0; 

    }

  
  for(j=Mcx1-1; j >=0; --j)
    for(i=Mcx-1; i >=0; --i){
  
      double r = (double)rand() / RAND_MAX;

      g_AMPA_CX_CX_MINICE = 0;
      g_AMPA_CX_IN_MINICE = 0;
      int Ni=i/Mcx_cluster_size; // i/50
      Ni++; //matlab indexing of ROIs
      double p_gKld2 = network_gkld_get_node_leak(Ni, p_gKld); // p_gKld if does not exists

      cx_cell[i][j].G_kl = p_gKld2;
      cx_cell[i][j].G_Nal = p_gNald;
      cx_cell[i][j].g_kl = p_gKls;
      cx_cell[i][j].g_Nal = p_gNals;

      cx_cell[i][j].CX_DEND::G_l = p_gld;
      cx_cell[i][j].CX_SOMA::g_l = 0.0;

      cx_cell[i][j].CX_DEND::Imax =20;
      cx_cell[i][j].CX_SOMA::Imax = 20;
 
      cx_cell[i][j].G_h = 0.1;
      cx_cell[i][j].G_KCa = 2.5;
      cx_cell[i][j].G_HVA = 0.013;
		
      cx_cell[i][j].CX_SOMA::G_KNa = 1.5;
      cx_cell[i][j].CX_SOMA::G_Kv = 200.0;

      cx_cell[i][j].CX_DEND::G_Na = 1.1;
      cx_cell[i][j].CX_SOMA::G_Na = 3500;

      cx_cell[i][j].CX_DEND::G_Nap = 3.5;
      cx_cell[i][j].CX_SOMA::G_Nap = 3.5;
      
      cx_cell[i][j].CX_DEND::G_Km = 0.01;
      cx_cell[i][j].ICa::Taur = 300;
      cx_cell[i][j].ICa::D = 0.85;
      cx_cell[i][j].rho = 165;
      cx_cell[i][j].I_Stim1 = 0;

      cx_cell[i][j].CX_SOMA::Na_sc = 500;
      cx_cell[i][j].CX_SOMA::K_sc = 500;

      cx_cell[i][j].CX_DEND::ICl::Cl_inf = 5.0;
      cx_cell[i][j].CX_DEND::ICl::Cl_inf1 = 5.0; 
      cx_cell[i][j].CX_DEND::ICl::Tau_inf = 0.08;
      cx_cell[i][j].CX_DEND::ICl::Cl_drive = p_cldrive;
      cx_cell[i][j].CX_DEND::ICl::Cl_Tauinf = p_cltau;

    }

  // --------------------------------

  f_cx = fopen("time_cx_new", "w");
  f_cx_syn = fopen("time_cx_syn", "w");
  f_cx_syn_inh = fopen("time_cx_syn_inh", "w");
  f_in = fopen("time_in_new", "w");
  f_py_rates = fopen("py_rates", "w");
  f_in_rates = fopen("in_rates", "w");
  g_py_py = fopen("E_py_py","w");
  g_py_in = fopen("py_in","w");
  g_in_py = fopen("in_py","w");
  myconnect = fopen("connectivity","w");
  f1000 = fopen("time_K", "w");
  fsp = fopen("spikes", "w");
  ///----------Connection matrix-------------------------------

  printf("\n Begin to set up connection matrix...");

  for(j=Mcx1-1; j>=0; --j)
    for(i=Mcx-1; i>=0; --i){
      k_CXCX[i][j]=0;
      for(j1=Mcx1-1; j1>=0; --j1)
        for(i1=Mcx-1; i1>=0; --i1){
          C_CXCX[i][j][i1][j1]=0;
          Cluster_strength_CXCX[i][j][i1][j1]=0;
        }
      
    }

  for(i=Min-1; i>=0; --i)
    for(j=Min1-1; j>=0; --j){
      k_ININ[i][j]=0;

      for(j1=Min1-1; j1>=0; --j1)
        for(i1=Min-1; i1>=0; --i1) C_ININ[i][j][i1][j1]=0;
    }

  for(i1=Min-1; i1>=0; --i1)
    for(j1=Min1-1; j1>=0; --j1){
      k_CXIN[i1][j1]=0;

      for(j=Mcx1-1; j>=0; --j)
        for(i=Mcx-1; i>=0; --i) {
          C_CXIN[i][j][i1][j1]=0;
          Cluster_strength_CXIN[i][j][i1][j1]=0;          
        }
    }

  for(i1=Mcx-1; i1>=0; --i1)
    for(j1=Mcx1-1; j1>=0; --j1){
      k_INCX[i1][j1]=0;
      for(j=Min1-1; j>=0; --j)
        for(i=Min-1; i>=0; --i)
          C_INCX[i][j][i1][j1]=0;
    }


  //   ****************** Code for cluster ****************
 
  //to ensure correct radius computation we need to calculate with double
  double Min_=Min, Min1_=Min1, Mcx_=Mcx, Mcx1_=Mcx1;
  //  Clusters in 1D
  for(int cli=0; cli<no_clusters; cli++) {
    printf("Cluster: %d ", cli);
  
    for(i1=cli*Mcx_cluster_size; i1<(cli+1)*Mcx_cluster_size; i1++){
      for(j1=0; j1<Mcx1; j1++){
        for(i=cli*Mcx_cluster_size; i<(cli+1)*Mcx_cluster_size; i++){
          for(j=0; j<Mcx1; j++){
            scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
            R = 0;
            if(scale <= MS_CX_CX_MAX){
              if(R <= 0.2){
                C_CXCX[i1][j1][i][j]=1;
                k_CXCX[i][j]=k_CXCX[i][j]+1;}
            }
	      
            if( (scale == 0) && (SELFINGcx == 0)  ){
              C_CXCX[i1][j1][i][j]=0;
              k_CXCX[i][j]=k_CXCX[i][j]-1;}

            C_CXCX[i1][j1][i][j] && fprintf(myconnect,"%d %d %d %d %d %d %d \n",1, cli, i1, j1, i, j, C_CXCX[i1][j1][i][j]);	//CX->CX 1
          }
        }
      }
    }
  
    for(i1=cli*Min_cluster_size; i1<(cli+1)*Min_cluster_size; i1++) //(i1=0; i1<Min; i1++)
      for(j1=0; j1<Min1; j1++)
        for(i=cli*Min_cluster_size; i<(cli+1)*Min_cluster_size; i++) //(i=0; i<Min; i++)
          for(j=0; j<Min1; j++){
            scale = sqrt((double) ((i1-i)*(i1-i) + (j1-j)*(j1-j)));
            if(scale <= MS_IN_IN_MAX){
              C_ININ[i1][j1][i][j]=1;
              k_ININ[i][j]=k_ININ[i][j]+1;}
            if( scale == 0 ){
              C_ININ[i1][j1][i][j]=0;
              k_ININ[i][j]=k_ININ[i][j]-1;}
            C_ININ[i1][j1][i][j] && fprintf(myconnect,"%d %d %d %d %d %d %d \n",2, cli, i1, j1, i, j, C_ININ[i1][j1][i][j]);	//IN->IN 2
          }

    for(i1=cli*Min_cluster_size; i1<(cli+1)*Min_cluster_size; i1++){ //(i1=0; i1<Min; i1++){
      for(j1=0; j1<Min1; j1++){
        for(i=cli*Mcx_cluster_size; i<(cli+1)*Mcx_cluster_size; i++){ //for(i=0; i<Mcx; i++){
          for(j=0; j<Mcx1; j++){
            scale = sqrt((double) ((i1*Mcx_/Min_-i)*(i1*Mcx_/Min_-i)
                                   + (j1*Mcx1_/Min1_-j)*(j1*Mcx1_/Min1_-j)));
            R = 0;
            if(scale <= MS_IN_CX_MAX){
              if(R <= 0.2){
                C_INCX[i1][j1][i][j]=1;
                k_INCX[i][j]=k_INCX[i][j]+1;}
              C_INCX[i1][j1][i][j] && fprintf(myconnect,"%d %d %d %d %d %d %d \n",3, cli, i1, j1, i, j, C_INCX[i1][j1][i][j]);	//IN->CX 3
            }
          }
        }
      }
    }
 
    for(i1=cli*Mcx_cluster_size; i1<(cli+1)*Mcx_cluster_size; i1++){ //for(i1=0; i1<Mcx; i1++){
      for(j1=0; j1<Mcx1; j1++){
        for(i=cli*Min_cluster_size; i<(cli+1)*Min_cluster_size; i++) {//for(i=0; i<Min; i++){
          for(j=0; j<Min1; j++){
            scale = sqrt((double) ((i1*Min_/Mcx_-i)*(i1*Min_/Mcx_-i)
                                   + (j1*Min1_/Mcx1_-j)*(j1*Min1_/Mcx1_-j)));
            R = 0;
            if(scale <= MS_CX_IN_MAX){
              if(R <= 0.25){
                C_CXIN[i1][j1][i][j]=1;
                k_CXIN[i][j]=k_CXIN[i][j]+1;}
            }
              C_CXIN[i1][j1][i][j] && fprintf(myconnect,"%d %d %d %d %d %d %d \n",4, cli, i1, j1, i, j, C_CXIN[i1][j1][i][j]);	//CX->IN 4
          }
        }
      }
    }

  }


  k_CXCXmax=k_CXCX[0][0];
  k_ININmax=k_ININ[0][0];
  k_CXINmax=k_CXIN[0][0];
  k_INCXmax=k_INCX[0][0];

  for(i=0; i<Mcx; i++)
    for(j=0; j<Mcx1; j++){
      if(k_CXCX[i][j] > k_CXCXmax) k_CXCXmax=k_CXCX[i][j];
      if(k_INCX[i][j] > k_INCXmax) k_INCXmax=k_INCX[i][j];
    }
  for(i=0; i<Min; i++)
    for(j=0; j<Min1; j++) {
      if(k_ININ[i][j] > k_ININmax) k_ININmax=k_ININ[i][j];
      if(k_CXIN[i][j] > k_CXINmax) k_CXINmax=k_CXIN[i][j];
    }


  //  *************  Cluster ************* 

  // Connect between cluster --- Only in 1D

  for(int cli=0; cli<no_clusters; cli++) {
    for(i1=cli*Mcx_cluster_size; i1<(cli+1)*Mcx_cluster_size; i1++){
      for(j1=0; j1<Mcx1; j1++){
        for(int clio=0; clio<no_clusters; clio++) {
          if (!(clio==cli)) {
            for(i=clio*Mcx_cluster_size; i<(clio+1)*Mcx_cluster_size; i++){
              if (Cluster_conn_matrix[cli][clio]>0.0){
              for(j=0; j<Mcx1; j++){
                scale = sqrt((double) ((i1-(i-(clio-cli)*Mcx_cluster_size))*(i1-(i-(clio-cli)*Mcx_cluster_size)) + (j1-j)*(j1-j)));
                R = rand()/(RAND_MAX + 1.0);
                if(scale <= Mcx_cluster_size){
                  if(R <= 0.15){
                    C_CXCX[i1][j1][i][j]=2; //2 is for long range conns, they don't have NMDA
                    Cluster_strength_CXCX[i1][j1][i][j]=Cluster_conn_matrix[cli][clio];
                    k_CXCX[i][j]=k_CXCX[i][j]+1; 
                    C_CXCX[i1][j1][i][j] && fprintf(myconnect,"%d %d %d %d %d %d %lf \n",11, cli, i1, j1, i, j, Cluster_strength_CXCX[i1][j1][i][j]*Long_range_cxcx_strength); //long CX->CX 11
                   } 
                 }
               }
               }
           }
         }
       }
     }
   }
 }


  // Add connection between CX->IN
  for(int cli=0; cli<no_clusters; cli++) {
    for(i1=cli*Mcx_cluster_size; i1<(cli+1)*Mcx_cluster_size; i1++){
      for(j1=0; j1<Mcx1; j1++){
        for(int clio=0; clio<no_clusters; clio++) {
          if (!(clio==cli)) {
            for(i=clio*Min_cluster_size; i<(clio+1)*Min_cluster_size; i++){
              if (Cluster_conn_matrix[cli][clio]>0.0){
                for(j=0; j<Min1; j++){
                  scale = sqrt((double) ((i1*Min_/Mcx_-(i-(clio-cli)*Min_cluster_size))*(i1*Min_/Mcx_-(i-(clio-cli)*Min_cluster_size)) + (j1-j)*(j1-j))); 

                  R = rand()/(RAND_MAX + 1.0);
                  if(scale <= MS_CX_IN_MAX){
                    if(R <= 0.25){
                      C_CXIN[i1][j1][i][j]=2;
                      Cluster_strength_CXIN[i1][j1][i][j]=Cluster_conn_matrix[cli][clio];
                      k_CXIN[i][j]=k_CXIN[i][j]+1;
                      C_CXIN[i1][j1][i][j] && fprintf(myconnect,"%d %d %d %d %d %d %lf \n",12, cli, i1, j1, i, j, Cluster_strength_CXIN[i1][j1][i][j]*Long_range_cxin_strength); //long CX->IN 12
                    } 
                  }
                }
              }
            }
         }
       }
     }
   }
 }

 fclose(myconnect);

 // Figure out maximum connections after connecting clusters
 for(i=0; i<Mcx; i++)
   for(j=0; j<Mcx1; j++){
     if(k_CXCX[i][j] > k_CXCXmax) k_CXCXmax=k_CXCX[i][j];
   }


/// print connections
for(j=Mcx1-1; j>=0; --j)
  for(i=Mcx-1; i>=0; --i){
    // k_CXCX[i][j]=0;
    for(j1=Mcx1-1; j1>=0; --j1)
      for(i1=Mcx-1; i1>=0; --i1)
        0 && printf("%d %d %d \n", i, j, C_CXCX[i][j][i1][j1]);
  }

printf("...done");

//----------------CALCULATION----------------------------------------


printf("\nCalculation : t= %lf: tmax= %lf\n", t,tmax);

ih = (int) (1/h);
TAU = h;

double print_time=0.0;

while( t < tmax){
  tmax1 = t + TAU;
  ii++;
  rk(N_EQ, fun, h, t, y_ini, f_ini, y1, y2);
  t = tmax1;

  if (t>print_time){
    printf("Running simulation time: %lf msec \n", t);
    print_time=t+1000;
  }

  //************************************************************************
  // Scale conductances

  if( (t>50) && (t<(50+1.5*h)) ){
    printf("Scale conductance with surface area...\n");

    if(I_CX == 1){
      for(j = 0; j < Mcx1; ++j)
        for(i = 0; i < Mcx; ++i){
          g_ampa_cx_cx[i][j] = g_AMPA_CX_CX / cx_cell[0][0].S_CX_DEND;
          g_gaba_a_in_cx[i][j] = g_GABA_A_IN_CX / cx_cell[0][0].S_CX_DEND;

	  // subnetwork AMPA selection
	  if (networks_ampa[0][0] != 0) {
	      if (j!=0) { printf("Code not ready for 2D.\n"); return 1; } //not to forget
	      int cluster = i/Mcx_cluster_size; // i/50
	      cluster++; //matlab indexing of ROIs
	      double ampa = network_get_node_ampa(cluster); // rets 1 by default
	      g_ampa_cx_cx[i][j] *= ampa;
	  } //subnetwork AMPA

        }
      g_ampa_cx_cx_MINICE = g_AMPA_CX_CX_MINICE / cx_cell[0][0].S_CX_DEND;
      g_nmda_cx_cx = g_NMDA_CX_CX / cx_cell[0][0].S_CX_DEND;
      g_ext_cx1 = g_Extern_ampa/cx_cell[0][0].S_CX_DEND;
    }

    if(I_IN == 1){
      for(j = 0; j < Min1; ++j)
        for(i = 0; i < Min; ++i){		
          g_ampa_cx_in[i][j] = g_AMPA_CX_IN / in_cell[0][0].S_CX_DEND;
          g_ampa_cx_in_MINICE = g_AMPA_CX_IN_MINICE / in_cell[0][0].S_CX_DEND;
          g_nmda_cx_in = g_NMDA_CX_IN / in_cell[0][0].S_CX_DEND;
          g_gaba_a_in_in = g_GABA_A_IN_IN / in_cell[0][0].S_CX_DEND;
          g_ext_in1 = g_Extern_ampa1/in_cell[0][0].S_CX_DEND;
        }
    }
    printf("...done.");
  }

  //************************************************************************
  // External input

  if( (t>5) && (t<(5+1.5*h)) ){
    for(i = 0; i < Min; ++i)
      for(j = 0; j < Min1; ++j){
        a_ext44[i][j].w = 0.12; //0.1;
      } 
    for(i = 0; i < Mcx; ++i)
      for(j = 0; j < Mcx1; ++j){
        a_ext33[i][j].w = 0.11; //0.1;
      }
  }

  /// Synaptic depression
  if(t>1000) {
    // SpikeCount
    for(j = 0; j < Mcx1; ++j)
      for(i = 0; i < Mcx; ++i){	    
        if (cx_cell[i][j].v_SOMA > 10 && (t-prev[i][j])> 2){
	  fprintf(fsp,"0 %d %d %lf\n",i,j,t);
          spikecount[i][j]= spikecount[i][j]+1;
          if ( ((spikecount[i][j])/(Check_time/1000)) < 100 )
            E[i][j] = 1 - (1 - E[i][j]*(1-0.07)) * exp(-(t-prev[i][j])/700);

          prev[i][j] = t;
        }
      }
    
    for(j = 0; j < Min1; ++j)
      for(i = 0; i < Min; ++i){
        if (in_cell[i][j].v_SOMA > 10 && (t-prev_in[i][j])> 2){
	  fprintf(fsp,"1 %d %d %lf\n",i,j,t);
          spikecount_in[i][j]= spikecount_in[i][j]+1;
          prev_in[i][j] = t;
        }
      }
  }
        

  if(cx_cell[Mcx/2][0].CX_DEND::Ko>6.0) {
    t=tmax+10;
    printf("\n Seizure !!!");
  }

  // Compute firing rate and export 
  if((t>10000) && (ii/(ih*Check_time))*(ih*Check_time) == ii) {
    average_rate_py = 0;
    average_rate_in = 0;

    fprintf(f_in_rates,"%lf \t",t);
    fprintf(f_py_rates,"%lf \t",t);

    for(j = 0; j < Min1; ++j)
      for(i = 0; i < Min; ++i){
        fprintf(f_in_rates,"%lf \t",spikecount_in[i][j]);
        average_rate_in = average_rate_in + spikecount_in[i][j];
        spikecount_in[i][j]=0;
      }
    average_rate_in = (average_rate_in/Min)/(Check_time/1000);
    fprintf(f_in_rates,"%lf \n",average_rate_in);

    for(j = 0; j < Mcx1; ++j)
      for(i = 0; i < Mcx; ++i){
        fprintf(f_py_rates,"%lf \t",spikecount[i][j]);
        average_rate_py = average_rate_py + spikecount[i][j];
        spikecount[i][j]=0;
      }

    average_rate_py = (average_rate_py/Mcx)/(Check_time/1000);
    fprintf(f_py_rates,"%lf \n",average_rate_py);

    fflush(f_py_rates);
    fflush(f_in_rates);
  }

  // Output 
  if((t > ttime) && ((ii/(2500))*(2500) == ii) ){

    fprintf(f_cx,"%lf \t",t);
    for(i = 0; i < Mcx; ++i) {
      fprintf(f_cx,"%lf \t", cx_cell[i][0].v_SOMA);
      fprintf(f_cx_syn,"%lf ", cx_cell[i][0].tot_syn_exc);
      cx_cell[i][0].tot_syn_exc=0; //total excitatory synaptic input summed for the whole dt
      fprintf(f_cx_syn_inh,"%lf ", cx_cell[i][0].tot_syn_inh);
      cx_cell[i][0].tot_syn_inh=0;

    } 
    fprintf(f_cx,"\n");
    fprintf(f_cx_syn,"\n");
    fprintf(f_cx_syn_inh,"\n");

    for(i = 0; i < Mcx; ++i) {
      fprintf(f1000,"%lf ", cx_cell[i][0].CX_DEND::Ko);
      fprintf(f1000,"%lf ", cx_cell[i][0].CX_DEND::Nai);
      fprintf(f1000,"%lf ", cx_cell[i][0].CX_DEND::Cli); 
      fprintf(f1000,"%lf ", y_ini[no_cx[0][0][1]]); //[Ca]
    }
    fprintf(f1000,"\n");

    for(i = 0; i < Min; ++i){
      fprintf(f_in,"%lf ", in_cell[i][0].v_SOMA);
      fprintf(f_in,"%lf ", in_cell[i][0].CX_DEND::Ko);
      fprintf(f_in,"%lf ", in_cell[i][0].CX_DEND::Nai);
    }

    fprintf(g_py_py,"\n");
    fprintf(f_in,"\n");
    fprintf(g_py_in,"\n");
  }

 }
//--------------------END CALCULATION-------------------------------
      
//-----------------close ALL files-----------------------------------
printf("\n GPY: %lf", GPY);
printf("\n GIN: %lf", GIN);
printf("\n GPYPY: %lf", GPYPY);
printf("\n GPYIN: %lf", GPYIN);
      
printf("\n excounter: %d", excounter);
fclose(f_cx);
fclose(f_cx_syn);
fclose(f_cx_syn_inh);
fclose(f_in);
fclose(f1000);
fclose(g_py_in);
fclose(log_file);
printf("\n");

return 0;
}

void index_loc(int index, int *type, int *x, int *y){  
  int i = 0;
  int j = 0;

  for(i = 0; i<Num_Types; i++){
    index = index - (cell_sizes[i][0]*cell_sizes[i][1]);
    if(index < 0){
      break;
    }
  }

  index = index + (cell_sizes[i][0]*cell_sizes[i][1]);
  for(j = 0; j<cell_sizes[i][0]; j++){
    index = index - cell_sizes[i][1];
    if(index <0){
      break;
    }
  }

  index = index + cell_sizes[i][1];
  *type = i;
  *x = j;
  *y = index;

}

void index_array(int ***array){

  int i = 0;
  int type;
  int x;
  int y;
  (*array) = new int*[Num_Cells];
  
  for(i = 0; i<Num_Cells; i++){
    (*array)[i] = new int[3];
    index_loc(i,&type,&x,&y);
    (*array)[i][0] = type;
    (*array)[i][1] = x;
    (*array)[i][2] = y;
  }
}

//+++++++++++ Function to calculate the right sides for ALL ODE +++++++++++++++  
void fun(unsigned neq, double x, double *y_ini, double *f_ini){
  int i, j, k, kmax, i1, j1, ii, jj;

  int tid;

#pragma omp parallel shared(x) private(tid,i,j)  num_threads(NUM_THREADS)
  {
    //========here the MAIN loop to calculate intrinsic and synaptic conductances========= 
    //--------(f_ini IS changed, y_ini IS NOT changed)-------------------------      
    //printf("num threads %d \n",omp_get_num_threads());                                                                                                                        
#pragma omp for schedule(guided) nowait
    for(i=0; i<Num_Cells; i++){
      calcit(i, neq,  x,  y_ini, f_ini);
    }

#pragma omp barrier

#pragma omp for schedule(guided) nowait
    for(i=0; i<Num_Cells; i++){
      calc_syn(i, neq,  x,  y_ini, f_ini);
    }

  }

}

//+++++++++++ Function to calculate the right sides for ALL ODE +++++++++++++++
void calcit(int index,unsigned neq, double x, double *y_ini, double *f_ini){
  //double scale;
  //int i, j, k, kmax, ii1, ii2, k1, i1, j1, ii, jj, nst, kk[Max][Max1], kk1[Max][Max1];
  int i, j, ii1, ii2;
  //int k_ampa, k_nmda;
  double Kos_tmp, Naos_tmp, Kod_tmp, Naod_tmp;
  double *tmpa,tmp;
  
  switch(index_loc_array[index][0]){

    ///////////////////// case cx
  case 0:

    i = index_loc_array[index][1];
    j = index_loc_array[index][2];

    // Assuming Mcx > 3
    ii1=i-1;
    ii2=i+1;
    if(ii1 == -1)ii1=0;
    if(ii2 == Mcx)ii2=Mcx-1;


    // If neuron is edge of the network
    ii1=i-1;
    ii2=i+1;
    if (Mcx==1) {
      ii1=i;
      ii2=i; }
    else if (Mcx==2) {
      if(ii1 == -1) {
        ii1=ii2; }
      else{
        ii2=ii1;}}
    else{
      if(ii1 == -1)ii1=0;
      if(ii2 == Mcx)ii2=Mcx-1;}
      
    if (i % Mcx_cluster_size == 0 & i>0){
      ii2=ii1; }

    if (i % (Mcx_cluster_size+1) == 0 & i>0){
      ii1=ii2; }
  
    Kod_tmp = ((y_ini[no_cx[ii1][j][14]]) + (y_ini[no_cx[ii2][j][14]]))/2;
    Naod_tmp =((y_ini[no_cx[ii1][j][12]]) + (y_ini[no_cx[ii2][j][12]]))/2;
    Kos_tmp = ((y_ini[N_DEND+no_cx[ii1][j][7]]) + (y_ini[N_DEND+no_cx[ii2][j][7]]))/2;
    Naos_tmp = ((y_ini[N_DEND+no_cx[ii1][j][5]])+ (y_ini[N_DEND+no_cx[ii2][j][5]]))/2;

    cx_cell[i][j].calc(x, y_ini+no_cx[i][j][0], f_ini+no_cx[i][j][0], 
                       Kos_tmp, Naos_tmp, Kod_tmp, Naod_tmp, 
                       ga_in_cx[i][j].I ); 

        break;
    
  case 1:
    i = index_loc_array[index][1];
    j = index_loc_array[index][2];
    ii1=i-1;
    ii2=i+1;

    // If neuron is edge of the network
    if (Min==1) {
      ii1=i;
      ii2=i; }
    else if (Min==2) {
      if(ii1 == -1) {
        ii1=ii2; }
      else{
        ii2=ii1;}}
    else{
      if(ii1 == -1)ii1=0;
      if(ii2 == Min)ii2=Min-1;}

    if (i % Min_cluster_size == 0 & i>0){
      ii2=ii1; }

    if (i % (Min_cluster_size+1) == 0 & i>0){
     ii1=ii2; }

    Kod_tmp = ((y_ini[no_in[ii1][j][14]]) + (y_ini[no_in[ii2][j][14]]))/2;
    Naod_tmp = ((y_ini[no_in[ii1][j][12]]) + (y_ini[no_in[ii2][j][12]]))/2;
    Kos_tmp = ((y_ini[N_DEND+no_in[ii1][j][7]]) + (y_ini[N_DEND+no_in[ii2][j][7]]))/2;
    Naos_tmp = ((y_ini[N_DEND+no_in[ii1][j][5]])+ (y_ini[N_DEND+no_in[ii2][j][5]]))/2;

    in_cell[i][j].calc(x, y_ini+no_in[i][j][0], f_ini+no_in[i][j][0], 
                       Kos_tmp, Naos_tmp, Kod_tmp, Naod_tmp, 
                       ga_in_in[i][j].I);

    break;

  }
}

void calc_syn(int index,unsigned neq, double x, double *y_ini, double *f_ini){
  
  int i, j, k, i1, j1;
  int k_ampa, k_nmda;

  
  i = index_loc_array[index][1];
  j = index_loc_array[index][2];

  double total_syn_input_exc = 0; //total excitatory/ihibitory input
  double total_syn_input_inh = 0;

  switch(index_loc_array[index][0]){

    ///////////////////// case cx
  case 0:

    if(I_CX == 1){
      for(i1 = 0, k_ampa = 0,k_nmda = 0; i1 < Mcx; ++i1){
        for(j1 = 0; j1 < Mcx1; ++j1){
          if((C_CXCX[i1][j1][i][j] ==1)){
            a_cx_cx[i][j].calc(g_ampa_cx_cx[i][j]/k_CXCX[i][j], g_ampa_cx_cx_MINICE,
                               x, y_ini[no_cx[i][j][0]],
                               cx_cell[i1][j1].v_SOMA,E[i1][j1],i1,j1);
            ++k_ampa; 
          }
	  
          if((C_CXCX[i1][j1][i][j] ==2)){
            a_cx_cx[i][j].calc(Long_range_cxcx_strength*Cluster_strength_CXCX[i1][j1][i][j]*g_ampa_cx_cx[i][j]/k_CXCX[i][j],
                               g_ampa_cx_cx_MINICE,
                               x, y_ini[no_cx[i][j][0]],
                               cx_cell[i1][j1].v_SOMA,E[i1][j1],i1,j1);
            ++k_ampa; 
          }

          if(C_CXCX[i1][j1][i][j] ==1){
            nmda_cx_cx[i][j].calc(g_nmda_cx_cx/k_CXCX[i][j], x, y_ini[no_cx[i][j][0]],
                                  cx_cell[i1][j1].v_SOMA,i1,j1);
            ++k_nmda; 
          }
        }
      }
      
      if(k_ampa > 0) {
        f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - a_cx_cx[i][j].I;
	total_syn_input_exc -= a_cx_cx[i][j].I;
      }

      if (i ==5) GPYPY = GPYPY + a_cx_cx[i][j].Rout;
      
      if (k_nmda > 0) {
        nmda_cx_cx[i][j].I = ( y_ini[no_cx[i][j][0]]/(1+exp(-(y_ini[no_cx[i][j][0]] - (-25))/12.5)) ) *nmda_cx_cx[i][j].I;
        f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - nmda_cx_cx[i][j].I;
        total_syn_input_exc -= nmda_cx_cx[i][j].I;
      }
    }
    
    if (x > 0) { // with external input
      firstflag = 1;
      if(I_CX == 1) {
	
        a_ext33[i][j].calc(g_ext_cx1, x);
        if (i == 5){
          GPY = GPY + a_ext33[i][j].g;
          excounter++;
        }
        f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - a_ext33[i][j].g * y_ini[no_cx[i][j][0]];
	total_syn_input_exc -= a_ext33[i][j].g * y_ini[no_cx[i][j][0]];
      }
    
      ///****** CONSTANT Ext inp ********///
      /* if( (x > 40000) && (x < 41000) ) {
      //for(i = 0; i < Mcx; ++i)
      // for(j = 0; j < Mcx1; ++j){
      // a_ext33[i][j].calc(g_ext_cx1, x);
      f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] + 0.35;
      } */

      if( (x > input_start_time) && (x < input_end_time) ) {
         //for(i = 0; i < Mcx; ++i)
         // for(j = 0; j < Mcx1; ++j){
         // a_ext33[i][j].calc(g_ext_cx1, x);
        if (i<Mcx_cluster_size) // apply input to one cluster
          f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] + p_extinp;
	  total_syn_input_exc += p_extinp;
      }  
      
    }else{
      firstflag = 0;
    }   
    
    //--------------GABA-A from IN to CX cells------------------------------------
    if(I_CX == 1 && I_IN == 1){
      for(i1 = 0, k = 0; i1 < Min; ++i1)
        for(j1 = 0; j1 < Min1; ++j1){
          if(C_INCX[i1][j1][i][j] > 0){
            my_E_GABA[i][j] = 26.639 * log(y_ini[no_cx[i][j][16]]/130.0); //inserted
            ga_in_cx[i][j].calc(g_gaba_a_in_cx[i][j]/k_INCX[i][j],x,y_ini[no_cx[i][j][0]],
                                in_cell[i1][j1].v_SOMA,i1,j1, my_E_GABA[i][j]);
            ++k;
          }
        }
      
      if(k > 0) {
        f_ini[no_cx[i][j][0]] = f_ini[no_cx[i][j][0]] - ga_in_cx[i][j].I;
	total_syn_input_inh -= ga_in_cx[i][j].I;
      }
    }

    cx_cell[i][j].tot_syn_exc += total_syn_input_exc;
    cx_cell[i][j].tot_syn_inh += total_syn_input_inh;
    break;
  
    //////////////////////////////    
  case 1:
    //--------------GABA-A from IN to IN cells------------------------------------
    if(I_IN == 1){
      for(i1 = 0, k = 0; i1 < Min; ++i1)   
        for(j1 = 0; j1 < Min1; ++j1){	 
          if(C_ININ[i1][j1][i][j] > 0){
            my_E_GABA1[i][j] = 26.639 * log(y_ini[no_in[i][j][16]]/130.0);
            ga_in_in[i][j].calc(g_gaba_a_in_in/k_ININ[i][j],x,y_ini[no_in[i][j][0]],
                                in_cell[i1][j1].v_SOMA, i1,j1, my_E_GABA1[i][j]);
            ++k;
          }
        }
      if(k >  0){
        f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - ga_in_in[i][j].I; 
      }

      if( (x > input_start_time) && (x < input_end_time) ) {
	f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] + p_extinp;
      }  

    }

    //-------------AMPA from CX to IN cells--------------------------------------
    if(I_CX == 1 && I_IN == 1){
      for(i1 = 0, k_ampa = 0, k_nmda = 0; i1 < Mcx; ++i1)
        for(j1 = 0; j1 < Mcx1; ++j1){
          if(C_CXIN[i1][j1][i][j] == 1){
            a_cx_in[i][j].calc(g_ampa_cx_in[i][j]/k_CXIN[i][j], g_ampa_cx_in_MINICE,
                               x, y_ini[no_in[i][j][0]],
                               cx_cell[i1][j1].v_SOMA,E[i1][j1],i1,j1);
            ++k_ampa; 
          }

          if(C_CXIN[i1][j1][i][j] == 2){
            a_cx_in[i][j].calc(Long_range_cxin_strength*Cluster_strength_CXIN[i1][j1][i][j]*g_ampa_cx_in[i][j]/k_CXIN[i][j],
                               g_ampa_cx_in_MINICE,
                               x, y_ini[no_in[i][j][0]],
                               cx_cell[i1][j1].v_SOMA,E[i1][j1],i1,j1);
            ++k_ampa; 
          }

          if(C_CXIN[i1][j1][i][j] == 1){
            nmda_cx_in[i][j].calc(g_nmda_cx_in/k_CXIN[i][j], x, y_ini[no_in[i][j][0]],
                                  cx_cell[i1][j1].v_SOMA,i1,j1);
            ++k_nmda; 
          }
        }
      
      
      if(k_ampa > 0) {
        f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - a_cx_in[i][j].I; 
      }
      if (i ==1){
        GPYIN = GPYIN + a_cx_in[i][j].Rout;
      }
      if (k_nmda >0){
        nmda_cx_in[i][j].I = ( y_ini[no_in[i][j][0]]/(1+exp(-(y_ini[no_in[i][j][0]] - (-25))/12.5))  )*nmda_cx_in[i][j].I;
        f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] - nmda_cx_in[i][j].I;
      }
    }
    if (x > 0) { // with external input
      firstflag = 1;
      if(I_IN == 1) {
        a_ext44[i][j].calc(g_ext_in1, x);
        f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] -
          a_ext44[i][j].g * y_ini[no_in[i][j][0]];
        if (i == 1) GIN = GIN + a_ext44[i][j].g;
      }
      /*
        if( (x > 50000) && (x < 55000) ) {
        //for(i = 0; i < Mcx; ++i)
        // for(j = 0; j < Mcx1; ++j){
        // a_ext33[i][j].calc(g_ext_cx1, x);
        f_ini[no_in[i][j][0]] = f_ini[no_in[i][j][0]] + 0.6;
        }*/

    }else{
      firstflag = 0;
    }    
    break;
  }
}

//***************************************************************************
// Solver
//***************************************************************************
void rk(unsigned n, void fun(unsigned, double, double*, double*),
        double h, double x, double* y, double* f, double* s, double* yk)
{
  int i;
  double xk;
  double h_half = h/2.;


  fun(n, x, y, f);

  for(i = n-1; i >=0; --i){
    s[i] = f[i];
    yk[i] = y[i] + h_half*f[i];
  }

  xk = x + h_half;
  fun(n, xk, yk, f);

  for(i = n-1; i >=0; --i){
    s[i] += 2.*f[i];
    yk[i] = y[i] + h_half*f[i];
  }
  fun(n, xk, yk, f);

  for(i = n-1; i >=0; --i){
    s[i] += 2.*f[i]; yk[i] = y[i] + h*f[i];
  }

  xk = x + h;
  fun(n, xk, yk, f);

  for(i = n-1; i >=0; --i){
    y[i] += (h/6.)*(s[i] + f[i]);
  }
}

//***************************************************************************
double gaussrand()
{
  static double U, V;
  static int phase = 0;
  double Z;

  if(phase == 0) {
    U = (rand() + 1.) / (RAND_MAX + 2.);
    V = rand() / (RAND_MAX + 1.);
    Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
  } else
    Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

  phase = 1 - phase;

  return Z;
}

