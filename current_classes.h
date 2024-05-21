//++++++++++++++ CURRENTs DESCRIPTION ++++++++++++++++++++++++++++++++++++++
#define time2(a) (1)


//-------------- Ca-dependent potassium current (RE cell) --------------------
class IKCa_RE {
 static double Alpha, Beta, Cels,Tad_inv, AlphaBeta;
 double m_inf, tau_m,  car;
 public:
 double iKCa, m0, m_KCa;
 double G_KCa, g_KCa;
 IKCa_RE(double cai) {
  G_KCa = 0.;
  car = Alpha*cai*cai/Beta;
  m0 = car/(car + 1);  }
 void calc(double, double&, double, double, double);
};

double IKCa_RE::Alpha = 48, IKCa_RE::Beta = 0.03, IKCa_RE::Cels = 36, IKCa_RE::Tad_inv = 1.0/pow(3,((Cels-22)/10)), IKCa_RE::AlphaBeta = Alpha/Beta;

void IKCa_RE::calc(double m, double &fm, double v, double cai, double eK){
 m_KCa=m;
 g_KCa= G_KCa*m*m;
 iKCa = g_KCa*(v - eK);
 car = AlphaBeta*cai*cai;
 m_inf = car/(car + 1);
 tau_m = (Tad_inv/(Beta*(car + 1)));
 if(tau_m < 0.1) tau_m = 0.1;
 fm = -(1/tau_m)*(m - m_inf);

 if (fm != fm)
  printf("\n m= %lf, tau_m= %f, cai = %lf",m,tau_m,cai);
}


//------------------Ca-dynamics------------------------------------
class ICa {
 static double Ca_inf;
 double drive, drive0;
 public:
 double Taur, D;
 ICa() {Taur = 5; D = 1.; //0.1;
  drive0 = 10.0/(2.*96489.); }
 void calc(double cai, double &fcai, double iT);
};

double ICa::Ca_inf = 2.4e-4;

void ICa::calc(double cai, double &fcai, double iT) {
 drive = -drive0 * iT / D;
 if(drive < 0) drive = 0;
 fcai = drive + (Ca_inf - cai)/Taur;
}


//---------------------Hight-threshold Ca2+ current (CX cell)----------------
class IHVA_CX {
 static double Ca_0, Cels, Qm, Qh, E_Ca;
 double m_inf, tau_m, h_inf, tau_h, Phi_m, Phi_h, a, b, vm;
 double ratio, eca0, eca;
 public:
 double iHVA, m0, h0;
 double G_HVA, g_HVA, m_HVA, h_HVA;
 IHVA_CX(double v) {
  G_HVA = 0.03;
  Phi_m = pow(Qm,((Cels-23)/10));
  Phi_h = pow(Qh,((Cels-23)/10));
  a = 0.055*(-27 - v)/(exp((-27-v)/3.8) - 1);
  b = 0.94*exp((-75-v)/17);
  m0 = a/(a+b);
  a = 0.000457*exp((-13-v)/50);
  b = 0.0065/(exp((-v-15)/28) + 1);
  h0 = a/(a+b);

 }
 void calc(double m, double h, double &fm, double &fh, double v, double cai);
};

double IHVA_CX::Ca_0 = 2, IHVA_CX::E_Ca = 140;
double IHVA_CX::Qm = 2.3, IHVA_CX::Qh = 2.3, IHVA_CX::Cels = 36;

void IHVA_CX::calc(double m, double h, double &fm, double &fh,
                   double v, double cai) {
 m_HVA=m;
 h_HVA=h;

 g_HVA = Phi_m * G_HVA * m*m*h;
 iHVA = g_HVA * (v - E_Ca);

 a = 0.055*(-27 - v)/(exp((-27-v)/3.8) - 1);
 b = 0.94*exp((-75-v)/17);
 tau_m = (1/(a+b))/Phi_m;
 m_inf = a/(a+b);

 a = 0.000457*exp((-13-v)/50);
 b = 0.0065/(exp((-v-15)/28) + 1);
 tau_h = (1/(a+b))/Phi_h;
 h_inf = a/(a+b);

 fm = -(m - m_inf)/tau_m;
 fh = -(h - h_inf)/tau_h;
}


//--------------------Potassium M-current (CX cell)-------------------------
class IKm_CX {
 static double E_Km, Ra, Rb, Cels, Q, tha, qa;
 double Tad, a, b;
 public:
 double iKm, m0;
 double G_Km, g_Km;
 IKm_CX(double v) {
  G_Km = 0;
  Tad = pow(Q,((Cels-23)/10));
  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  m0 = a/(a+b);
 }
 void calc(double, double&, double, double);
};

double IKm_CX::E_Km = -90, IKm_CX::Q = 2.3;
double IKm_CX::tha = -30, IKm_CX::qa = 9;
double IKm_CX::Ra = 0.001, IKm_CX::Rb = 0.001, IKm_CX::Cels = 36;

void IKm_CX::calc(double m, double &fm, double v, double eK){
 g_Km= Tad * G_Km * m; 
 iKm = g_Km * (v - eK);
 a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
 b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
 fm = -(Tad*(a+b))*(m - (a/(a+b)));
}

//--------------------Fast potassium current (CX cell)-------------------
class IKv_CX {
 static double Ra, Rb, Cels, Q, tha, qa;
 double Tad, a, b;
 public:
 static double E_Kv;
 double iKv, g_Kv, m0;
 double G_Kv;
 IKv_CX(double v) {
  G_Kv = 150;
  Tad = pow(Q,((Cels-23)/10));

  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
  m0 = a/(a+b);
 }
 void calc(double, double&, double, double);
};

double IKv_CX::E_Kv = -90, IKv_CX::Q = 2.3;
double IKv_CX::tha = 25, IKv_CX::qa = 9;
double IKv_CX::Ra = 0.02, IKv_CX::Rb = 0.002, IKv_CX::Cels = 36;

void IKv_CX::calc(double m, double &fm, double v,double eK){
 g_Kv = Tad * G_Kv * m;
 iKv = g_Kv * (v - eK);
 a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa));
 b = -Rb * (v - tha) / (1 - exp((v - tha)/qa));
 fm = -(Tad*(a+b))*(m - (a/(a+b)));
}


//----------------Fast sodium current (CX cell)--------------------------
class INa_CX {
 static double Ca_0, Cels, Qm, Qh;
 static double tha, qa, Ra, Rb, thi1, thi2, qi, thinf, qinf, Rg, Rd;
 double h_inf, Phi_m, Phi_h, a, b, vm;
 inline double trap0(double v, double th, double a, double q) {
  if (fabs(v/th) > 1.0e-6) {
   return ( a * (v - th) / (1 - exp(-(v - th)/q)) );
  } else {
   return (a * q ); }
 }
 public:
 static double E_Na;
 double iNa, g_Na, m0, h0, m_Na, h_Na;
 double G_Na;
 INa_CX(double v) {
  G_Na = 3000;
  Phi_m = pow(Qm,((Cels-23)/10));
  Phi_h = pow(Qh,((Cels-23)/10));
  vm = v -10;
  a = trap0(vm,tha,Ra,qa);
  b = trap0(-vm,-tha,Rb,qa);
  m0 = a/(a+b);

  a = trap0(vm,thi1,Rd,qi);
  b = trap0(-vm,-thi2,Rg,qi);
  h0 = 1/(1+exp((vm-thinf)/qinf));
 }
 void calc(double, double, double&, double&,
           double, double);
};

double INa_CX::E_Na = 50;
double INa_CX::Qm = 2.3, INa_CX::Qh = 2.3, INa_CX::Cels = 36;
double INa_CX::tha = -35, INa_CX::qa = 9;
double INa_CX::Ra = 0.182,INa_CX::Rb = 0.124;
double INa_CX::thi1 = -50, INa_CX::thi2 = -75, INa_CX::qi = 5;
double INa_CX::thinf = -65, INa_CX::qinf = 6.2;
double INa_CX::Rg = 0.0091, INa_CX::Rd = 0.024;

void INa_CX::calc(double m, double h, double &fm, double &fh,
                  double v, double eNa) {
 m_Na=m;
 h_Na=h;
  
 g_Na = Phi_m * G_Na * m*m*m*h;
 iNa = g_Na * (v - eNa);
 vm = v -10;

 a = trap0(vm,tha,Ra,qa);
 b = trap0(-vm,-tha,Rb,qa);
 fm = -(m - a/(a+b))*(a+b)*Phi_m;

 a = trap0(vm,thi1,Rd,qi);
 b = trap0(-vm,-thi2,Rg,qi);
 h_inf = 1/(1+exp((vm-thinf)/qinf));

 fh = -(h - h_inf)*(a+b)*Phi_h;
}

//-------------------Persist Sodium current (CX cell)-------------------
class INap_CX {
 double Tet, Sig, f, Cels, Q10, Phi_m, tau_m, m_inf;
 public:
 double m0, iNap, G_Nap, g_Nap;
 INap_CX(double v) {
  G_Nap = 2;
  Tet = -42;
  Sig = 5;
  f = 0.02;
  Cels = 36;
  Q10 = 2.7;
  Phi_m = pow(Q10,((Cels-22)/10));
  tau_m = 0.8/Phi_m;
  m0 = f/(1 + exp(-(v - Tet)/Sig));
 }
 void calc(double, double&, double, double);
 double cond();
};

void INap_CX::calc(double m, double &fm, double v, double eNa){
 g_Nap = G_Nap * m;
 iNap = g_Nap * (v - eNa); //EE_Na);
 m_inf = f/(1 + exp(-(v - Tet)/Sig));
 fm = -(m - m_inf)/tau_m;
}
double INap_CX::cond(){
 return(m_inf);
}

//----------------- h-current (CX cell) ---------------------------------
class Ih_CX {
 static double Cels;
 double h_inf, tau_s1, tau_s2, alpha, beta, Phi;
 public:
 double ih, m01, G_h, E_h, m_h;

 Ih_CX(double v) {
  G_h = 0.02;
  E_h = -40;
  Phi = pow(3,((Cels-36)/10));
  h_inf = 1/(1 + exp((v + 78)/7));
  tau_s1 = 38;
  m01=h_inf;
 }
 void calc(double, double&, double, double);
};

double Ih_CX::Cels = 36;

void Ih_CX::calc(double m, double &fm1,
                 double v, double eH) {
 m_h=m;
 ih = G_h*m*(v - eH); //E_h);
 h_inf = 1/(1 + exp((v + 82)/7));
 fm1 = -(m - h_inf)/tau_s1;
}


//------------------Na-K-dynamics------------------------------------
class IKNa {
 double driveNai, driveKi, drive0, A, dr;
 public:
 double D, DD, DD1, DD2, DD3, I_Kpump, I_Napump, Bo,  k1, k2, Koc2, Imax, k1N, k2N, ex_space;
 static double Nai0, Nao0, Ki0, Ko0, Ki1, Nai1, Nao1;
 IKNa() {D = 1.;
  drive0 =  10.0/(96489.);
  k1=0.0008;
  Bo=500;
  k1N = 1.05;
  k2N = 1.05;
  Koc2 = 15;
  Imax = 20; //5;
  ex_space=0.15;
  DD=6.0e-6/0.0001; // K diffusion outside 6.0e-6/0.0001;
  DD1=6.0e-6/0.0001; // K, Na - intracompartmental outside
  DD3=6.0e-6/0.0001; // Na diffusion outside
  DD2=6.0e-6/0.0001; // K, Na - soma <--> dendrite intracellular
        
  Nai0 = 20; 
  Nao0 = 130;
  Ko0=2.5; 
  Ki0=130;
 }
 void calc(double, double, double, double, double,
           double&, double&, double&, double&, double&,
           double, double, double, 
           double, double, double, double, double, double);
};

double IKNa::Nai0 = 20, IKNa::Nao0 = 130, IKNa::Ko0=2.5, IKNa::Ki0=130;

void IKNa::calc(double Nai, double Nao, double Ki, double Ko, double B,
                double &fNai, double &fNao, double &fKi, double &fKo, double &fB,
                double iNatotal, double iKtotal, double x, 
                double Ki1, double Nai1, double Ko1, double Ko2, double Nao1, double Nao2) {

 A = (1/((1+Ko0/Ko)*(1+Ko0/Ko))) * (1/((1+Nai0/Nai)*(1+Nai0/Nai)*(1+Nai0/Nai)));
 I_Kpump = -2*Imax*A;
 I_Napump = 3*Imax*A;

 driveNai = -drive0 * (iNatotal + I_Napump) / D; 
 fNai = driveNai + DD2*(Nai1-Nai);
 fNao = -driveNai/ex_space + DD1*(Nao1-Nao) + DD3*(Nao2-Nao);

 k2= k1/(1+exp((Ko-Koc2)/(-1.05)));

 driveKi = -drive0 * (iKtotal + I_Kpump) / D; 
 fKi = driveKi + DD2*(Ki1-Ki);

 fKo = -1*driveKi/ex_space + k1*(Bo-B)-k2*Ko*B + DD1*(Ko1-Ko) + DD*(Ko2-Ko); 
 fB = (k1*(Bo-B) - k2*Ko*B); 
}


//----------------- Na sensitive K-current (CX cell) ---------------------------------

class IK_Na {
 public:
 double G_KNa, ikna, g_KNa, Nai;
 IK_Na() {
  G_KNa = 1.33; 
 }
 void calc(double Nai, double v, double eK, double x);
};

void IK_Na::calc(double Nai, double v, double eK, double x) {
 g_KNa = G_KNa * 0.37 / ( 1 + (pow((38.7*2/Nai),3.5) ) );
 ikna =  g_KNa * (v - eK);
}


class ICl {
 double drive, drive0, taucl;
 public:
 double Cl_inf, Cl_inf1, Tau_inf, Cl_Tauinf, Cl_drive;
 double Taur, D;
 ICl() {Taur = 100000; D = 1.;
  drive0 = Cl_drive/(1.*96489.); }
 void calc(double cli, double &fcli, double iCl, double x);
};

void ICl::calc(double cli, double &fcli, double iCl, double Ko) {
 taucl= (100 + Cl_Tauinf/(1+exp((Cl_inf1-Ko)/Tau_inf)));
 drive0 = Cl_drive/(1.*96489.);
 drive = drive0*iCl/D + (Cl_inf-cli)/(taucl);
 fcli = drive;
}


//-------------------CX CELL------------------------------------------------
//-------------------CX CELL (DENDRITE)-------------------------------------
class CX_DEND: public IHVA_CX, public IKCa_RE, public IKm_CX, public INa_CX, public INap_CX, public ICa, public Ih_CX, public ICl, public IKNa 
{
 double ratio, e0, eNa1, eK1, eH1, eLK1, A;
 int index;
 public:
 double iDEND, I_Stim1, E_l, G_l, G_kl, G_Nal, E_K, Nai, Nao, Ki, Ko, Ki1, Nai1, Nao1, eNa, eK, B, iNatotal, iKtotal, IKl, INal, eLK, eH, mycalc, Cli;
 CX_DEND(double V0, double Cai0) :IHVA_CX(V0), IKCa_RE(Cai0), IKm_CX(V0),
  INa_CX(V0), INap_CX(V0), ICa(), Ih_CX(V0), ICl(), IKNa() {
  E_l = -74.273;
  E_K = -95;
  G_kl = 0;
  IKNa::D = 1;
  INa_CX::G_Na = 0.0;
  I_Stim1 = 0;
  ICa::Taur = 500;
  INap_CX::G_Nap = 0.0;
  IKm_CX::G_Km = 0;
  IKCa_RE::G_KCa = 0.0;
  IHVA_CX::G_HVA = 0.0;
  Ih_CX::G_h = 0;
  G_l = 1.0e3/30000; 
  e0 = 1000*8.31441*(273.15 + 36)/(96489);
  index=0;

  IKNa::Koc2 = 15;
  IKNa::ex_space=0.15;
  Imax = 20;

  k1N =1.0;
  k2N = k1N;
 }
 void init(double V0, double Cai0, double *y) {
  y[0] = V0;
  y[1] = Cai0;
  y[2] = IHVA_CX::m0;
  y[3] = IHVA_CX::h0;
  y[4] = IKCa_RE::m0;
  y[5] = IKm_CX::m0;
  y[6] = INa_CX::m0;
  y[7] = INa_CX::h0;
  y[8] = INap_CX::m0;
  y[9] = Ih_CX::m01;
  //y[10] = Ih_CX::m02;
  y[11] = IKNa::Nai0;
  y[12] = IKNa::Nao0;
  y[13] = IKNa::Ki0;
  y[14] = IKNa::Ko0;
  y[15] = IKNa::Bo;
  y[16] = 5.0; //8

 }
 void calc(double, double*, double*, double, double, double, double, double, double, double);
};


void CX_DEND::calc(double x, double *y, double *f, 
                   double Ki1, double Nai1, double Ko1, double Ko2, double Nao1, double Nao2, 
                   double TotalInhib){
 // Ionic concentration dynamics
 Nai=y[11];
 Nao=y[12];
 Ki=y[13];
 Ko=y[14];
 B=y[15];
 Cli=y[16];

 ratio = Nao/Nai;
 eNa = e0 * log(ratio);
 if(ratio <= 0.)
  printf("\n LOG ERROR: Na-DEND: t=%lf Nai=%lf, Ko=%lf iKpump=%lf", x, y[11], y[14], I_Kpump);
  
 ratio = Ko/Ki; 
 eK = e0 * log(ratio);
 if(ratio <= 0.)
  printf("\n LOG ERROR: K-DEND: t=%lf Ki=%lf, Ko=%lf iKpump=%lf", x, y[13], y[14], I_Kpump);
  
 ratio = (Ko + 0.2*Nao)/(Ki + 0.2*Nai);
 eH = e0 * log(ratio);
 if(ratio <= 0.)
  printf("\n LOG ERROR: H-DEND: t=%lf Ki=%lf, Ko=%lf iKpump=%lf", x, y[13], y[14], I_Kpump);
  
 ratio = (Cli)/(130.0);
 eLK =  e0 * log(ratio);
 if(ratio <= 0.)
  printf("\n LOG ERROR: LK-DEND: t=%lf Ki=%lf, Ko=%lf iKpump=%lf", x, y[13], y[14], I_Kpump);

 ratio = (Ko + 0.085*Nao + 0.1*Cli)/(Ki + 0.085*Nai + 0.1*130.0);
 E_l =  e0 * log(ratio); 
  
 A = (1/((1+Ko0/Ko)*(1+Ko0/Ko))) * (1/((1+Nai0/Nai)*(1+Nai0/Nai)*(1+Nai0/Nai)));
 I_Kpump = -2*Imax*A;
 I_Napump = 3*Imax*A;
  
 TotalInhib = G_l*(y[0]-eLK)+TotalInhib;

 if (G_HVA>0){
  IHVA_CX::calc(y[2], y[3], f[2], f[3], y[0], y[1]);}
 else{
  iHVA=0.0;}

 if (G_KCa>0){
  IKCa_RE::calc(y[4], f[4], y[0], y[1], eK);}
 else{
  iKCa=0.0;}
  
 if (G_KCa>0){
  ICa::calc(y[1], f[1], iHVA);}
 else{
  iHVA=0.0;}

 if (G_Nap>0){
  INap_CX::calc(y[8], f[8], y[0],eNa);}
 else{
  iNap=0.0;}  

 if (G_h>0){
  Ih_CX::calc(y[9], f[9], y[0], eH);}
 else{
  ih=0.0;}

 IKm_CX::calc(y[5], f[5], y[0], eK);
 INa_CX::calc(y[6], y[7], f[6], f[7], y[0], eNa);
 ICl::calc(y[16], f[16], TotalInhib, Ko);

 mycalc = y[1];
  
 iNatotal = G_Nal*(y[0]-eNa) + iNap + iNa;
 iKtotal =  G_kl*(y[0]-eK) + iKm + iKCa;

 IKNa::calc(y[11], y[12], y[13], y[14], y[15], f[11], f[12], f[13], f[14],
            f[15], iNatotal, iKtotal, x, 
            Ki1, Nai1, Ko1, Ko2, Nao1, Nao2);

 iDEND = -G_l*(y[0]-eLK) -G_kl*(y[0]-eK) -G_Nal*(y[0]-eNa) -I_Kpump -I_Napump -ih -iNap -iKm -iNa -iHVA -iKCa;

}

//-------------------CX CELL (SOMA)-------------------------------------
class CX_SOMA: public IKv_CX, public INa_CX, public INap_CX, public IKNa, public IK_Na{
 double e0, ratio, eNa1, eK1, eK2, A;
 public:
 double v_SOMA, iSOMA, g1_SOMA, g2_SOMA, I_Stim2, eNa, eK, Nai, Nao, Ki, Ko, B, iNatotal, iKtotal, g_l, g_Nal, eLK, g_kl, Na_sc, K_sc;
 CX_SOMA(double V0, double Cai0) :IKv_CX(V0), INa_CX(V0), INap_CX(V0), IKNa(), IK_Na()
 {
  I_Stim2 = 0;
  IKv_CX::G_Kv = 0.0;
  INa_CX::G_Na = 0.0;
  INap_CX::G_Nap = 0.0;
  e0 = 1000*8.31441*(273.15 + 36)/(96489);
  g_l = 0.05;
  g_kl = 0.05;

  IKNa::ex_space=0.15;
  IKNa::Koc2 = 15;
  Imax = 20;

  k1N = 1.0;
  k2N = k1N;
 }
 void init(double V0, double Cai0, double *y) {
  v_SOMA = V0;
  y[0] = IKv_CX::m0;
  y[1] = INa_CX::m0;
  y[2] = INa_CX::h0;
  y[3] = INap_CX::m0;
  y[4] = IKNa::Nai0;
  y[5] = IKNa::Nao0;
  y[6] = IKNa::Ki0;
  y[7] = IKNa::Ko0;
  y[8] = IKNa::Bo;
 }
 void calc(double, double*, double*, double, double, double, double, double, double);
};

void CX_SOMA::calc(double x, double *y, double *f, 
                   double Ki1, double Nai1, double Ko1, double Ko2, double Nao1, double Nao2){
 // Ion concentration dynamics
 Nai=y[4];
 Nao=y[5];
 Ki=y[6];
 Ko=y[7];
 B=y[8];
        
 ratio = Nao/Nai;
 eNa = e0 * log(ratio);
 if(ratio <= 0.)
  printf("\n LOG ERROR: Na: t=%lf Nai=%lf Nao=%lf iNa=%lf iNap=%lf iNapump=%lf ", x, y[4], y[5], iNa, iNap, I_Napump);

 ratio = Ko/Ki;
 eK = e0 * log(ratio);
 if(ratio <= 0.)
  printf("\n LOG ERROR: K: t=%lf Ki=%lf, Ko=%lf iKv=%lf iKpump=%lf",x, y[6], y[7], iKv, I_Kpump);

 A = (1/((1+Ko0/Ko)*(1+Ko0/Ko))) * (1/((1+Nai0/Nai)*(1+Nai0/Nai)*(1+Nai0/Nai)));
    
 I_Kpump = -2*Imax*A;
 I_Napump = 3*Imax*A;

 IKv_CX::calc(y[0], f[0], v_SOMA, eK);
 INa_CX::calc(y[1], y[2], f[1], f[2], v_SOMA, eNa);
 INap_CX::calc(y[3], f[3], v_SOMA, eNa);
 IK_Na::calc(Nai, v_SOMA, eK, x);

 iNatotal= g_Nal*(v_SOMA-eNa) + iNap + (iNa)/Na_sc;
 iKtotal= g_kl*(v_SOMA-eK) + ikna + (iKv)/K_sc;

 IKNa::calc(y[4], y[5], y[6], y[7], y[8], f[4], f[5], f[6], f[7], f[8],
            iNatotal, iKtotal, x, Ki1, Nai1, Ko1, Ko2, Nao1, Nao2);

 g1_SOMA = g_Nal + g_Nap + g_Na + g_kl + g_Kv + g_KNa;  
 g2_SOMA = (g_Na + g_Nal + g_Nap)*eNa + (g_Kv + g_KNa +g_kl)*eK - I_Kpump - I_Napump; 
}

//------------CX CELL (connect DENDRITE and SOMA)---------------------------
class CX: public CX_DEND, public CX_SOMA {
 static double Cai0, V0, C_inv;
 public:
 double kappa, rho, S_CX_SOMA, S_CX_DEND, kampa, knmda;
 double tot_syn_exc, tot_syn_inh; //total synaptic input
 CX() :CX_DEND(V0,Cai0), CX_SOMA(V0,Cai0) {
  kappa = 10.0e3; 
  rho = 165;
 }
 void init(double *y) {
  V0 = -68.0;
  CX_DEND::init(V0, Cai0, y);
  CX_SOMA::init(V0, Cai0, y+N_DEND);
  S_CX_SOMA = 1.0e-6;
  S_CX_DEND = S_CX_SOMA * rho;
  kampa = 0;
  knmda = 0;
  tot_syn_exc = tot_syn_inh = 0;
 }
 void calc(double, double*, double*, double, double, double, double, double);
};

double CX::Cai0 = 0.00024, CX::V0 = -68.0;
double CX::C_inv = 1.0/0.75; 

void CX::calc(double x, double *y, double *f, 
              double Ko2s, double Nao2s, double Ko2d, double Nao2d,
              double TotalInhibit){

 CX_SOMA::calc(x, y+N_DEND, f+N_DEND, 
               CX_DEND::Ki, CX_DEND::Nai, CX_DEND::Ko, Ko2s, CX_DEND::Nao, Nao2s);
    
 v_SOMA =(y[0] + kappa * S_CX_SOMA * g2_SOMA) /
	(1 + kappa*S_CX_SOMA * g1_SOMA);
  
 CX_DEND::calc(x, y, f, 
               CX_SOMA::Ki, CX_SOMA::Nai, CX_SOMA::Ko, Ko2d, CX_SOMA::Nao, Nao2d,
               TotalInhibit);

 f[0] = C_inv * (iDEND  +  1.0/(kappa*S_CX_DEND) * (v_SOMA - y[0]));

}


/////////////////////////////////////////////////////////////////////////////
/////// SYNAPSES

//------second order kiner model (including G-proteins) for GABA-B synapce----
class Gaba_B {
 static double E_GABA, Cmax, Deadtime, Prethresh;
 static double Kd, n;
 double Gn, q;
 public:
 double C, lastrelease;
 double I, r0, g0, Gn1, Cdur, K1, K2, K1K2, K3, K4;
 Gaba_B() {
  Cdur = 0.3;
  K1 = 0.52;
  K2 = 0.0013;
  K3 = 0.098;
  K4 = 0.033;
  lastrelease = -10000000;
  C = 0, r0 = 0, g0 = 0;
 }
 void calc(double r, double g, double &fr, double &fg,
           double g_GABA_B, double x, double y_pre, double y_post);
};
double Gaba_B::E_GABA = -95, Gaba_B::Cmax = 0.5, Gaba_B::Deadtime = 1;
double Gaba_B::Prethresh = 0;
double Gaba_B::Kd = 100, Gaba_B::n = 4;



void Gaba_B::calc(double r, double g, double &fr, double &fg,
                  double g_GABA_B, double x, double y_pre, double y_post) {
 Gn = pow(g,n);
 Gn1 = Gn/(Gn + Kd);
 I = g_GABA_B * Gn1 * (y_pre - E_GABA);

 q = ((x - lastrelease) - Cdur);
 if (q > Deadtime) {
  if (y_post > Prethresh) {
   C = Cmax;
   lastrelease = x; }
 } else if (q < 0) {
 } else if (C == Cmax) {   C = 0;  }
 fr = K1 * C * (1 - r) - r * K2;
 fg = K3 * r - K4 * g;
}

class GB: public Gaba_B {
 public:
 GB() :Gaba_B() { }
 void init(double *y){
  lastrelease = -10000000;
  C = 0;
  y[0] = 0;
  y[1] = 0;
 }
 void calc(double g_GABA_B, double x, double *y, double *f,
           double y_pre, double y_post){
  Gaba_B::calc(y[0], y[1], f[0], f[1], g_GABA_B, x, y_pre, y_post);
 }
};



// Excitatory synaptic transmission
// AMPAergic synapse (with depression)

class AMPA_D2 {
 static double E_AMPA;
 static double Cdur, Cmax, Deadtime, Prethresh, Cdel;
 static double Alpha,Beta;
 int spike;
 double R, Rinf, q;
 double lastrelease;
 double prev[Mcx][Mcx1];
 double exptable(double z)
 {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return( 0 );
 }
 public:
 double I, g1,Rout;
 double E[Mcx][Mcx1];
 AMPA_D2() {
  R = 0, Rout =0;
  lastrelease = -10000;
  Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
  for(int j=Mcx1-1; j>=0; --j)
   for(int i=Mcx-1; i>=0; --i){
    prev[i][j] = 0;
   }
 }
 void iii(unsigned int seek) {srand(seek);}
 void calc(double, double, double, double, double, double, int, int);
};
double AMPA_D2::E_AMPA = 0, AMPA_D2::Cdur = 0, AMPA_D2::Deadtime = 2;
double AMPA_D2::Cdel = 0, AMPA_D2::Cmax = 0.5;
double AMPA_D2::Prethresh = 0, AMPA_D2::Beta = 0.18,AMPA_D2::Alpha = 0.94;
void AMPA_D2::calc(double g_AMPA, double g_AMPAmin, double x, double y_post,
                   double y_pre, double my_E, int i, int j) {

 q = ((x - prev[i][j]) - Cdur);

 if(q > Deadtime) {
  if(y_pre > Prethresh) {
   g1 = g_AMPA;
   lastrelease = x;
   spike = 1;
   prev[i][j] = x;
  }
 }

 if (spike == 1){
  Rout = Rout+(g1*my_E);
  R = Rout;
  spike = 0;
 } else {
  Rout = R * exptable (-Beta * (x - (lastrelease)));
 }

 I =  Rout * y_post;
}

//------------first order kiner model for NMDA synapce WITH depression------
class NMDA_D1 {
 static double E_NMDA;
 static double Cdur, Cmax, Deadtime, Prethresh, Cdel;
 static double Alpha, Beta;
 int s;
 int spike;
 double prev[Mcx][Mcx1];
 double R, C, Rold, R1;
 double lastrelease, lastspike;
 double q, Rinf, Rtau, fn;
 double exptable(double z)
 {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return( 0 );
 }
 public:
 double I, Rout, g;
 int myflag;
 NMDA_D1() {
  R = 0, C = 0, Rold = 0, R1 = 0, Rout = 0;
  lastrelease = -100;
  lastspike = -100;
  s = 1;
  Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
  Rtau = 1 / ((Alpha * Cmax) + Beta);
  for(int j=Mcx1-1; j>=0; --j)
   for(int i=Mcx-1; i>=0; --i){
    prev[i][j] = 0;
   }
 }
 void iii(unsigned int seek) {srand(seek);}
 void calc(double g_NMDA, double x, double y_pre, double y_post, int i, int j);
};
double NMDA_D1::E_NMDA = 0, NMDA_D1::Cdur = 0, NMDA_D1::Cmax = 0.3;
double NMDA_D1::Deadtime = 2;
double NMDA_D1::Prethresh = 0, NMDA_D1::Alpha = 1, NMDA_D1::Beta = 0.0067;

void NMDA_D1::calc(double g_NMDA, double x, double y_pre, double y_post, int i, int j) {

 q = ((x - prev[i][j]) - Cdur);

 if(q > Deadtime) {
  if(y_post > Prethresh) {

   prev[i][j] = x;
   lastrelease = x;
   Rold = Rout;}
 }

 if ((x - lastrelease) < Cmax) {
  Rout = Rold + Rinf + (- Rinf) * exptable (-(x - lastrelease) / (Rtau));
  R1 = Rout;
 } else {
  Rout = R1 * exptable (-Beta * (x - (lastrelease+Cmax)));
 }

 g = g_NMDA * Rout;
 I = g_NMDA * Rout;
}




//---first order kiner model for GABA-A synapce with DEPRESSION & spont IPSPs--
class Gaba_A_D2 {
 static double Cdur, Cmax, Deadtime, Prethresh;
 static double Alpha, Beta;
 double R;
 double lastrelease, lastrelease1, lastrandom, newrelease, Use, Tr, UseSlow, TrSlow;
 double q, q1, Rinf, Rtau;
 int spike, id;
 double Tau, Period, S, SS, factor;
 double prev[Min][Min1];
 double exptable(double z)
 {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return( 0 );
 }
 public:
 double I, E_GABA,g, Rout;
 Gaba_A_D2() {
  E_GABA = -70;
  R = 0, Rout = 0;
  lastrelease = -10000;

  Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
  Rtau = 1 / ((Alpha * Cmax) + Beta);
  g = 0;
  for(int j=Min1-1; j >=0; --j)
   for(int i=Min-1; i >=0; --i){
    prev[i][j] =0;
   }


 }
 void iii(unsigned int seek) {srand(seek);}
 void calc(double g_GABA_A, double x, double y_pre, double y_post, int i, int j, double R_Pot);
};
double Gaba_A_D2::Cdur = 0, Gaba_A_D2::Cmax = 0.5, Gaba_A_D2::Deadtime = 2;
double Gaba_A_D2::Prethresh = 0, Gaba_A_D2::Alpha = 10, Gaba_A_D2::Beta =0.25;
void Gaba_A_D2::calc(double g_GABA_A, double x, double y_pre,
                     double y_post, int i, int j, double R_Pot) {

 q = ((x - prev[i][j]) - Cdur);

 if ((q > Deadtime)) {
  if (y_post > Prethresh) {
   prev[i][j] = x;
   lastrelease = x;
   spike = 1;
  }
 }

 if (spike == 1){
  Rout = Rout + g_GABA_A;
  R = Rout;
  spike = 0;
 }

 else {
  Rout = R * exptable (-Beta * (x - (lastrelease)));
 }

 g = Rout;
 I = Rout * (y_pre - R_Pot);
}



//-----first order kiner model for AMPA synapce used for external stimulation----
class Extern_ampa1 {
 static double Cdur, Cmax, Deadtime, Prethresh;
 double R, C, R0, R1;
 double lastrelease;
 int myflag, spike;
 double q, Rinf, Rtau;
 double wom, RRR;
 drand48_data *rand_buffer;  

 double exptable(double z)
 {
  if((z > -10) && (z < 10)) return( exp(z) );
  else return( 0 );
 }
 public:
 double g, Alpha, Beta, TR, Rout, w;
 Extern_ampa1() {
  Alpha = 0.94;
  Beta = 0.18;
  R = 0, C = 0, R0 = 0, R1 = 0, Rout =0;
  lastrelease = -100;
  myflag = 0;
  spike = 0;
  Rinf = Cmax*Alpha / (Cmax*Alpha + Beta);
  Rtau = 1 / ((Alpha * Cmax) + Beta);
  TR = 25;

  rand_buffer = new drand48_data;
  static int seed = time2(NULL);
  seed = seed + 100;
  srand48_r(seed,rand_buffer);
  w=0.2;
 }
 void iii(unsigned int seek) {/*srand(seek);*/ srand48_r(seek,rand_buffer);}
 void calc(double g_Extern_ampa, double x);
};


double Extern_ampa1::Cdur = 0.3, Extern_ampa1::Cmax = 0.5, Extern_ampa1::Deadtime = 1;
double Extern_ampa1::Prethresh = 0;
void Extern_ampa1::calc(double g_Extern_ampa, double x) {

 q = ((x - lastrelease) - Cdur);
 if (q > Deadtime) {
  if ((x - lastrelease) > TR) {
   C = Cmax;
   R0 = R;
   lastrelease = x;
   drand48_r(rand_buffer,&RRR);
   if(RRR < 0.000001) RRR = 0.000001;
   TR = -(log(RRR))/(w);
   spike =1;
  }
 }
 if (spike == 1){
  Rout = Rout + g_Extern_ampa;
  R1 = Rout;
  spike = 0;
 }
   
 else {
  Rout = R1 * exptable (-Beta * (x - (lastrelease)));
 }
 g =  Rout;
}

