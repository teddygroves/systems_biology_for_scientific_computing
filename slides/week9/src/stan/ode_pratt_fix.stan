/* Functions for calculating ODEs in the whole body model. */

real mm_one_substrate(real vmax, real km, real sub){
  return vmax * sub / (km + sub);
}

real mass_action_one_substrate(real k, real x){
  return k * x;
}

real mass_action_two_substrates(real k, real x1, real x2){
  return k * x1 * x2;
}

real inhibition(real ki, real x){
  return 1 / (1 + ki * x);
}

real get_diet(real t, vector breaks, real k, real d){
  real out = 0;
  for (b in breaks){
    if (t > b){
      out += k * (t - b) * exp(-(t - b) / d);
    }
  }
  return out;
}

vector get_flux_pratt_fix(real t, vector x, vector p){
  // parameters
  real M_LPL = p[1];
  real A_LPL = p[2];
  real k11 = p[3];
  real k22 = p[4];
  real nu = p[5];
  real cc = p[6];
  real kd = p[7];
  real kr = p[8];
  real nu9 = p[9];
  real k9 = p[10];
  real k9a = p[11];
  real kgi = p[12];
  real kgm = p[13];
  real kt = p[14];
  real dba = p[15];
  real kga = p[16];
  real mu1 = p[17];
  real kbm = p[18];
  real kcm = p[19];
  real kna = p[20];
  real kba = p[21];
  real kcl = p[22];
  real kgl2 = p[23];
  real kbl = p[24];
  real klp = p[25];
  real ka = p[26];
  real kai = p[27];
  real nuLG = p[28];
  real kLG = p[29];
  real nuLH = p[30];
  real kLH = p[31];
  real krep = p[32];
  real kyl = p[33];
  real lmax = p[34];
  real cl0 = p[35];
  real kp = p[36];
  real mub = p[37];
  real PtF = p[38];
  real kal = p[39];
  real v6 = p[40];
  real k6 = p[41];
  real v8 = p[42];
  real k8 = p[43];
  real kgl = p[44];
  real betaL =p[45];
  real kdl = p[46];
  real yl0 = p[47];
  real k7 = p[48];
  real k5 = p[49];
  real v10 = p[50];
  real k10 = p[51];
  real k6L = p[52];
  real beta6 = p[53];
  real kp6 = p[54];
  real mu_amp = p[55];
  real mu4 = p[56];
  real FED = p[57];
  real mu3 = p[58];
  real CED = p[59];
  real mue = p[60];
  real mus = p[61];
  real nuMH = p[62];
  real kMH = p[63];
  real betam = p[64];
  real ym0 = p[65];
  real kdy = p[66];
  real kym = p[67];
  real mmax = p[68];
  real cm0 = p[69];
  real k6p = p[70];
  real kgm2 = p[71];
  real kpp = p[72];
  real kaa = p[73];
  real kgp = p[74];
  real betaf = p[75];
  real kft = p[76];
  // States:
  real P_insulin = x[1];
  real P_glucose = x[2];
  real P_TG = x[3];
  real P_NEFA = x[4];
  real P_chylo = x[5];
  real L_glucose = x[6];
  real L_glycogen = x[7];
  real L_G6P = x[8];
  real L_pyruvate = x[9];
  real L_TG_secr = x[10];
  real L_TG_stor = x[11];
  real L_FFA = x[12];
  real M_glucose = x[13];
  real M_glycogen = x[14];
  real M_G6P = x[15];
  real M_pyruvate = x[16];
  real M_TG = x[17];
  real M_FFA = x[18];
  real A_glucose = x[19];
  real A_glycerol = x[20];
  real A_TG = x[21];
  real A_FFA = x[22];
  // calculated
  real diet_TG = get_diet(t, [0, 300, 600]', 0.00185859360751066, 180);
  real diet_glucose = get_diet(t, [0, 300, 600]', 0.139081199343437, 45);
  real M_AMP = mu_amp / (mu4 * FED * M_FFA + mu3 * CED * M_pyruvate * P_insulin);
  return [
    k11+k22*tanh((P_glucose-nu)/(0.8*cc)),
    kd*P_insulin,
    diet_TG,
    kr*P_TG,
    nu9*L_TG_stor/(k9+L_TG_stor),
    k9a*L_TG_secr,
    (1+kgi*P_insulin)*kgm*P_glucose,
    kt*P_TG*M_LPL,
    dba*(1+kga*P_insulin)*P_glucose,
    mu1,
    kbm*P_NEFA,
    kcm*P_chylo,
    kna*P_NEFA,
    kba*P_TG*A_LPL,
    kcl*P_chylo,
    kgl2*P_glucose,
    kbl*P_NEFA,
    (1.0-klp)*ka*(1.0+kai*P_insulin)*P_chylo*A_LPL,
    klp*ka*(1.0+kai*P_insulin)*P_chylo*A_LPL,
    diet_glucose,
    nuLG*L_glucose/(kLG+L_glucose)+nuLH*L_glucose/(kLH+L_glucose)*(1.0/(1.0+krep*L_G6P)),
    kyl*P_insulin*L_G6P*(1.0/2.0)*(1.0+tanh((lmax-L_glycogen)/cl0)),
    kp*L_G6P*P_insulin,
    mub,
    PtF*kal*L_pyruvate*P_insulin,
    v6*L_FFA/(k6+L_FFA),
    v8*L_FFA/(k8+L_FFA),
    kgl*L_glucose,
    betaL/(1+kdl*P_insulin)*L_glycogen/(L_glycogen+yl0),
    k7*L_FFA/(1+k5*P_insulin),
    v10*L_TG_stor/(k10+L_TG_stor),
    k6L*L_G6P,
    beta6*L_pyruvate/(1+kp6*P_insulin),
    mue,
    mus*P_insulin*M_FFA,
    nuMH*(M_glucose/(kMH+M_glucose))*(1.0/(1.0+krep*M_G6P)),
    betam*(M_glycogen/(M_glycogen+ym0))/(1.0+kdy*P_insulin),
    kym*P_insulin*M_G6P*(1.0/2.0)*(1.0+tanh((mmax-M_glycogen)/cm0)),
    k6p*M_G6P*P_insulin,
    (1.0+kgi*P_insulin)*kgm2*M_glucose,
    mu3*M_pyruvate*P_insulin*M_AMP,
    mu4*M_FFA*M_AMP,
    kpp*M_pyruvate,
    kaa*P_insulin*A_glucose*A_FFA,
    kgp*A_glycerol,
    betaf/(1.0+kft*pow(P_insulin,2.0)),
    dba*(1.0+kga*P_insulin)*A_glucose
  ]';
}
vector get_dxdt_pratt_fix(real t, vector x, vector p, matrix S){
  return S * get_flux_pratt_fix(t, x, p);
}

