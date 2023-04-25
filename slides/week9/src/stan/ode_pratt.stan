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

vector get_flux_pratt(real t, vector x, vector p){
  real diet_tg = get_diet(t, [0, 300, 600]', 0.00185859360751066, 180);
  real diet_glucose = get_diet(t, [0, 300, 600]', 0.139081199343437, 180);
  real m_amp = p[55] / (p[56] * p[57] * x[18] + p[58] * p[59] * x[16] * x[1]);
  return [
    p[3] + p[4] *tanh((x[2]- p[5])/(0.8*p[6])),     // vP_insulin_in 
    mass_action_one_substrate(p[7], x[1]),          // vP_insulin_out
    diet_tg,                                        // vP_chylo_in
    mass_action_one_substrate(p[8], x[3]),          // vP_TG_to_L_FFA
    p[9] * x[11] / (p[10] + x[11]),                 // vL_TG_stor_to_P_TG
    mass_action_one_substrate(p[11], x[12]),        // vL_TG_secr_to_P_TG
    (1 + p[12] * x[1]) * p[13] * x[2],              // vP_glucose_to_M_glucose
    mass_action_two_substrates(p[14], x[3], p[1]),  // vP_TG_to_M_FFA
    p[15] * (1 + p[16] * x[1]) * x[2],              // vP_glucose_to_A_glucose; p[16] i.e. kga is hardcoded to zero??
    p[17],                                          // vP_glucose_out
    mass_action_one_substrate(p[18], x[4]),         // vP_NEFA_to_M_FFA
    mass_action_one_substrate(p[19], x[5]),         // vP_chylo_to_M_FFA
    mass_action_one_substrate(p[20], x[4]),         // vP_NEFA_to_A_FFA
    mass_action_two_substrates(p[21], x[3], p[2]),  // vP_TG_to_A_FFA
    mass_action_one_substrate(p[22], x[5]),         // vP_chylo_to_L_FFA
    mass_action_one_substrate(p[23], x[2]),         // vP_glucose_to_L_glucose
    mass_action_one_substrate(p[24], x[4]),         // vP_NEFA_to_L_FFA
    (1 - p[25]) * p[26] * (1 + p[27] * x[1]) * x[5] * p[2], //vP_chylo_to_A_FFA
    p[25] * p[26] * (1 + p[27] * x[1]) * x[5] * p[2],  // vP_chylo_to_P_NEFA
    diet_glucose,                                   // vL_glucose_in
    p[28] * x[6] / (p[29] + x[6]) + p[30] * x[6] / (p[31] + x[6]) * (1 / (1 + p[32] * x[8])), // vL_glucose_to_L_G6P
    p[33] * x[1] * x[8] * 0.5 * (1 + tanh((p[34] - x[7])/ p[35])), // vL_G6P_to_L_glycogen
    mass_action_two_substrates(p[36], x[8], x[1]),  // vL_G6P_to_L_pyruvate
    p[37],                                          // vL_pyruvate_in
    mass_action_two_substrates(p[38] * p[39], x[9], x[1]), // vL_pyruvate_to_L_FFA
    mm_one_substrate(p[40], p[41], x[12]),          // vL_FFA_to_L_TG_secr
    mm_one_substrate(p[42], p[43], x[12]),          // vL_FFA_to_L_TG_stor
    mass_action_one_substrate(p[44], x[6]),         // vL_glucose_to_P_glucose
    p[45] / (1 + p[46] * x[1]) * x[7] / (x[7] + p[47]), // vL_glycogen_to_L_G6P
    p[48] * x[12] / (1 + p[49] * x[1]),             // vL_FFA_out
    mm_one_substrate(p[50], p[51], x[11]),          // vL_TG_stor_to_L_FFA
    mass_action_one_substrate(p[52], x[8]),         // vL_G6P_to_L_glucose
    p[53] * x[9] / (1 + p[54] * x[1]),              // vL_pyruvate_to_L_G6P
    p[60],                                          // vM_TG_to_M_FFA
    mass_action_two_substrates(p[61], x[1], x[18]), // vM_FFA_to_M_TG
    mm_one_substrate(p[62], p[63], x[13]) * inhibition(p[32], x[15]), // vM_glucose_to_M_G6P
    mm_one_substrate(p[64], p[65], x[14]) * inhibition(p[66], x[1]), // vM_glycogen_to_M_G6P
    p[67] * x[1] * x[15] * 0.5 * (1 + tanh((p[68] - x[14]) / p[69])), // vM_G6P_to_M_glycogen 
    mass_action_two_substrates(p[70], x[15], x[1]), // vM_G6P_to_M_pyruvate
    (1 + p[12] * x[1]) * p[71] * x[13],             // vM_glucose_to_P_glucose        
    p[58] * x[16] * x[1] * m_amp,                   // vM_pyruvate_out                
    mass_action_two_substrates(p[56], x[18], m_amp),// vM_FFA_out                     
    mass_action_one_substrate(p[72], x[16]),        // vM_pyruvate_to_L_pyruvate      
    p[73] * x[1] * x[19] * x[22],                   // vA_FFA_and_A_glucose_to_A_TG   
    mass_action_one_substrate(p[74], x[20]),        // vA_glycerol_to_L_G6P           
    p[75] / (1 + p[76] * x[1]^2),                   // vA_TG_to_P_NEFA_and_A_glycerol 
    p[15] * (1 + p[16] * x[1]) * x[19]              // vA_glucose_to_P_glucose
  ]';
}
vector get_dxdt_pratt(real t, vector x, vector p, matrix S, matrix S_cpt, vector vol){
  vector[cols(S)] flux = get_flux_pratt(t, x, p);
  vector[rows(x)] a = (S * flux);
  vector[rows(x)] b = (S_cpt * vol);
  return a ./ b;
}
