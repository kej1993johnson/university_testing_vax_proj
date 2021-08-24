
#2 means without the isolation abandonment rate

# Functions for cost-effectiveness of testing strategy in a university
solve_transmission_eqns<- function(par_table, vacc_inf_correlation = 'independent', test_type = "antigen"){
  par_names<-colnames(par_table)
  for (i in 1:length(par_names)){
    assign(par_names[i], as.double(par_table[i]))
  }
  
  if (vacc_inf_correlation =='independent'){
  # Assume vaccination is independent of being previously infected
  # P(R|V) = P(R|!V) = P(R)
  init_ru = init_rec- init_rec*init_vacc
  init_rv = init_rec*init_vacc
  init_sv = (1-init_rec)*init_vacc 
  }
  else if (vacc_inf_correlation == 'correlated'){
    # P(R|!V) = x*P(R)
    x = odds_inf_not_vac
    prob_rec_given_vac = init_rec/(init_vacc + x*(1-init_vacc))
    prob_rec_given_not_vac = x*prob_rec_given_vac
    init_ru = prob_rec_given_not_vac*(1-init_vacc) #P(R|!V)*P(!V)
    init_su = (1-prob_rec_given_not_vac)*(1-init_vacc) #P(!R|!V)*P(!V)
    init_rv = prob_rec_given_vac*init_vacc # P(R|V)*P(V) 
    init_sv = (1-prob_rec_given_vac)*init_vacc # P(!R|V)*P(V)
    test<-init_ru+init_rv+init_sv+init_su
    p_rec = init_ru + init_rv
    p_vacc = init_rv + init_sv
  }
  #VE = 1-odds_case_vacc/odds_pop_vacc
  odds_pop_vacc = init_vacc/(1-init_vacc)
  odds_case_vacc = odds_pop_vacc*sigma_v
  prob_cases_vacc = odds_case_vacc/(1+odds_case_vacc)
  
  psym_v = sym_red/sigma_v
  avg_time_asym = psym*t_pre_sym + (1-psym)*t_inf
  avg_time_asym_v = psym_v*t_pre_sym + (1-psym_v)*t_inf
  avg_time_sym_undetected = (1-is_sym)*(t_sym_state) + is_sym*t_iso
  delta_au = 1/avg_time_asym
  delta_av = 1/avg_time_asym_v
  delta_s = 1/t_sym_state
  delta_s = 1/avg_time_sym_undetected
  delta_q = 1/t_quarantine
  beta_a = R0*1/t_inf #transmission rate of asymptomatic individuals
  beta_s = R0*1/t_inf # transmission rate of symptomatic individuals 
  intro_rate = init_prev/7 #  p = new infections per day per capita in surrounding community. 
  tdays<-seq(from = 0, to = duration_school_year, by = 1)
  dates<-seq(from = as.Date("2021-08-20"), by = "day", length.out = duration_school_year+1)
  # define initial conditions of compartments
  Squ =0
  Equ = 0
  Iaqu = 0
  Isqu = 0
  Rqu = 0
  Su = N*(1-init_ru - init_rv - init_sv - init_prev)# (unvaccinated susceptibles = 1-unvaccinated recovered - vaccinated susceptible - currently infected)
  Eu = 0
  Iau = 0.5*N*init_prev*(1-sigma_v)*(1-prob_cases_vacc) # % pre/asymptomatic * percent infected* percent unvaccinated
  Isu =0.5*N*init_prev*(1-prob_cases_vacc) # % symptomatic * percent infeted * percent unvaccinated
  Ru = N*(init_ru)
  Sqv = 0
  Eqv = 0
  Iaqv = 0
  Isqv = 0
  Rqv = 0
  Sv = N*(init_sv)
  Ev = 0
  Iav = 0.5*N*init_prev*prob_cases_vacc
  Isv = 0.5*N*init_prev*prob_cases_vacc
  Rv = N*(init_rv)
  Ntest = Su+Eu+Iau+Isu+Ru+Sv+Ev+Iav+Isv+Rv # should equal N
  I0 = Iau + Isu + Iaqu + Isqu
  Q0 = Squ + Equ + Iaqu + Isqu + Sqv + Eqv + Iaqv + Isqv 
  Nu = Su+Eu+Iau+Isu +Ru
  Nv = Sv + Ev + Iav + Isv + Rv
  RTs = 0
  PCR = 0
  cumI = 0
  cumsI = 0
  
  # Account for differences in timing of the different test types
  if (test_type == 'antigen'){
    inf_rem = 1
  }
  if (test_type == 'PCR'){
    pct_inf_removed = get_pct_infectiousness_removed(V0=3, Vf=6, delta_add=0, t0topeak=2, Vpeak = 9,
    tpeaktof=7, LOD=5, LOI=6, dt=0.01, t_test = 1, test_sensitivity = 0.95)
    inf_rem<-pct_inf_removed$eff_percent_removed
  }
  
  y <- c(Squ = Squ, Equ=Equ, Iaqu=Iaqu, Isqu=Isqu, Rqu = Rqu, Su=Su, Eu=Eu, Iau=Iau, Isu=Isu, Ru=Ru, Sqv=Sqv, 
         Eqv=Eqv, Iaqv=Iaqv, Isqv=Isqv, Rqv = Rqv, Sv=Sv, Ev=Ev, Iav=Iav, Isv=Isv, Rv=Rv, I=I0,Isymdet = Isqu+Isqv, Idet =Isqu+Isqv+Iaqu+Iaqv,
         Q = Q0, Nu = Nu, Nv = Nv, RTs =0, PCR = 0, cumI =0, cumIvax=0, cumIunvax=0, cumsI=0, FPs=0, TPs=0)
  params<-c(delta_au=delta_au, delta_av = delta_av, delta_s=delta_s, delta_q=delta_q, beta_a=beta_a, beta_s=beta_s, gamma=gamma,
            psym = psym, e_v=e_v, sigma_v=sigma_v, Sp=Sp, Se=Se, 
            w_v =w_v, w_u=w_u, f_v=f_v, f_u=f_u, is=is, k=k, N=N, intro_rate = intro_rate,
            init_vacc= init_vacc, inf_rem = inf_rem)
  model<-function(tdays, y, params){
    with(as.list(c(y,params)), {
      X = Su+Eu+Iau+Isu+Ru+Sv+Ev+Iav+Isv+Rv # active population
      beta_star = beta_a*(Iau + e_v*Iav) + beta_s*(Isu + e_v*Isv) # force of infection
      # quarantined unvaccinated
      dSqu<- Su*w_u*(1-Sp)*is*f_u - Squ*k
      dEqu<-Eu*w_u*(1-Sp)*is*f_u - Equ*k
      dIaqu<-Iau*w_u*Se*is*f_u*inf_rem - delta_q*Iaqu 
      dIsqu<-Isu*w_u*Se*is*f_u*inf_rem + is_sym*delta_s*Isu  - delta_q*Isqu 
      dRqu<- Ru*w_u*(1-Sp)*is*f_u - Rqu*k
      # unvaccinated
      dSu<- - beta_star*Su/X - Su*w_u*(1-Sp)*is*f_u + Squ*(k) - (1-prob_cases_vacc)*intros_per_week/7
      dEu<- beta_star*Su/X- gamma*Eu - Eu*w_u*(1-Sp)*is*f_u + (1-prob_cases_vacc)*intros_per_week/7
      dIau<- gamma*Eu - delta_au*Iau  - Iau*w_u*Se*is*f_u*inf_rem + Equ*k 
      dIsu<- psym*delta_au*Iau - delta_s*Isu - Isu*w_u*Se*is*f_u*inf_rem 
      dRu<- (1-is_sym)*delta_s*Isu + (1-psym)*delta_au*Iau - Ru*w_u*(1-Sp)*is*f_u + Rqu*k  + delta_q*Isqu +delta_q*Iaqu
      # quarantined vaccinated
      dSqv<-Sv*w_v*(1-Sp)*is*f_v - Sqv*k
      dEqv<- Ev*w_v*(1-Sp)*is*f_v - Eqv*k
      dIaqv<- Iav*w_v*Se*is*f_v*inf_rem - delta_q*Iaqv 
      dIsqv<-Isv*w_v*Se*is*f_v*inf_rem +is_sym*delta_s*Isv  - delta_q*Isqv
      dRqv<-Rv*w_v*(1-Sp)*is*f_v - Rqv*k
      # vaccinated
      dSv<- -sigma_v*beta_star*Sv/X - Sv*w_v*(1-Sp)*is*f_v + Sqv*(k)- prob_cases_vacc*intros_per_week/7
      dEv<- sigma_v*beta_star*Sv/X - gamma*Ev - Ev*w_v*(1-Sp)*is*f_v + prob_cases_vacc*intros_per_week/7
      dIav<- gamma*Ev - delta_av*Iav  - Iav*w_v*Se*is*f_v*inf_rem + Eqv*(k) 
      dIsv<-psym_v*delta_av*Iav - delta_s*Isv - Isv*w_v*Se*is*f_v*inf_rem 
      dRv<- delta_s*(1-is_sym)*Isv + (1-psym_v)*delta_av*Iav - Rv*w_v*(1-Sp)*is*f_v + Rqv*k + delta_q*Isqv +delta_q*Iaqv
      dI<- dIau + dIsu + dIaqu + dIsqu + dIav + dIsv + dIaqv + dIsqv
      dIsymdet<-dIsqu + dIsqv
      dIdet<-dIsqu + dIsqv + dIaqu + dIaqv
      dQ<-dSqu + dEqu + dIaqu + dIsqu + dRqu + dSqv + dEqv + dIaqv + dIsqv + dRqv
      dNu<-dSu+dEu+dIau+dIsu +dRu + dSqu+dEqu+dIaqu+dIsqu +dRqu
      dNv<- dSv + dEv + dIav + dIsv + dRv + dSqv + dEqv + dIaqv + dIsqv + dRqv
      dRTs<-Nu*f_u*w_u + Nv*f_v*w_v
      dPCR<- Su*w_u*(1-Sp)*is*f_u + Eu*w_u*(1-Sp)*is*f_u + Iau*w_u*Se*is*f_u + Isu*w_u*Se*is*f_u + is*Isu+
        Ru*w_u*(1-Sp)*is*f_u + Sv*w_v*(1-Sp)*is*f_v + Ev*w_v*(1-Sp)*is*f_v + Iav*w_v*Se*is*f_v + Isv*w_v*Se*is*f_v+
        is*Isv + Rv*w_v*(1-Sp)*is*f_v# all new quarantines
      dcumI<-beta_star*Su/X + sigma_v*beta_star*Sv/X # all new infections 
      dcumIvax<-sigma_v*beta_star*Sv/X
      dcumIunvax<-beta_star*Su/X
      dcumsI<- Isu*w_u*Se*is*f_u + is_sym*Isu + Isv*w_v*Se*is*f_v +is_sym*Isv# all new symptomatic detected infections
      dFP<-Su*w_u*(1-Sp)*is*f_u + Eu*w_u*(1-Sp)*is*f_u + Ru*w_u*(1-Sp)*is*f_u + Sv*w_v*(1-Sp)*is*f_v +
        Ev*w_v*(1-Sp)*is*f_v + Rv*w_v*(1-Sp)*is*f_v  #-> new false positives?
      dTP<- Iau*w_u*Se*f_u + Isu*w_u*Se*f_u + is_sym*delta_s*Isu + Iav*w_v*Se*f_v + Isv*w_v*Se*f_v +
        is_sym*delta_s*Isv#-> new true positives (this includes some potential for double counting of the people who test positive 
      #but don't quarantine, since they stay in the active compartment)
      #dTP<-dIaqu + dIsqu + dIaqv + dIsqv
      #dFP<-dSqu + dEqu + dRqu + dSqv + dEqv + dRqv 
      
      
      list(c(dSqu, dEqu, dIaqu, dIsqu, dRqu, dSu, dEu, dIau, dIsu, dRu, dSqv, dEqv, dIaqv, dIsqv, dRqv, dSv,
             dEv, dIav, dIsv, dRv, dI,dIsymdet, dIdet, dQ, dNu, dNv, dRTs, dPCR, dcumI, dcumIvax, dcumIunvax, dcumsI, dFP, dTP))
    })
  }
  
  out<-ode(y, tdays,model, params)
  #plot(out)
  Ntot<-rowSums(out[,2:21])
  #plot(tdays, Ntot)
  out_df<-data.frame(out)
  # add some columns
  init_imm<-1-(1-init_vacc)*(1-init_rec)
  pct_vacc<-rep(init_vacc, length(tdays))
  pct_imm<-rep(init_imm, length(tdays))
  testing_freq<-rep(f_u, length(tdays))
  prevalence<-rep(init_prev, length(tdays))
  
  out_df<-cbind(dates, out_df,  pct_vacc, testing_freq, prevalence)
  
  return(out_df)
}

testing_model<-function(par_table, cost_table, risk_tolerance, vacc_inf_correlation = 'independent', test_type = 'antigen'){
  out_df<-solve_transmission_eqns(par_table, vacc_inf_correlation, test_type)

  # Gets the outputs of interest for the cost-effectiveness model
  # input cost parameters
  
  cost_names<-colnames(cost_table)
  for (i in 1:length(cost_names)){
    assign(cost_names[i], as.double(cost_table[i]))
  }
  
  if (risk_tolerance =='CDC yellow' ){
    close_thres = 10/100000
  }
  
  if (risk_tolerance == 'CDC orange'){
    close_thres = 50/100000
  }
  
  if (risk_tolerance == 'CDC red'){
    close_thres = 100/100000
  }
  
  if (risk_tolerance == 'UT Jan 2021'){
    close_thres = (20 + 0.7*32)/30000 # week ending on 1/16 --> forced school remote in Jan
  }
  
  thres_yellow = 10/100000
  thres_orange = 50/100000
  thres_red = 100/100000
  thres_UTJan = (20 + 0.7*32)/30000
  
  
  # outputs
  N<-par_table$N
  n_weeks<-par_table$duration_school_year/7
  n_inf = out_df$cumI[nrow(out_df)] # total number infected at end of the year
  n_inf_vax = out_df$cumIvax[nrow(out_df)]
  n_inf_unvax = out_df$cumIunvax[nrow(out_df)]
  max_symdet<-max(out_df$Isymdet)
  n_sym_det<-out_df$cumsI[nrow(out_df)]
  pct_inf<-n_inf/N
  pct_detected = out_df$TPs[nrow(out_df)]/out_df$cumI[nrow(out_df)]
  testing_freq<-par_table$f_u
  n_pos<-out_df$TPs[nrow(out_df)]
  days_of_online = as.numeric(sum((out_df$Isymdet)/N>close_thres))
  days_above_yellow = as.numeric(sum((out_df$Isymdet)/N>thres_yellow))
  days_above_orange = as.numeric(sum((out_df$Isymdet)/N>thres_orange))
  days_above_red = as.numeric(sum((out_df$Isymdet)/N>thres_red))
  days_above_UTJan = as.numeric(sum((out_df$Isymdet)/N>thres_UTJan))
  cross_yellow <- days_above_yellow>0
  cross_orange <- days_above_orange>0
  cross_red <-days_above_red>0
  cross_UTJan<-days_above_UTJan>0
  DLL= sum(out_df$Q) # integral of Q (days of learning lost from quarantine only)
  #DLL_tot = sum(out_df$Q) + sum(N-out_df$Q[out_df$Isymdet/N>close_thres]) # plus days online days*
  # the population not in quarantine 
  cost_of_online = days_of_online*online_price
  cost_DLL = DLL*DLL_price # resources spent by students in response to COVID
  cost_RT<- out_df$RTs[nrow(out_df)]*RT_price
  cost_PCR <-out_df$PCR[nrow(out_df)]*PCR_price
  cost_testing<-cost_PCR+cost_RT
  cost_isofac<-out_df$TPs[nrow(out_df)]*pct_pos_isofac*isofac_price
  cost_sequencing<-out_df$TPs[nrow(out_df)]*sequencing_price
  # calculate the number of kiosks needed:
  pct_per_week_tested<-par_table$f_u*7*par_table$w_u + par_table$f_v*7*par_table$w_v
  pct_per_week_tested_uv<-par_table$f_u*7*par_table$w_u
  pct_per_week_tested_v<-par_table$f_v*7*par_table$w_v
  n_tests_per_week<- (1-par_table$init_vacc)*par_table$N*pct_per_week_tested_uv + 
    (par_table$init_vacc)*par_table$N*pct_per_week_tested_v +out_df$PCR[nrow(out_df)]/n_weeks
  n_kiosks<-n_tests_per_week/n_rapid_tests_per_kiosk_per_week
  cost_kiosks<-n_kiosks*price_per_kiosk
 
  cost_contact_tracing<-out_df$TPs[nrow(out_df)]*contact_tracing_price

  n_PCR<-out_df$PCR[nrow(out_df)]
  n_RT<-out_df$RTs[nrow(out_df)]
  cost_to_UT<- cost_testing + cost_isofac + cost_contact_tracing + cost_sequencing + cost_kiosks + cost_of_online
  #testing, contact-tracing, isolation facilities, cost to maintain a PCR lab/sequencing, 
    #cost for staff/faculty/administrator time
  
  
  # to find things averted, need to run with no testing at all
  par_table$f_u<-0
  par_table$f_v<-0
  out_nt_df<-solve_transmission_eqns(par_table)
  
  DLL_nt<-sum(out_nt_df$Q)
  #DLL_tot_nt = sum(out_nt_df$Q) + sum(N-out_nt_df$Q[out_nt_df$Isymdet/N>close_thres]) # plus days online days*
  n_inf_nt<-out_nt_df$cumI[nrow(out_df)]
  pct_inf_nt<-n_inf_nt/N
  days_of_online_nt<-as.numeric(sum((out_nt_df$Isymdet)/N>close_thres))
  #DLL_averted<-DLL_tot_nt-DLL_tot
  DLL_averted<-DLL_nt- DLL
  inf_averted<-n_inf_nt-n_inf
  pct_inf_averted<- (n_inf_nt-n_inf)/n_inf_nt
  days_of_online_averted<-days_of_online_nt - days_of_online
  pct_detected_nt = out_nt_df$TPs[nrow(out_nt_df)]/out_df$cumI[nrow(out_nt_df)]
  
  cost_nt<-out_nt_df$PCR[nrow(out_nt_df)]*PCR_price
  n_PCR_nt<-out_nt_df$PCR[nrow(out_nt_df)]
  cost_intervention<- (cost_testing - cost_nt) 
  incr_cost_testing<- (cost_testing - cost_nt)
  incr_cost_DLL<- DLL_averted*DLL_price
  incr_cost_inf<- inf_averted*(pct_pos_isofac*isofac_price + contact_tracing_price)
  cost_per_DLL_averted<-cost_intervention/DLL_averted
  cost_per_DO_averted<-cost_intervention/days_of_online_averted
  cost_per_inf_averted<-cost_intervention/inf_averted
  
  
  
  pct_vacc<-par_table$init_vacc
  pct_imm<- 1-(1-init_rec)*(1-init_vacc)
  init_prev<-par_table$init_prev
  
  summary_df<-data.frame(pct_vacc, testing_freq,init_prev, n_inf, n_inf_vax, n_inf_unvax, max_symdet, pct_inf, pct_detected, 
                         n_pos, DLL, days_of_online, cost_of_online, cost_RT, cost_PCR,
                         cost_DLL, cost_testing,cost_isofac, cost_contact_tracing,cost_sequencing,
                         cost_kiosks, n_kiosks, cost_to_UT,
                         n_PCR, n_RT, n_inf_nt, pct_inf_nt, DLL_nt, days_of_online_nt, DLL_averted,
                         inf_averted, pct_inf_averted, days_of_online_averted, cost_nt, n_PCR_nt, pct_detected_nt,
                         cost_intervention, incr_cost_testing, incr_cost_DLL, incr_cost_inf,
                         cost_per_DLL_averted, cost_per_DO_averted, cost_per_inf_averted, n_tests_per_week,
                         days_above_yellow, days_above_orange, days_above_red, days_above_UTJan,
                         cross_yellow, cross_orange, cross_red, cross_UTJan)
  return(summary_df)
}
testing_model_w_uncertainty<-function(par_table, cost_table, risk_tolerance, nsamps, par_bounds, vacc_inf_correlation = 'independent',
                                      test_type = 'antigen'){
  par_names<-colnames(par_table)
  for (i in 1:length(par_names)){
    assign(par_names[i], as.double(par_table[i]))
  }
  
  par_bounds_names<-colnames(par_bounds)
  for (i in 1:length(par_bounds_names)){
    assign(par_bounds_names[i],par_bounds[i])
  }
  

  
  
  
  
  # Make the distributions of parameters to sample from:
  #TRIANGULAR
  R0s<-rtri(n = nsamps, min = min(R0vals), max = max(R0vals), mode = R0) 
  init_recs<-rtri(n= nsamps, min = min(init_rec_vals), max = max(init_rec_vals), mode = init_rec)
  psyms<-rtri(n = nsamps, min = min(psym_vals), max = max(psym_vals), mode = psym) 
  sigma_vs<-rtri(n = nsamps, min = min(sigma_v_vals), max = max(sigma_v_vals), mode = sigma_v)
  sym_reds<-rtri(n = nsamps, min = min(sym_red_vals), max = max(sym_red_vals), mode = sym_red) # #UNIFORM
  iss<-runif(n=nsamps, min = min(is_vals), max = max(is_vals)) 
  is_syms<-runif(n= nsamps, min = min(is_sym_vals), max = max(is_sym_vals)) 
  intros<-runif(n = nsamps, min = min(intro_vals), max = max(intro_vals))
  init_prevs<-rtri(n = nsamps, min = min(init_prev_vals), max = max(init_prev_vals), mode = init_prev)
  
  
  # Run the model for each drawn sample
  par_tablei<-par_table
  for (i in 1:nsamps){
    
    par_tablei$R0<-R0s[i]
    par_tablei$init_rec<-init_recs[i]
    par_tablei$psym<-psyms[i]
    par_tablei$sigma_v<-sigma_vs[i]
    par_tablei$sym_red<-sym_reds[i]
    par_tablei$is<-iss[i]
    par_tablei$is_sym<-is_syms[i]
    par_tablei$intros_per_week<-intros[i]
    par_tablei$init_prev<-init_prevs[i]
    
    if (i==1){
      df_t<-solve_transmission_eqns(par_tablei, vacc_inf_correlation, test_type)
      df<-testing_model(par_tablei, cost_table, risk_tolerance, vacc_inf_correlation, test_type)
      sample<-rep(i, nrow(df_t))
      df_t$sample<-sample
      df$sample<- i
    }
    else{
      df_ti<-solve_transmission_eqns(par_tablei, vacc_inf_correlation, test_type)
      dfi<-testing_model(par_tablei, cost_table, risk_tolerance, vacc_inf_correlation, test_type)
      sample<-rep(i, nrow(df_ti))
      df_ti$sample<-sample
      dfi$sample<- i
      df_t<-rbind(df_t, df_ti)
      df<-rbind(df, dfi)
    }
    
  }
  # df_t is the stacked nsamps*length(t)xncols full simulation data frame
  # df is the stacked nsampsxncols summary otuput data frame
  
  
  
  # MAKE AGGREGATED DATAFRAMES
  # aggr_df_t is now length(t)x 3*length(colnames(df_t) because has header for ub, median, lb
  time<-seq(from = 0, to = duration_school_year, by = 1)
  dates<-seq(from = as.Date("2021-08-20"), by = "day", length.out = duration_school_year+1)
  aggr_df_t<-data.frame(time, dates)
  
  
  col_names<-colnames(df_t)
  
  for (i in 3:(length(col_names)-1)){
    df_long<-df_t%>%dplyr::select(col_names[i], time, sample)
    df_wide<-df_long%>%pivot_wider(names_from = sample, values_from = col_names[i])
    val_mat<-as.matrix(df_wide[,2:ncol(df_wide)])
    med<-apply(val_mat , 1 , quantile , probs = 0.5 , na.rm = TRUE )
    lb<-apply(val_mat , 1 , quantile , probs = 0.025 , na.rm = TRUE )
    ub<-apply(val_mat , 1 , quantile , probs = 0.975 , na.rm = TRUE )
    # put these into data frame 
    med_name<-paste0(col_names[i], '_med')
    lb_name<-paste0(col_names[i], '_lb')
    ub_name<-paste0(col_names[i], '_ub')
    aggr_df_t[[med_name]]<-med
    aggr_df_t[[lb_name]]<-lb
    aggr_df_t[[ub_name]]<-ub
  }
  
  
  
  aggr_df<-data.frame() # all the things in df but with ub, lb, and median + the probability of exceeding each threshold
  run<-1
  aggr_df<-data.frame(run) 
  col_names<-colnames(df)
  for (i in (1:(length(col_names)-5))){
    vec<-df[,i]
    med<-quantile(vec, probs = 0.5, na.rm = TRUE)
    lb<-quantile(vec, probs = 0.05, na.rm = TRUE)
    ub<-quantile(vec, probs = 0.95, na.rm = TRUE)
    # put these into data frame 
    med_name<-paste0(col_names[i], '_med')
    lb_name<-paste0(col_names[i], '_lb')
    ub_name<-paste0(col_names[i], '_ub')
    aggr_df[[med_name]]<-med
    aggr_df[[lb_name]]<-lb
    aggr_df[[ub_name]]<-ub
  }
  
  for (i in ((length(col_names)-4):(length(col_names)-1))){
    vec<-df[,i]
    prob<-sum(vec)/nrow(df)
    # put these into data frame 
    prob_name<-paste0(col_names[i], '_prob')
    aggr_df[[prob_name]]<-prob
  }
  out_list<-list(df, df_t, aggr_df, aggr_df_t)
  return(out_list)
}

get_min_testing_per_vacc<-function(df, threshold_prob, vacc_coverage){

for(j in 1:length(vacc_coverage)){
  df_pv<-df[df$pct_vacc_med == vacc_coverage[j],]
  
  if(j==1){
    
    # find all the testing frequencies that have a probability of exceeding the threshold less than 0.05
    df_prob_orange<-df_pv[df_pv$cross_orange_prob<threshold_prob,]
    df_prob_red<-df_pv[df_pv$cross_red_prob<threshold_prob,]
    df_prob_UTJan<-df_pv[df_pv$cross_UTJan_prob<threshold_prob,]
    
    
    # find the minim
    df_min_orange<-df_prob_orange[which.min(df_prob_orange$testing_freq_med),]
    df_min_red<-df_prob_red[which.min(df_prob_red$testing_freq_med),]
    df_min_UTJan<-df_prob_UTJan[which.min(df_prob_UTJan$testing_freq_med),]
    df_no_test<-df_pv[df_pv$testing_freq_med == 0,]
  }
  else{
    df_prob_orangei<-df_pv[df_pv$cross_orange_prob<threshold_prob,]
    df_prob_redi<-df_pv[df_pv$cross_red_prob<threshold_prob,]
    df_prob_UTJani<-df_pv[df_pv$cross_UTJan_prob<threshold_prob,]
    
    
    df_min_orangei<-df_prob_orangei[which.min(df_prob_orangei$testing_freq_med),]
    df_min_redi<-df_prob_redi[which.min(df_prob_redi$testing_freq_med),]
    df_min_UTJani<-df_prob_UTJani[which.min(df_prob_UTJani$testing_freq_med),]
    df_no_testi<-df_pv[df_pv$testing_freq_med == 0,]
    
    df_min_orange<-rbind(df_min_orange, df_min_orangei)
    df_min_red<-rbind(df_min_red, df_min_redi)
    df_min_UTJan<-rbind(df_min_UTJan, df_min_UTJani)
    df_no_test<-rbind(df_no_test, df_no_testi)
    
    
  }
}
  df_no_test<-df_no_test%>%mutate(thres = "none")
  df_min_orange<-df_min_orange%>%mutate(thres = "orange")
  df_min_red<-df_min_red%>%mutate(thres = "red")
  df_min_UTJan<-df_min_UTJan%>%mutate(thres = "UT Jan")
  df_comb<-rbind( df_min_orange, df_min_red, df_min_UTJan, df_no_test)
  
  return(df_comb)
}

get_pct_infectiousness_removed<-function(V0=3, Vf=6, delta_add=0, t0topeak=2, Vpeak = 9,
                                         tpeaktof=7, LOD=5, LOI=6, dt=0.01, t_test, test_sensitivity){
  
  
  t<-seq(0, 20, dt)
  
  
  # Viral load trajectory
  v<-numeric(length(t))
  growth_rate<-(Vpeak-V0)/t0topeak
  decline_rate<-(Vpeak-Vf)/tpeaktof
  v[t<=t0topeak]<-V0+ growth_rate*(t[t<=t0topeak])
  v[t>t0topeak]<- Vpeak-decline_rate*(t[t>t0topeak]-t0topeak)
  v_i<-v[v>LOI]
  t_rel<-t[v>LOI]
  t_inf<-t_rel- t_rel[1]
  duration_inf<-t_inf[length(t_inf)]
  prop_inf_ind_removed<-(duration_inf-t_test)/duration_inf
  
  v_inf<-v_i-LOI 
  # Find the sum of the AUCs at each stage of infetion (here assuming infectious individuals equally likely to be in any stage)
  AUC_remaining<-rep(0, length(v_inf))
  for (i in 1:length(v_inf)){
    t_stage<-t_inf[i]
    v_rem<-v_inf[t_inf>t_stage]
    AUC_remaining[i]<-sum((v_rem)*dt)
  }
  sum_AUC_all_stages<-sum(AUC_remaining)
  
  # Find the remaining AUC at each stage of infection only including those who would have been caught from a test
  # sum them up, this is the numerator
  v_caught<-v_i[t_inf>t_test] - LOI
  t_caught<-t_inf[t_inf>t_test]
  AUC_remaining_test<-rep(0, length(t_caught))
  for (i in 1:length(v_caught)){
    t_stage<-t_caught[i]
    v_rem_test<-v_caught[t_caught>t_stage]
    AUC_remaining_test[i]<-sum((v_rem_test)*dt)
  }
  sum_AUC_post_test<-sum(AUC_remaining_test)
  pct_removed<-sum_AUC_post_test/sum_AUC_all_stages
  
  eff_percent_removed = test_sensitivity*pct_removed
  df<-data.frame(eff_percent_removed, pct_removed, prop_inf_ind_removed)
  return(df)
  
}
  
