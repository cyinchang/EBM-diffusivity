%% Truncated EBM solutions based on an alternative formulation of the linearized MSE 
 
  % compare the numerically solved versus truncated solutions 
  % numerical solution: as in Figs. 1-2 in Merlis and Henry (2018)
  % truncated solution: derived in may28_analysis 

  %% copied and modified from the codes in driver_script.m to solve EBM

  N=60; Nstr = ['N=' num2str(N)];

  ebmparams.B = 1.8;
  ebmparams.forcing = 3.6;
  ebmparams.rh = 0.8; 
  ebmparams.albo = 0.68;
  ebmparams.albi = -0.2;
  ebmparams.A = 210;
  ebmparams.D0 = 0.3;
  
  Dstr = ['D0=' num2str(ebmparams.D0)];
  Bstr = ['B=' num2str(ebmparams.B)];
  rhstr = ['rh=' num2str(ebmparams.rh)];
  albostr = ['albo=' num2str(ebmparams.albo)];
  albistr = ['albi=' num2str(ebmparams.albi)];

  ebmparams.gamma = 0.0;
  gammastr = ['gamma=' num2str(ebmparams.gamma)];

  % control EBM solution
  Arefstr = ['A=' num2str(ebmparams.A)];
  [t,x,Tebm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
			       'do_D_T0=1',gammastr);
  mse_ebm = calc_mse(Tebm,ebmparams.rh);
 

  % perturbed EBM solution
  Awarm = ebmparams.A - ebmparams.forcing; Awarmstr = ['A=' num2str(Awarm)];
  [t,x,Twarmebm,Fwarm,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
				    'alb_P2=1',albostr,albistr,...
				    'do_D_T0=1',gammastr);
  mse_ebm_warm = calc_mse(Twarmebm,ebmparams.rh);

  % compute the numerically solved T0 and T2 using the function compute_T0T2
  [T0,T2]=compute_T0T2(x, Tebm);
  [T0warm,T2warm]=compute_T0T2(x, Twarmebm);

  [h0,h2]=compute_T0T2(x, mse_ebm);
  [h0warm,h2warm]=compute_T0T2(x, mse_ebm_warm);


  %% compute the estimated T(x) and h(x) from the truncated EBM solution 

  Q=1360;
  Sa2=ebmparams.albi-0.482*ebmparams.albo;

  [tmp, q2] = saturation_thermodynamics(T0+273.15+0.1,1e5);
  [tmp, q1] = saturation_thermodynamics(T0+273.15-0.1,1e5);
  dqsdTref = (q2 - q1)/0.2;
  fref = H*L/cp.*dqsdTref;
  eta_est_T0 = (1+fref);
  T2_est=0.25*Q*Sa2./(6*ebmparams.D0*eta_est_T0+ebmparams.B);

  [tmp, q2] = saturation_thermodynamics(T0warm+273.15+0.1,1e5);
  [tmp, q1] = saturation_thermodynamics(T0warm+273.15-0.1,1e5);
  dqsdTwarm = (q2 - q1)/0.2;
  fwarm = H*L/cp.*dqsdTwarm;
  eta_warm_est_T0 = (1+fwarm);
  T2_warm_est=0.25*Q*Sa2/(6*ebmparams.D0*eta_warm_est_T0+ebmparams.B);

  T0_est=T0;
  T0_warm_est=T0_est+ebmparams.forcing/ebmparams.B;

  T_est=T0_est+P2*T2_est;
  Twarm_est=T0_warm_est+P2*T2_warm_est;

  mse_est= calc_mse(T_est,ebmparams.rh);
  mse_warm_est= calc_mse(Twarm_est,ebmparams.rh);

  [h0_est,h2_est]=compute_T0T2(x, mse_est);
  [h0_warm_est,h2_warm_est]=compute_T0T2(x, mse_warm_est);

  %% plot the results

  lat=asin(x)*180/pi;
  P2=(3*x.*x-1)/2;

  figure(1);
  subplot(2,1,1);
  plot(lat,Twarmebm-Tebm,'k-');hold on;
  %plot(lat,(T0warm-T0)+P2*(T2warm-T2),'k--');
  plot(lat,Twarm_est-T_est,'r-');
  %legend('numerical','numerical: P2 comp.','truncated')
  legend('numerical','truncated')
  xlim([-90,90]);ylim([0,5]);
  xlabel('Latitude');ylabel('\Delta T (K)');
  text(-80,4.6,['numerical T_2 = ',num2str(T2,'%3.1f'),' -> ',num2str(T2warm,'%3.1f')])
  text(-80,4.2,['truncated T_2 =',num2str(T2_est,'%3.1f'),' -> ',num2str(T2_warm_est,'%3.1f')]);

  subplot(2,1,2);
  plot(lat,mse_ebm_warm-mse_ebm,'k-');hold on;
  plot(lat,mse_warm_est-mse_est,'r');
  legend('numerical','truncated with')
  xlim([-90,90]);ylim([0,9]);
  xlabel('Latitude');ylabel('\Delta h (K)'); 
  text(-80,8,['numerical h_2 = ',num2str(h2,'%3.1f'),' -> ',num2str(h2warm,'%3.1f')])
  text(-80,7.5,['truncated h_2 =',num2str(h2_est,'%3.1f'),' -> ',num2str(h2_warm_est,'%3.1f')]);
  text(-80,7,['truncated \tilde{\tilde{h}}_2 =',num2str(T2_est*eta_est_T0,'%3.1f'),' -> ',num2str(T2_warm_est*eta_warm_est_T0,'%3.1f')]);
