
  N=60; Nstr = ['N=' num2str(N)];


  ebmparams.B = 1.8;
  ebmparams.forcing = 3.6;
  ebmparams.rh = 0.8;
%  ebmparams.albo = 0.32;
%  ebmparams.albi = 0.32;

  % formerly 0.32 
  ebmparams.albo = 0.68;
  ebmparams.albi = -0.2;
  ebmparams.D0 = 0.3;
  ebmparams.A = 210;

  % this is near the critical value using the polar/equator metric
  ebmparams.gamma = -0.027;

  % default case
  ebmparams.gamma = 0.0;

  Dstr = ['D0=' num2str(ebmparams.D0)];
  Bstr = ['B=' num2str(ebmparams.B)];
  rhstr = ['rh=' num2str(ebmparams.rh)];
  albostr = ['albo=' num2str(ebmparams.albo)];
  albistr = ['albi=' num2str(ebmparams.albi)];
  gammastr = ['gamma=' num2str(ebmparams.gamma)];
  
  % old way of calling EBM has same result as before
  if 0
  % control EBM solution
  Arefstr = ['A=' num2str(ebmparams.A)];
  [t,x,Tebm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr);
  mse_ebm = calc_mse(Tebm,ebmparams.rh);

  % perturbed EBM solution
  Awarm = ebmparams.A - ebmparams.forcing; Awarmstr = ['A=' num2str(Awarm)];
  [t,x,Twarmebm,Fwarm,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,'alb_P2=1',albostr,albistr);
  mse_ebm_warm = calc_mse(Twarmebm,ebmparams.rh);
  else
  % new way of calling EBM test with do_D_T0

%  int_D_str = ['do_D_T0=1']
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
  end
  
  figure;  clf;
  plot(x,Twarmebm-Tebm)


