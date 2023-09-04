%% EBM solutions with the diffusivity as a function of the global mean temperature
 
if 0

  gamma = 2:-1:-7; % gamma = d\ln D/ dT_0 in %/K
  
  n = length(gamma);
  T0f = zeros(n,1); % Fit of 0th Legendre polynomial coefficient
  T2f = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient
  T0fw = zeros(n,1); % Fit of 0th Legendre polynomial coefficient for warm climate
  T2fw = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient for warm climate
  h0f = zeros(n,1); % Fit of 0th Legendre polynomial coefficient
  h2f = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient
  h0fw = zeros(n,1); % Fit of 0th Legendre polynomial coefficient for warm climate
  h2fw = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient for warm climate

  %% numerical solutions

  N=180; Nstr = ['N=' num2str(N)];
  %N=60; Nstr = ['N=' num2str(N)];

  dT_ebm_all = zeros(n,N); % stored numerical EBM results

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

  for i=1:1:n % for loop for different gamma values 

	  i

      % control EBM solution
      ebmparams.gamma = 0.0;
      gammastr = ['gamma=' num2str(ebmparams.gamma)];
      Arefstr = ['A=' num2str(ebmparams.A)];
      [t,x,Tebm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
			       'do_D_T0=1',gammastr);
      mse_ebm = calc_mse(Tebm,ebmparams.rh);
 

      % perturbed EBM solution
      ebmparams.gamma = 0.0+0.01*gamma(i);
      gammastr = ['gamma=' num2str(ebmparams.gamma)];
      Awarm = ebmparams.A - ebmparams.forcing; Awarmstr = ['A=' num2str(Awarm)];
      [t,x,Twarmebm,Fwarm,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
				       'alb_P2=1',albostr,albistr,...
				       'do_D_T0=1',gammastr);
      mse_ebm_warm = calc_mse(Twarmebm,ebmparams.rh);

      % compute the numerically solved T0 and T2 using the function compute_T0T2
      [T0f(i),T2f(i)]=compute_T0T2(x, Tebm);
      [T0fw(i),T2fw(i)]=compute_T0T2(x, Twarmebm);
      [h0f(i),h2f(i)]=compute_T0T2(x, mse_ebm);
      [h0fw(i),h2fw(i)]=compute_T0T2(x, mse_ebm_warm);

      % store dT EBM for plotting
      dT_ebm_all(i,:) = Twarmebm-Tebm;
      dh_ebm_all(i,:) = mse_ebm_warm-mse_ebm;

  end
  
  
  %% compute the theoretical estimates

  P2=1/2*(3*x.^2-1);

  % constants
  L = 2.5*10^6;
  H = 0.8;
  cp = 1004.6;
  B = 1.8;
  D0=ebmparams.D0;

  [T0ref,T2ref]=compute_T0T2(x, Tebm);
  dT0=ebmparams.forcing/B;

  [tmp, qp1] = saturation_thermodynamics(T0ref+273.15+0.1,1e5);
  [tmp, q0] = saturation_thermodynamics(T0ref+273.15+0.0,1e5);
  [tmp, qm1] = saturation_thermodynamics(T0ref+273.15-0.1,1e5);
  dqsdTref = (qp1 - qm1)/0.2;
  d2qsdTref2 = (qp1 - 2*q0 + qm1)/0.1^2;
  fref = H*L/cp.*dqsdTref;
  T0warm=T0ref+dT0;
  gamma_c = -L*H/cp*d2qsdTref2/(1+fref);
  [tmp, q0warm] = saturation_thermodynamics(T0warm+273.15+0.0,1e5);
  dh0=dT0+H*L/cp.*(q0warm-q0);

  dT2dT0 = (-6*D0*(L*H/cp*d2qsdTref2+gamma*0.01*(1+fref))./(6*D0*(1+fref)+B))*T2ref;
  dh2dT0 = T2ref*L*H/cp*d2qsdTref2*(6*D0*(1+fref)*gamma*0.01/gamma_c+B)./(6*D0*(1+fref)+B);

  lat = asin(x)*180/pi;

   save('D_as_ftn_of_T0.mat')
else
   load D_as_ftn_of_T0.mat
end

%%
  figure; init_figure;

  fwidth          = 6.5;
  fheight         = fwidth/1.618*.95;

  loff = .07;
  roff = .015;
  boff = .10;
  toff = .02;
  hdist = 0.08;
  vdist = 0.105;
  panel_width = (1 - loff - roff -hdist)/2;
  panel_height = (1-boff-toff-1*vdist)/2;
  popts =[];

  set(gcf, 'Units', 'inches', 'PaperUnits', 'inches')
  psize           = get(gcf, 'PaperSize');
  fpos            = get(gcf, 'Position');
  ppos            = [(psize(1) - fwidth)/2, (psize(2) - fheight)/2,...
				      fwidth, fheight];
  fpos            = [fpos(1) fpos(2) fpos(3) fpos(3)*fheight/fwidth];
  set(gcf, 'PaperPosition', ppos, ...
    'Position', fpos, ...
    'DefaultLineLineWidth', PARS('line_width', popts)*1.5, ...
    'DefaultLineMarkerSize', PARS('marker_size'), ...
    'DefaultAxesLineWidth', PARS('axes_line_width', popts)*1.5, ...
    'DefaultAxesFontSize',  PARS('font_size', popts)*1.3, ...
    'DefaultTextFontSize',  PARS('font_size', popts)*1.3, ...
    'DefaultTextFontName', 'Helvetica', ...
    'DefaultAxesFontName', 'Helvetica', ...
    'DefaultAxesLayer', 'top', ...
    'DefaultAxesTickDir', 'out', ...
    'DefaultAxesTickLength', PARS('tick_length', popts), ...
    'DefaultAxesXMinorTick', 'off', ...
    'DefaultAxesXMinorTick', 'off', ...
    'DefaultAxesYMinorTick', 'off', ...
    'DefaultAxesPosition', [.135 .17 .85 .815] )


  % panel a:                                                          
  subplot('Position',[loff boff+1*(panel_height+vdist) ...
                      panel_width panel_height]);


  gamma_inds = [1:2:9]
  
  for i=1:length(gamma_inds)

    if i==1 plot_col(i,:) = [.15 .15 1]; end
    if i==2 plot_col(i,:) = [0 0 0]; end
    if i>2
      plot_col(i,:) = [1.0 1-0.16*i 1-0.18*i];
    end
    legend_list{i} = [num2str(gamma(gamma_inds(i))) '% K^{-1}']

    plot(lat,dT_ebm_all(gamma_inds(i),:),'color',plot_col(i,:));
    hold on;

  end

  for i=1:length(gamma_inds)
    plot(lat,dT0+dT2dT0(gamma_inds(i))*P2*dT0,'--','color',plot_col(i,:));
  end

  legend(legend_list,'location','N'); legend boxoff;
  set(gca,'box','off');
  xlabel('Lat (deg)','fontweight','bold')
  ylabel('\Delta T (K)','fontweight','bold')
  text(-79,4.8,'a','fontweight','bold')
  axis([-85 85 0 5]);

  % panel b:                                                          
  subplot('Position',[loff boff+0*(panel_height+vdist) ...
                      panel_width panel_height]);


  gamma_inds = [1:2:9]
  
  for i=1:length(gamma_inds)

    if i==1 plot_col(i,:) = [.15 .15 1]; end
    if i==2 plot_col(i,:) = [0 0 0]; end
    if i>2
      plot_col(i,:) = [1.0 1-0.16*i 1-0.18*i];
    end
    legend_list{i} = [num2str(gamma(gamma_inds(i))) '% K^{-1}']

    plot(lat,dh_ebm_all(gamma_inds(i),:),'color',plot_col(i,:));
    hold on;

  end

  for i=1:length(gamma_inds)
    plot(lat,dh0+dh2dT0(gamma_inds(i))*P2*dT0,'--','color',plot_col(i,:));
  end

  legend(legend_list,'location','S'); legend boxoff;
  set(gca,'box','off');
  xlabel('Lat (deg)','fontweight','bold')
  ylabel('\Delta h (K)','fontweight','bold')
  text(-79,11.2,'c','fontweight','bold')
  axis([-85 85 -4 12]);

  % panel c:                                                          
  subplot('Position',[loff+panel_width+hdist boff+1*(panel_height+vdist) ...
			 panel_width panel_height]);

  plot(gamma,T2fw-T2f,'ko-');hold on;                                     
  plot(gamma,dT2dT0.*dT0,'k*:');
  plot(gamma,zeros(size(T2fw)),'k:')

  set(gca,'box','off');
  legend('numerical','est Eq. (12)', ...
	 'location','southeast'); legend boxoff;
  ylabel('\Delta T_2 (K)','fontweight','bold');
  xlabel('\gamma=dln D/dT_0 (% K  )','fontweight','bold')
  text(-1.0,-2.62,'-1','fontweight','bold','fontsize',6)
  axis([-7 2 -2 2.5]);

  text(-6.7,2.3,'b','fontweight','bold')

  % panel d:
  subplot('Position',[loff+panel_width+hdist boff+0*(panel_height+vdist) ...
			 panel_width panel_height]);


  plot(gamma,h2fw-h2f,'ko-');hold on;
  plot(gamma,dh2dT0.*dT0,'k*:');
  plot(gamma,zeros(size(T2fw)),'k:')
  set(gca,'box','off');
  ylabel('\Delta h_2 (K)','fontweight','bold');
  xlabel('\gamma=dln D/dT_0 (% K  )','fontweight','bold')
  text(-1,-10.4,'-1','fontweight','bold','fontsize',6)
  legend('numerical','est Eq. (16)', ...
      'location','southeast'); legend boxoff;
  axis([-7 2 -9 1]);
  text(-6.7,.55,'d','fontweight','bold')

  print -depsc fig2.eps


  return;
