%% EBM solutions with different values of the control diffusivity

if 0
 
  factor = [0:0.25:2]; % factor times the control diffusivity

  n = length(factor);

  T0f = zeros(n,1); % Fit of 0th Legendre polynomial coefficient
  T2f = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient
  T0fw = zeros(n,1); % Fit of 0th Legendre polynomial coefficient for warm climate
  T2fw = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient for warm climate
  h0f = zeros(n,1); % Fit of 0th Legendre polynomial coefficient
  h2f = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient
  h0fw = zeros(n,1); % Fit of 0th Legendre polynomial coefficient for warm climate
  h2fw = zeros(n,1); % Fit of 2nd Legendre polynomial coefficient for warm climate

  T0f_L = T0f; T2f_L = T2f; T0fw_L = T0fw; T2fw_L = T2fw;
  h0f_L = h0f; h2f_L = h2f; h0fw_L = h0fw; h2fw_L = h2fw;

  %% numerical solutions

  N=180; Nstr = ['N=' num2str(N)];
  %N=60; Nstr = ['N=' num2str(N)];

  % default EBM formulation
  ebm_all = zeros(n,N); % stored numerical EBM results
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
 
  for i=1:1:n % for loop for different factor values 

	  i

      % control EBM solution
      ebmparams.D0 = 0.3*factor(i);
      Dstr = ['D0=' num2str(ebmparams.D0)];
      Arefstr = ['A=' num2str(ebmparams.A)];
      [t,x,Tebm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr);
      mse_ebm = calc_mse(Tebm,ebmparams.rh);

      % perturbed EBM solution
      ebmparams.D0 = 0.3*factor(i);
      Dstr = ['D0=' num2str(ebmparams.D0)];
      Awarm = ebmparams.A - ebmparams.forcing; Awarmstr = ['A=' num2str(Awarm)];
      [t,x,Twarmebm,Fwarm,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
				       'alb_P2=1',albostr,albistr);
      mse_ebm_warm = calc_mse(Twarmebm,ebmparams.rh);

      % compute the numerically solved T0 and T2 using the function compute_T0T2
      [T0f(i),T2f(i)]=compute_T0T2(x, Tebm);
      [T0fw(i),T2fw(i)]=compute_T0T2(x, Twarmebm);
      [h0f(i),h2f(i)]=compute_T0T2(x, mse_ebm);
      [h0fw(i),h2fw(i)]=compute_T0T2(x, mse_ebm_warm);

      % store dT EBM for plotting
      ebm_all(i,:) = Tebm;
      dT_ebm_all(i,:) = Twarmebm-Tebm;

  end
  
  % EBM formulation with linearized MSE approximation
  ebm_all_L = zeros(n,N); % stored numerical EBM results
  dT_ebm_all_L = zeros(n,N); % stored numerical EBM results

  ebmparams.B = 1.8;
  ebmparams.forcing = 3.6;
  ebmparams.rh = 0.8; 
  ebmparams.albo = 0.68;
  ebmparams.albi = -0.2;
  ebmparams.A = 210;
  ebmparams.D0 = 0.3;

  % linearizing CC equation is equivalent to rescaling diffusivity
  ebmparams.rh = 0; 
  T0ref=15.4283;
  [tmp, qp1] = saturation_thermodynamics(T0ref+273.15+0.1,1e5);
  [tmp, q0] = saturation_thermodynamics(T0ref+273.15+0.0,1e5);
  [tmp, qm1] = saturation_thermodynamics(T0ref+273.15-0.1,1e5);
  dqsdTref = (qp1 - qm1)/0.2;
  d2qsdTref2 = (qp1 - 2*q0 + qm1)/0.1^2;
  L = 2.5*10^6;H = 0.8;cp = 1004.6;
  fref = H*L/cp.*dqsdTref;
  T0=T0ref+ebmparams.forcing/ebmparams.B;
  [tmp, qp1] = saturation_thermodynamics(T0+273.15+0.1,1e5);
  [tmp, q0] = saturation_thermodynamics(T0+273.15+0.0,1e5);
  [tmp, qm1] = saturation_thermodynamics(T0+273.15-0.1,1e5);
  dqsdT = (qp1 - qm1)/0.2;
  d2qsdT2 = (qp1 - 2*q0 + qm1)/0.1^2;
  f = H*L/cp.*dqsdT;

  Dstr = ['D0=' num2str(ebmparams.D0)];
  Bstr = ['B=' num2str(ebmparams.B)];
  rhstr = ['rh=' num2str(ebmparams.rh)];
  albostr = ['albo=' num2str(ebmparams.albo)];
  albistr = ['albi=' num2str(ebmparams.albi)];
 
  for i=1:1:n % for loop for different factor values 

	  i

      % control EBM solution
      ebmparams.D0 = 0.3*factor(i)*(1+fref); % effective diffusivity
      Dstr = ['D0=' num2str(ebmparams.D0)];
      Arefstr = ['A=' num2str(ebmparams.A)];
      [t,x,Tebm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr);
      mse_ebm = calc_mse(Tebm,ebmparams.rh);
 
      % perturbed EBM solution
      ebmparams.D0 = 0.3*factor(i)*(1+f); % effective diffusivity
      Dstr = ['D0=' num2str(ebmparams.D0)];
      Awarm = ebmparams.A - ebmparams.forcing; Awarmstr = ['A=' num2str(Awarm)];
      [t,x,Twarmebm,Fwarm,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
				       'alb_P2=1',albostr,albistr);
      mse_ebm_warm = calc_mse(Twarmebm,ebmparams.rh);

      % compute the numerically solved T0 and T2 using the function compute_T0T2
      [T0f_L(i),T2f_L(i)]=compute_T0T2(x, Tebm);
      [T0fw_L(i),T2fw_L(i)]=compute_T0T2(x, Twarmebm);
      [h0f_L(i),h2f_L(i)]=compute_T0T2(x, mse_ebm);
      [h0fw_L(i),h2fw_L(i)]=compute_T0T2(x, mse_ebm_warm);

      % store dT EBM for plotting
      ebm_all_L(i,:) = Tebm;
      dT_ebm_all_L(i,:) = Twarmebm-Tebm;

  end
  
  
  %% theoretical estimates

  % constants
  L = 2.5*10^6;
  H = 0.8;
  cp = 1004.6;
  B = 1.8;
  
  % control EBM formulation
  dT2dT0 = zeros(n,1); % stored theoretical estimates
  dh2dT0 = zeros(n,1); % stored theoretical estimates

  for i=1:1:n % for loop for different factor values

  D0=0.3*factor(i);
  Tebm = ebm_all(i,:)';
  [T0ref,T2ref]=compute_T0T2(x, Tebm);
  dT0=ebmparams.forcing/B;

  [tmp, qp1] = saturation_thermodynamics(T0ref+273.15+0.1,1e5);
  [tmp, q0] = saturation_thermodynamics(T0ref+273.15+0.0,1e5);
  [tmp, qm1] = saturation_thermodynamics(T0ref+273.15-0.1,1e5);
  dqsdTref = (qp1 - qm1)/0.2;
  d2qsdTref2 = (qp1 - 2*q0 + qm1)/0.1^2;
  fref = H*L/cp.*dqsdTref;
  T0warm=T0ref+dT0;

  dT2dT0(i) = (-6*D0*(L*H/cp*d2qsdTref2)./(6*D0*(1+fref)+B))*T2ref;
  dh2dT0(i) = T2ref*L*H/cp*d2qsdTref2*B./(6*D0*(1+fref)+B);

  end

  % EBM formulation with linearized MSE approximation
  dT2dT0_L = zeros(n,1); % stored theoretical estimates
  dh2dT0_L = zeros(n,1); % stored theoretical estimates

  for i=1:1:n % for loop for different factor values

  D0=0.3*factor(i);
  Tebm = ebm_all_L(i,:)';
  [T0ref,T2ref]=compute_T0T2(x, Tebm);
  dT0=ebmparams.forcing/B;

  [tmp, qp1] = saturation_thermodynamics(T0ref+273.15+0.1,1e5);
  [tmp, q0] = saturation_thermodynamics(T0ref+273.15+0.0,1e5);
  [tmp, qm1] = saturation_thermodynamics(T0ref+273.15-0.1,1e5);
  dqsdTref = (qp1 - qm1)/0.2;
  d2qsdTref2 = (qp1 - 2*q0 + qm1)/0.1^2;
  fref = H*L/cp.*dqsdTref;
  T0warm=T0ref+dT0;

  dT2dT0_L(i) = (-6*D0*(L*H/cp*d2qsdTref2)./(6*D0*(1+fref)+B))*T2ref;
  dh2dT0_L(i) = T2ref*L*H/cp*d2qsdTref2*B./(6*D0*(1+fref)+B);

  end

  P2=1/2*(3*x.^2-1);
  lat = asin(x)*180/pi;

  save('sensitivity_to_control_D.mat')
 else
   load sensitivity_to_control_D.mat
 end

  %%
  %%
  figure; init_figure;

  fwidth          = 6.5;
  fheight         = fwidth/1.618*.95;

  loff = .07;
  roff = .015;
  boff = .1;
  toff = .02;
  hdist = 0.08;
  vdist = 0.12;
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


  inds = [1:2:9];

  % panel a:                                                          
  subplot('Position',[loff boff+1*(panel_height+vdist) ...
                      panel_width panel_height]);
  
%  colormap(autumn)
%plot_col = summer(5);
plot_col = spring(6);
  for i=1:length(inds)

    if i==1 plot_col(i,:) = [0.5 0.5 0.5]; end
    %if i==2 plot_col(i,:) = [.15 .15 1]; end
    if i==3 plot_col(i,:) = [0 0 0]; end
    if i>3
      %plot_col(i,:) = [1.0 1-0.16*i 1-0.18*i];
    end
    legend_list{i} = [num2str(factor(inds(i))*0.3) ' W m^{-2} K^{-1}']

    plot(lat,dT_ebm_all(inds(i),:),'color',plot_col(i,:));
    %plot(lat,dT_ebm_all(inds(i),:));
    hold on;

  end

  for i=1:length(inds)
    plot(lat,dT0+dT2dT0(inds(i))*P2*dT0,'--','color',plot_col(i,:));
    %plot(lat,dT0+dT2dT0(inds(i))*P2*dT0,'--');
  end

  legend(legend_list,'location','N'); legend boxoff;
  set(gca,'box','off');
  xlabel('Lat (deg)','fontweight','bold')
  ylabel('\Delta T (K)','fontweight','bold')
  text(-79,4.8,'a','fontweight','bold')
  %text(0,2.3,'$h$: Eq. (1)','fontweight','bold','HorizontalAlignment', 'center','interpreter','latex')
%  text(-17,2.3,'$h$','fontweight','bold','interpreter','latex')
%  text(-10,2.3,': Eq. (1)')
%  text(-16,5.,'$h$','fontweight','bold','interpreter','latex')
%  text(-10,5.,': Eq. (1)')
%  text(-36,5.,'EBM with ','fontweight','bold')
  text(-20.5,5.,'EBM with ','fontweight','bold')
  text(+16.7,4.99,'$h$','fontweight','bold','interpreter','latex')

%  text(0,2.3,'$h$','fontweight','bold','HorizontalAlignment', 'center','interpreter','latex')
%  text(0,2.6,'h: Eq. (1)','fontweight','bold','HorizontalAlignment', 'center')
  axis([-85 85 0 5]);

  % panel c:                                                          
  subplot('Position',[loff boff+0*(panel_height+vdist) ...
                      panel_width panel_height]);
  
  for i=1:length(inds)

    plot(lat,dT_ebm_all_L(inds(i),:),'color',plot_col(i,:));

    hold on;

  end

  for i=1:length(inds)
    plot(lat,dT0+dT2dT0_L(inds(i))*P2*dT0,'--','color',plot_col(i,:));
  end

  legend(legend_list,'location','N'); legend boxoff;
  set(gca,'box','off');
  xlabel('Lat (deg)','fontweight','bold')
  ylabel('\Delta T (K)','fontweight','bold')
  text(-79,4.8,'c','fontweight','bold')
%  text(0,2.3,'$\tilde{h}$: Eq. (5)','fontweight','bold','HorizontalAlignment', 'center','interpreter','latex')
%  text(-16,5,'$\tilde{h}$','fontweight','bold','interpreter','latex')
%  text(-10,5,': Eq. (4)')
  text(-20.5,5.,'EBM with ','fontweight','bold')
  text(+16.7,5.02,'$\tilde{h}$','fontweight','bold','interpreter','latex')
  axis([-85 85 0 5]);

  % panel b:
  subplot('Position',[loff+panel_width+hdist boff+1*(panel_height+vdist) ...
			 panel_width panel_height]);
  
  plot(factor*0.3,T2fw-T2f,'ko-');hold on;                                     
  plot(factor*0.3,dT2dT0.*dT0,'k*:');

  set(gca,'box','off');
  legend('numerical','est Eq. (7)', ...
      location='southeast'); legend boxoff;
  ylabel('\Delta T_2 (K)','fontweight','bold');
%xlabel('$\overline{D}$ (Wm^{-2}K^{-1})','fontweight','bold','interpreter','latex')
%text(.25,-0.4,'$\overline{D}$ (Wm^{-2}K^{-1})','interpreter','latex')
%xlabel(.25,-0.4,'\overline{D} (W m^{-2} K^{-1})','fontweight','bold')
xlabel('$\overline{D} (\mathrm{W m^{-2} K^{-1}})$','fontweight','bold','interpreter','latex')
  text(0.02,1.95,'b','fontweight','bold')
  axis([0 0.6 0 2]);


  % panel d:
  subplot('Position',[loff+panel_width+hdist boff+0*(panel_height+vdist) ...
			 panel_width panel_height]);

  plot(factor*0.3,h2fw-h2f,'ko-');hold on;                                     
  plot(factor*0.3,dh2dT0.*dT0,'k*:');

  set(gca,'box','off');
  legend('numerical','est Eq. (10)', ...
      location='southeast'); legend boxoff;
  ylabel('\Delta h_2 (K)','fontweight','bold');
%  text(.25,-0.4,'$\overline{D}$ (Wm^{-2}K^{-1})','interpreter','latex')
%  xlabel('$\overline{D}$ (Wm^{-2}K^{-1})','fontweight','bold')
xlabel('$\overline{D} (\mathrm{W m^{-2} K^{-1}})$','fontweight','bold','interpreter','latex')
  text(0.02,-1,'d','fontweight','bold')
  axis([0 0.6 -30 0]);

  print -depsc fig1.eps

  return;

