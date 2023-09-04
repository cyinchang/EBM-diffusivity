%% EBM solutions with the diffusivity as a function of T2 and h2

if 0

  k=[0,1,1.5,2,3]; % exponent of T_2 and h_2, i.e., (T_2/T_2ref)^k*(h_2/h_2ref)^k

  n = length(k);
  T0mf_B = zeros(n,2); % Fit of 0th Legendre polynomial coefficient 
  T2mf_B = zeros(n,2); % Fit of 2nd Legendre polynomial coefficient 
  T0mf_T_B = zeros(n,2); % Fit of 0th Legendre polynomial coefficient 
  T2mf_T_B = zeros(n,2); % Fit of 2nd Legendre polynomial coefficient 
  T0mf_h_B = zeros(n,2); % Fit of 0th Legendre polynomial coefficient 
  T2mf_h_B = zeros(n,2); % Fit of 2nd Legendre polynomial coefficient 
  h0mf_B = zeros(n,2); % Fit of 0th Legendre polynomial coefficient 
  h2mf_B = zeros(n,2); % Fit of 2nd Legendre polynomial coefficient 
  h0mf_T_B = zeros(n,2); % Fit of 0th Legendre polynomial coefficient 
  h2mf_T_B = zeros(n,2); % Fit of 2nd Legendre polynomial coefficient 
  h0mf_h_B = zeros(n,2); % Fit of 0th Legendre polynomial coefficient 
  h2mf_h_B = zeros(n,2); % Fit of 2nd Legendre polynomial coefficient 

  %% numerical solutions

  %N=60; Nstr = ['N=' num2str(N)];
  N=180; Nstr = ['N=' num2str(N)];

  ebmparams.gamma = 0.0;
  ebmparams.forcing = 3.6;
  ebmparams.albo = 0.68;
  ebmparams.albi = -0.2;
  ebmparams.A = 210;
  ebmparams.D0 = 0.3;
  
  albostr = ['albo=' num2str(ebmparams.albo)];
  albistr = ['albi=' num2str(ebmparams.albi)];
  Dstr = ['D0=' num2str(ebmparams.D0)];

  % reference moist climate
  ebmparams.rh = 0.8; % moist EBM 
  rhstr = ['rh=' num2str(ebmparams.rh)];
  ebmparams.B = 1.8;
  Bstr = ['B=' num2str(ebmparams.B)];
  Arefstr = ['A=' num2str(ebmparams.A)];
  [t,x,Trefm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr);

  for i=1:1:n % for loop for different k values 

      % moist climates
      ebmparams.rh = 0.8; % moist EBM
      rhstr = ['rh=' num2str(ebmparams.rh)];
      Arefstr = ['A=' num2str(ebmparams.A)];
      Awarm = ebmparams.A - ebmparams.forcing; Awarmstr = ['A=' num2str(Awarm)];

      % D as a function of T2 and h2
      tic
      [t,x,Tebm,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
                   ['do_D_T2=',num2str(k(i))],['do_D_h2=',num2str(k(i))],...
                   'D_Tref_input=1', Trefm);
      toc    
      mse = calc_mse(Tebm,ebmparams.rh);
      tic
      [t,x,Tebmwarm,F,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
                   ['do_D_T2=',num2str(k(i))],['do_D_h2=',num2str(k(i))],...
                   'D_Tref_input=1', Trefm);
      toc    
      msewarm = calc_mse(Tebmwarm,ebmparams.rh);

      % D as a function of T2
      tic
      [t,x,Tebm_T,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
                   ['do_D_T2=',num2str(k(i))], 'D_Tref_input=1', Trefm);
      toc     
      mse_T = calc_mse(Tebm_T,ebmparams.rh);
      tic
      [t,x,Tebmwarm_T,F,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
                   ['do_D_T2=',num2str(k(i))], 'D_Tref_input=1', Trefm);
      toc     
      msewarm_T = calc_mse(Tebmwarm_T,ebmparams.rh);

      % D as a function of h2
      tic
      [t,x,Tebm_h,F,SW,alb]=moistEBM(Arefstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
			       ['do_D_h2=',num2str(k(i))], 'D_Tref_input=1', Trefm);
      toc                    
      mse_h = calc_mse(Tebm_h,ebmparams.rh);
      tic
      [t,x,Tebmwarm_h,F,SW,alb]=moistEBM(Awarmstr, Dstr,Bstr,rhstr,Nstr,...
			       'alb_P2=1',albostr,albistr,...
			       ['do_D_h2=',num2str(k(i))], 'D_Tref_input=1', Trefm);
      toc                    
      msewarm_h = calc_mse(Tebmwarm_h,ebmparams.rh);

      % compute the numerically solved T0 and T2 using the function compute_T0T2
      [T0mf_B(i,1),T2mf_B(i,1)]=compute_T0T2(x, Tebm);
      [T0mf_T_B(i,1),T2mf_T_B(i,1)]=compute_T0T2(x, Tebm_T);
      [T0mf_h_B(i,1),T2mf_h_B(i,1)]=compute_T0T2(x, Tebm_h);
      [h0mf_B(i,1),h2mf_B(i,1)]=compute_T0T2(x, mse);
      [h0mf_T_B(i,1),h2mf_T_B(i,1)]=compute_T0T2(x, mse_T);
      [h0mf_h_B(i,1),h2mf_h_B(i,1)]=compute_T0T2(x, mse_h); 
      [T0mf_B(i,2),T2mf_B(i,2)]=compute_T0T2(x, Tebmwarm);
      [T0mf_T_B(i,2),T2mf_T_B(i,2)]=compute_T0T2(x, Tebmwarm_T);
      [T0mf_h_B(i,2),T2mf_h_B(i,2)]=compute_T0T2(x, Tebmwarm_h);
      [h0mf_B(i,2),h2mf_B(i,2)]=compute_T0T2(x, msewarm);
      [h0mf_T_B(i,2),h2mf_T_B(i,2)]=compute_T0T2(x, msewarm_T);
      [h0mf_h_B(i,2),h2mf_h_B(i,2)]=compute_T0T2(x, msewarm_h); 

      dTebm_B(i,:) = Tebmwarm-Tebm;
      dTebm_T_B(i,:) = Tebmwarm_T-Tebm_T;
      dTebm_h_B(i,:) = Tebmwarm_h-Tebm_h;
  end

%% compute the estimated solutions

% constants
L = 2.5*10^6;
H = 0.8;
cp = 1004.6;
B = 1.8;

Q=1360;
Sa2=ebmparams.albi-0.482*ebmparams.albo;
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

n_h=k;n_T=k;
dT2dT0_B = (-6*(D0*L*H/cp*d2qsdTref2*(n_h+1))./(6*D0*(1+fref)*(n_h+n_T+1)+B))*T2ref;
dh2dT0_B = (1+fref)*dT2dT0_B+L*H/cp*d2qsdTref2*T2ref;
dlnDdT0_B = -gamma_c./(6*D0*(1+fref)*(n_h+n_T+1)+B).*(n_h*B-n_T*6*D0*(1+fref));

n_h=k*0;n_T=k;
dT2dT0_T_B = (-6*(D0*L*H/cp*d2qsdTref2*(n_h+1))./(6*D0*(1+fref)*(n_h+n_T+1)+B))*T2ref;
dh2dT0_T_B = (1+fref)*dT2dT0_T_B+L*H/cp*d2qsdTref2*T2ref;
dlnDdT0_T_B = -gamma_c./(6*D0*(1+fref)*(n_h+n_T+1)+B).*(n_h*B-n_T*6*D0*(1+fref));

n_h=k;n_T=k*0;
dT2dT0_h_B = (-6*(D0*L*H/cp*d2qsdTref2*(n_h+1))./(6*D0*(1+fref)*(n_h+n_T+1)+B))*T2ref;
dh2dT0_h_B = (1+fref)*dT2dT0_h_B+L*H/cp*d2qsdTref2*T2ref;
dlnDdT0_h_B = -gamma_c./(6*D0*(1+fref)*(n_h+n_T+1)+B).*(n_h*B-n_T*6*D0*(1+fref));

Dm_B=(T2mf_B(:,2)./T2mf_B(:,1).*h2mf_B(:,2)./h2mf_B(:,1))'.^k;
Dm_T_B=(T2mf_T_B(:,2)./T2mf_T_B(:,1))'.^k;
Dm_h_B=(h2mf_h_B(:,2)./h2mf_h_B(:,1))'.^k;     

lat = asin(x)*180/pi;
P2=1/2*(3*x.^2-1);

 save('D_as_ftn_of_T2h2.mat')
else
  load D_as_ftn_of_T2h2.mat
end

% in text number n=3/2 T2 change vs. n=0
% calculation below incomplete
diff(T2mf_T_B,1,2)


figure; init_figure;

fwidth          = 6.5;
fheight         = fwidth/2/1.618*3;

loff = .065;
roff = .015;
boff = .075;
toff = .03;
hdist = 0.07;
vdist = 0.09;
panel_width = (1 - loff - roff -hdist)/2;
panel_height = (1-boff-toff-2*vdist)/3;
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

  brown = [0.5 0.22 0.08];
  purple = [0.5 0.2 1.0];
  orange = [1.0 0.5 0.0];


% panel a: dT vs lat for T2 dependent diffusivity
subplot('Position',[loff boff+2*(panel_height+vdist) ...
			 panel_width panel_height]);

% fill in dT vs lat here
k_inds = [1,3,5];
plot_sty = {'-','-','--'};

% want D fn to be all of same color and linestyle to vary with exponent
  plot(lat,dTebm_B(k_inds(1),:),plot_sty{1},'color','k');
  hold on;
for i=2:length(k_inds)
  plot(lat,dTebm_T_B(k_inds(i),:),plot_sty{i},'color',orange);
end

for i=2:length(k_inds) 
  plot(lat,dT0+dT2dT0_T_B(k_inds(i))*P2*dT0,plot_sty{i},'color',[0.6 0.6 0.6]);
end

%legend(legend_list,'location','N'); legend boxoff;
legend({'k=0','k=3/2','k=3','est Eq. (19)'},'location','SE'); legend boxoff;
text(36, 4.0, 'D(T_2)','fontweight','bold','fontsize',10.5,'color',orange)
set(gca,'box','off');
xlabel('Lat (deg)','fontweight','bold')
ylabel('\Delta T (K)','fontweight','bold')

text(3,3.7,'a','fontweight','bold')
%axis([-85 85 0 5]);
axis([0 85 0 4]);


% panel c: dT vs lat for h2 dependent diffusivity
subplot('Position',[loff boff+1*(panel_height+vdist) ...
			 panel_width panel_height]);

% fill in dT vs lat here
k_inds = [1,3,5];
plot_sty = {'-','-','--'};

% want D fn to be all of same color and linestyle to vary with exponent
  plot(lat,dTebm_B(k_inds(1),:),plot_sty{1},'color','k');
  hold on;
for i=2:length(k_inds)
  plot(lat,dTebm_h_B(k_inds(i),:),plot_sty{i},'color',purple);
end

for i=2:length(k_inds)
  plot(lat,dT0+dT2dT0_h_B(k_inds(i))*P2*dT0,plot_sty{i},'color',[.6 .6 .6]);
end

%title('D(h_2)','fontweight','bold')
text(36, 4.0, 'D(h_2)','fontweight','bold','fontsize',10.5,'color',purple)
legend({'k=0','k=3/2','k=3','est Eq. (19)'},'location','SE'); legend boxoff;
set(gca,'box','off');
xlabel('Lat (deg)','fontweight','bold')
ylabel('\Delta T (K)','fontweight','bold')

text(3,3.7,'c','fontweight','bold')
%axis([-85 85 0 5]);
axis([0 85 0 4]);


% panel e: dT vs lat for T2,h2 dependent diffusivity
subplot('Position',[loff boff+0*(panel_height+vdist) ...
			 panel_width panel_height]);

% fill in dT vs lat here
k_inds = [1,3,5];
plot_sty = {'-','-','--'};


% want D fn to be all of same color and linestyle to vary with exponent
  plot(lat,dTebm_B(k_inds(1),:),plot_sty{1},'color','k');
  hold on;
for i=2:length(k_inds)
  plot(lat,dTebm_B(k_inds(i),:),plot_sty{i},'color',brown);
end

for i=2:length(k_inds)
  plot(lat,2+dT2dT0_B(k_inds(i))*P2*2,plot_sty{i},'color',[0.6 0.6 0.6]);
end

text(33, 4.0, 'D(T_2,h_2)','fontweight','bold','fontsize',10.5,'color',brown)
legend({'k=0','k=3/2','k=3','est Eq. (19)'},'location','SE'); legend boxoff;
%legend(legend_list,'location','N'); legend boxoff;
set(gca,'box','off');
xlabel('Lat (deg)','fontweight','bold')
ylabel('\Delta T (K)','fontweight','bold')

text(3,3.7,'e','fontweight','bold')
axis([0 85 0 4]);
%axis([-85 85 0 5]);


% panel b:
subplot('Position',[loff+panel_width+hdist boff+2*(panel_height+vdist) ...
			 panel_width panel_height]);

plot(k,diff(T2mf_T_B,1,2),'o-','color',orange);hold on;                  
plot(k,dT2dT0_T_B*dT0,'*:','color',orange);hold on;                      
plot(k,diff(T2mf_B,1,2),'o-','color',brown);hold on;                   
plot(k,dT2dT0_B*dT0,'*:','color',brown);hold on;                         
plot(k,diff(T2mf_h_B,1,2),'o-','color',purple);hold on;                 
plot(k,dT2dT0_h_B*dT0,'*:','color',purple);hold on;                        

ylim([0,2]);                                 

% percent errors
(diff(T2mf_h_B,1,2)-dT2dT0_h_B'*dT0)./diff(T2mf_h_B,1,2)
(diff(T2mf_T_B,1,2)-dT2dT0_T_B'*dT0)./diff(T2mf_T_B,1,2)
(diff(T2mf_B,1,2)-dT2dT0_B'*dT0)./diff(T2mf_B,1,2)

%legend('numerical D(T_2^kh_2^k)','numerical D(T_2^k)','numerical D(h_2^k)', ...
%       'estimated from #18',...
%       location='southwest', NumColumns=2);                                   
%legend boxoff;                                                   
 set(gca,'box','off');                                                      
  ylabel('\Delta T_2 (K)','fontweight','bold');                                    
  xlabel('k','fontweight','bold')                        
%  axis([-7 2 -2 2.5]); 
  text(0.1,1.85,'b','fontweight','bold')                                     
  legend('numerical','est Eq. (19)','location','SW'); legend boxoff;

% panel d:
subplot('Position',[loff+panel_width+hdist boff+1*(panel_height+vdist) ...
			 panel_width panel_height]);

plot(k,diff(h2mf_h_B,1,2),'o-','color',purple);hold on;                      
plot(k,dh2dT0_h_B*dT0,'*:','color',purple);hold on;            
plot(k,diff(h2mf_B,1,2),'o-','color',brown);hold on;           
plot(k,diff(h2mf_T_B,1,2),'o-','color',orange);hold on;              

plot(k,dh2dT0_B*dT0,'*:','color',brown);hold on;         
plot(k,dh2dT0_T_B*dT0,'*:','color',orange);hold on;                     

% percent errors
(diff(h2mf_h_B,1,2)-dh2dT0_h_B'*dT0)./diff(h2mf_h_B,1,2)
(diff(h2mf_T_B,1,2)-dh2dT0_T_B'*dT0)./diff(h2mf_T_B,1,2)
(diff(h2mf_B,1,2)-dh2dT0_B'*dT0)./diff(h2mf_B,1,2)

ylim([-5,0]);                                 
%legend('numerical D(T_2^kh_2^k)','numerical D(T_2^k)','numerical D(h_2^k)', ...
%       'estimated from #18',...
%       location='southwest', NumColumns=2);                                   
%legend boxoff;                                                   
 set(gca,'box','off');                                                      
  ylabel('\Delta h_2 (K)','fontweight','bold');                                    
  xlabel('k','fontweight','bold')                        
%  axis([-7 2 -2 2.5]); 
  text(.1,-.375,'d','fontweight','bold')                                     
  legend('numerical','est Eq. (20)','location','SW'); legend boxoff;

% panel f: diffusivity
subplot('Position',[loff+panel_width+hdist boff+0*(panel_height+vdist) ...
			 panel_width panel_height]);

  plot(k,(Dm_B-1)*100/2,'o-','color',brown);hold on;                           
  plot(k,dlnDdT0_B*dT0*100/2,'*:','color',brown);hold on;                     
  plot(k,(Dm_T_B-1)*100/2,'o-','color',orange);hold on;                   
  plot(k,(Dm_h_B-1)*100/2,'o-','color',purple);hold on;                

  plot(k,dlnDdT0_T_B*dT0*100/2,'*:','color',orange);hold on;               
  plot(k,dlnDdT0_h_B*dT0*100/2,'*:','color',purple);hold on;                
  plot(k,zeros(size(k)),'k:');
  %ylabel('\Delta D/D (% K^{-1})','fontweight','bold');
  ylabel('\Delta    /    (% K^{-1})','fontweight','bold');
%  ylabel('\Delta   /   (% K  )','fontweight','bold');
%  htt=text(-.38,1.06,'-1','fontweight','bold','fontsize',PARS('font_size', popts));
%  set(htt,'Rotation',90);
  ht=text(-.465,-.45,'\_','fontweight','bold');
  set(ht,'Rotation',90);
  ht1=text(-.317,-1.2,'$\mathcal{D}$','fontweight','bold','interpreter','latex');
  set(ht1,'Rotation',90);
  ht2=text(-.317,-.54,'$\mathcal{D}$','fontweight','bold','interpreter','latex');
  set(ht2,'Rotation',90);
  ylim([-3,3]);         
%  legend('numerical D(T_2^kh_2^k)','numerical D(T_2^k)','numerical D(h_2^k)',...                                                                          

%      'estimated from #30 m=n=k','estimated from #30 n=k','estimated from #30 m=k', ...                                                                   
%         location='northwest', NumColumns=2);                             
                                                                          
           
%ylim([-5,0]);                                 
%legend('numerical D(T_2^kh_2^k)','numerical D(T_2^k)','numerical D(h_2^k)', ...
%       'estimated from #18',...
%       location='southwest', NumColumns=2);                                   
%legend boxoff;                                                   
 set(gca,'box','off');                                                      
%  ylabel('D','fontweight','bold');                                    
  xlabel('k','fontweight','bold')                        
  text(.1,5.1,'f','fontweight','bold')                                     
  legend('numerical','est Eq. (21)','location','SW'); legend boxoff;
   
print -depsc fig3.eps                                                
return


