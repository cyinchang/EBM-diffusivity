function [t_o,x_o,T_o,dfn_o,SW_o,alb_o]=moistEBM(varargin)
%   [t,x,T,dfn,SW,alb_o]=albfeedback_EBM(arguments)
% Sellers EBM with surface albedo feedback
% Ex: N=15; [t,x,T]=sellers_EBM(['N=' num2str(N)],'dur=30');
%     subplot(1,2,1), plot(t,mean(T,2)), subplot(1,2,2), plot(x,T(end,:))
% based on code from Till Wagner and Ian Eisenman

global A B SW D0 rh c N x albo albi dfn SW_input dfn_input alb_P2 albref Tref Tref_input alb do_D_T0 do_D_T2 do_D_h2 gamma D_Tref_input D_Tref

% === Parameters ===
% Measuring time in years for d/dt, using W/m^2 for fluxes and then
% multiplying d/dt by number of seconds in year.
N=120; % number of spatial grid points
dur=30; % duration of integration

%c=6.3; % heat capacity of 50m ocean mixed layer (W/m^2*yr)
c=6.342;
%D=0.3; % diffusivity (W/m2/C)
D0=0.15; % diffusivity (W/m2/C)
%B=1.75;
B=1.75;
A=203.3;
rh = 0.8;

% albedos
albo = 0.32;
%albi = 0.62;
albi = 0.32;

% sensitivity of diffusivity to global mean T
gamma = 0.0;

% flags for input arguments
SW_input = 0;
dfn_input = 0;
T0_input = 0;
Tref_input = 0;
alb_P2 = 0;
albref = 0;
do_D_T0 = 0;
do_D_T2 = 0;
do_D_h2 = 0;
D_Tref_input = 0;


T0=15; % initial value of T (degC) (same for all points)
Tref=0; % value of T for linear CC

RelTol=1e-8; % for ODE solver
AbsTol=1e-10; % for ODE solver

RelTol=1e-7; % for ODE solver
AbsTol=1e-9; % for ODE solver

%RelTol=1e-4; % for ODE solver
%AbsTol=1e-4; % for ODE solver

%RelTol=1e-2; % for ODE solver
%AbsTol=1e-2; % for ODE solver

% === EDIT PARAMETERS BELOW HERE ===
% parameter value changes as input to function (batch mode)
% assumed form for input arguments is flag variable
% followed immediately by that variable's input array: 'SW_input=1',SWref
if nargin>0 
 for j=1:length(varargin)
  if isstr(varargin{j})
   eval([varargin{j} ';'])
   if findstr(varargin{j},'SW_input=1')
     SW = varargin{j+1};
   end
   if findstr(varargin{j},'dfn_input=1')
     dfn = varargin{j+1};
   end
   if findstr(varargin{j},'alb_P2=1')
     albref = varargin{j+1};
   end
   if findstr(varargin{j},'T0_input=1')
     T0 = varargin{j+1};
   end
   if findstr(varargin{j},'Tref_input=1')
     Tref = varargin{j+1};
   end
   if findstr(varargin{j},'D_Tref_input=1')
     D_Tref = varargin{j+1};
   end
  end
 end
end

% Tref used for linear albedo
%Tref = T0;

%x = [0:1/N:1.0-1/N/2]';
%x = [-1.0+1/N/2:2/N:1.0-1/N/2]';
x = [-1.0+1/N:2/N:1.0-1/N]';

% === integration ===
if length(T0)==1, T0=repmat(T0,N,1); end % initial value of T
%if length(D)==1, D=repmat(D,N+1,1); end % diffusivity array
% run integration one yr at time until max duration
[t,T]=ode45(@dTdt,[0 dur],T0,odeset('RelTol',RelTol,'AbsTol',AbsTol));

% === plotting or output ===
if nargout>0 % output
    t_o=t; T_o=T(end,:)'; x_o=x; dfn_o=dfn; SW_o=SW; alb_o=alb;
else % plotting
    save tmp.mat t x T N
    % plot evolution of model state: entire time series
    figure(2), clf
    subplot(1,2,1)
    plot(t,T), title('T')
    axis tight, grid on
    ylabel('T (degC)')
    subplot(1,2,2)
    hold on
    plot(x,T(end,:))
    %disp(mean(T(end,:)))
end


% ===================


% === model equations ===
function F=dTdt(~,T)
global A B SW D0 rh c N x albo albi dfn SW_input dfn_input alb alb_P2 albref Tref Tref_input do_D_T0 do_D_T2 do_D_h2 gamma D_Tref_input D_Tref

if ~SW_input
 if alb_P2
  [SW,alb]=SW_P2(x,albo,albi);
 else
  [SW,alb]=SW_T(x,T,albo,albi);
%  [SW,alb]=SW_py(x,T,albo,albi);
 end
end

%Tk = T+273.15;
%[tmp,qs]=saturation_thermodynamics(Tk,1e5);
%mse = PARS('cp')*Tk + PARS('l_cond')*qs*.8;

% MSE in units of T
%mse = Tk + PARS('l_cond')*qs*.8/PARS('cp');

if Tref_input
  mse = calc_mse(T,rh,Tref_input,Tref); 
else
  mse = calc_mse(T,rh); 
end

% interactive diffusivity
if D_Tref_input
  D = calc_D(D0,T,D_Tref,x,do_D_T0,do_D_T2,do_D_h2,gamma,rh);
  D = repmat(D,N+1,1);
else
  % detailed numbers taken from control D0=0.3 mebm
  %Tref = 15.4253 - 29.32*(1/2*(3*x.^2-1));
  D_Tref = 15.4283 - 29.3316*(1/2*(3*x.^2-1));
  
  D = calc_D(D0,T,D_Tref,x,do_D_T0,do_D_T2,do_D_h2,gamma,rh);
  D = repmat(D,N+1,1);
end

% dry test
% this should be equivalent to input rh=0
%mse = Tk + 0.*PARS('l_cond')*qs*.8/PARS('cp');

% linearization as in Rose 14
%[tmp,qs2]=saturation_thermodynamics(Tk+0.1,1e5);
%dqsdT = (qs2-qs)/.1;
%f = PARS('l_cond')*dqsdT*.8/PARS('cp');
%mse2 = Tk.*(1+ f);
%plot(x,mse,x,mse2)
%return;
%sqrt(sum((mse-mse2).^2))

if ~dfn_input
 dx = 2/N;
 lam = (1-[-1:dx:1]'.^2)/dx^2;
 lam = D.*lam;

 M = zeros(N,1)*nan;
 for i=2:N-1
%  M(i) = ( lam(i+1)*T(i+1)-lam(i+1)*T(i)...
%           -lam(i)*T(i)+lam(i)*T(i-1) );
  M(i) = ( lam(i+1)*mse(i+1)-lam(i+1)*mse(i)...
           -lam(i)*mse(i)+lam(i)*mse(i-1) );
 end
% M(1) = lam(2)*(T(2)-T(1));
% M(end) = lam(end-1)*(T(end-1)-T(end));
 M(1) = lam(2)*(mse(2)-mse(1));
 M(end) = lam(end-1)*(mse(end-1)-mse(end));
 dfn = M;

end

F = ( dfn - A - B*T + SW )/c;
