function mse = calc_mse(T,rh,do_linCC,Tref)
% input variables:
% T in K or C
% rh in fraction 0-1, assumed scalar not array
% do_linCC 1/0, NB always assumes ps = 1e5 Pa
% Tref for linear CC; same size as T
%
% output mse in K

  error(nargchk(2, 4, nargin))

  if nargin == 2
    do_linCC = 0;
  elseif nargin == 3;
    error('need Tref input for linearized CC');
  end

  % make sure T is in K
  if max(T)<100
    T = T+273.15;
  end

  [tmp,qs]=saturation_thermodynamics(T,1e5);

  if do_linCC
    if max(Tref)<100
      Tref = Tref+273.15;
    end
    [tmp,qs0]=saturation_thermodynamics(Tref,1e5);
    [tmp, q2] = saturation_thermodynamics(Tref+0.1,1e5);
    [tmp, q1] = saturation_thermodynamics(Tref-0.1,1e5);
    dqsdTref = (q2 - q1)/0.2;
    qs = dqsdTref.*(T-Tref) + qs0;
  end

  mse = T + PARS('l_cond')*qs*rh/PARS('cp');
