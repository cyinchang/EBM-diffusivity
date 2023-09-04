function D = calc_D(D0,T,Tref,x,do_T0,do_T2,do_h2,gamma,rh)
% input variables:
% D0 in W/m2/K
% T(x) and Tref(x) in K or C
% x is sin(lat)  
% do_T0 -> D = D0( 1 + gamma (T0-T0ref) )
% do_T2 -> D = D0 * (T2/T2ref)^do_T2
% do_h2 -> D = D0 * (h2/h2ref)^do_h2
%   there could be coefficients here like beta^-2
%   following Chang and Held 2022, but that would also have a T2^3 dependence
%   vs. the analytic theory May 1 (3)  
%   could have an adhoc gamma 
%
% output D in W/m2/K

  error(nargchk(4, 9, nargin))

  if nargin == 4
    do_T0 = 0;
    do_T2 = 0;
    do_h2 = 0;
    gamma = 0.0;
    rh = 0.0;
  elseif nargin == 5;
    do_T2 = 0;
    do_h2 = 0;
    gamma = 0.0;
    rh = 0.0;
  elseif nargin > 5;
    if do_T0 == 1 & do_T2 >= 1
      error('cannot have both do_T0 and do_T2 active');
    end
    if do_T0 == 1 & do_h2 >= 1
      error('cannot have both do_T0 and do_h2 active');
    end
  end

  % for do_T0 == 0, do_T2 == 0, and do_h2 == 0, return D = D0
  D = D0;
  
  if do_T0
    % hard coded gamma, should be input argument
    %    gamma = -0.03;
    %gamma = -0.0275;
    % assume uniform grid in x for global mean
    T0 = mean(T);
    T0ref = mean(Tref);
    D = D0*( 1 + gamma*(T0-T0ref) );
  end

  if do_T2 >= 1   
    lat = asin(x)*180/pi;   phi_rad = lat*pi/180;
    Pmax = 4; do_cos_weight = 0;
    
%     [LPcoef, LPvecall, field_trunc] = ...
%       fit_legendreP_coef(T, phi_rad, Pmax, do_cos_weight);
%     T2 = LPcoef(3);
%     [LPcoef, LPvecall, field_trunc] = ...
%       fit_legendreP_coef(Tref, phi_rad, Pmax, do_cos_weight);
%     T2ref = LPcoef(3);

    [T0,T2]=compute_T0T2(x, T);
    [T0ref,T2ref]=compute_T0T2(x, Tref);

    D = D*(T2/T2ref)^do_T2;
  end

  if do_h2 >= 1   
    lat = asin(x)*180/pi;   phi_rad = lat*pi/180;
    Pmax = 4; do_cos_weight = 0;
    
    h = calc_mse(T,rh);
    href = calc_mse(Tref,rh);

%     [LPcoef, LPvecall, field_trunc] = ...
%       fit_legendreP_coef(h, phi_rad, Pmax, do_cos_weight);
%     h2 = LPcoef(3);
%     [LPcoef, LPvecall, field_trunc] = ...
%       fit_legendreP_coef(href, phi_rad, Pmax, do_cos_weight);
%     h2ref = LPcoef(3);

    [h0,h2]=compute_T0T2(x, h);
    [h0ref,h2ref]=compute_T0T2(x, href);

    D = D*(h2/h2ref)^do_h2;
  end
  
