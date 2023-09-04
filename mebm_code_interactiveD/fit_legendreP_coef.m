function [LPcoef, LPvecall, field_trunc] = ...
      fit_legendreP_coef(field, phi_rad, Pmax, do_cos_weight)
% output
% LPcoef: coeficients of Legendre polynomials
% LPvecall: matrix of Legendre polynomial vectors
% field_trunc: representation of field truncated at Pmax
%
% input
% field: to fit
% phi_rad: latitude in radians
% Pmax highest Legendre polynomial to fit
%
% This works for phi_rad uniformly spaced.
% Caveat emptor for other grid spacings.

 if nargin < 4
    do_cos_weight = 1;
 end

 x = sin(phi_rad);
 field_trunc = zeros(size(field));

 for p=1:Pmax+1
   LPvec = legendreP(p-1,x);

   LPvecall(p,:) = LPvec;

   % area weighting of these "integrals"
   if do_cos_weight
     LPcoef(p) = LPvec' * ((field-field_trunc).*cos(phi_rad)) ...
	 / (LPvec' * (LPvec.*cos(phi_rad)));
   else
   % leaving out cos weighting seems desirable if 
   % field uniformly spaced in x not phi
     LPcoef(p) = LPvec' * ((field-field_trunc)) ...
	 / (LPvec' * (LPvec));
   end

   field_trunc = LPcoef(p)*LPvec + field_trunc;
 end

