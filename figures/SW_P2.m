function varargout=SW_T(x,albo,albi)
%
% albo is coabledo globalmean e.g., 0.7
% albi is coalbedo P2 e.g., -0.2

  error(nargchk(3, 3, nargin))

  if albo < 0.5
   warning('Use co-albedo values')
  end

  Q=1360/4;
  P2=1/2*(3*x.^2-1);
  S=1-0.482*P2;
  coalb=albo*ones(size(x))+albi*P2;
  SW=Q*S.*coalb;      

  if nargout > 0
    varargout{1} = SW;
    varargout{2} = 1-coalb;
  end
