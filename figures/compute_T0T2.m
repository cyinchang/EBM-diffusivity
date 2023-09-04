% compute_T0T2.m

function [T0, T2] = compute_T0T2(x, T)

T0=mean(T);
P2=(3*x.*x-1)/2;
T2=sum(T.*P2)/sum(P2.*P2);

end