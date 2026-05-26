function [Km,Lm] = controlgain(muy,n,m,xg,K0s,K1s,K2s,L0s,L1s,L2s)
%CONTROLGAIN Samples n+m control gains based on 4-D continuum kernels
%   Detailed explanation goes here

syms x z y h
% same grids as in the finite difference approximation of the n+m system
yy = linspace(1/n,1,n); % grid for y
hh = linspace(1/m,1,m); % grid for eta
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % x grids for u and v

Km = zeros(m,xg*n); Lm = zeros(m,xg*m);
% sample control gains pointwise from the respective continuum kernels 
% depending on the segment P0, P1, or P2
for l = 1:m
  for k = 1:n
    if hh(l) <= yy(k)
      xh1 = x1(x1 <= muy(1,yy(k))/muy(1,hh(l)));
      x1h = x1(x1 > muy(1,yy(k))/muy(1,hh(l)));
      Km(l,(k-1)*xg+1:(k-1)*xg+numel(xh1)) = ...
        subs(subs(K0s,[x h y], [1 hh(l) yy(k)]),z,xh1);
      Km(l,k*xg-numel(x1h)+1:k*xg) = ...
        subs(subs(K1s,[x h y], [1 hh(l) yy(k)]),z,x1h);
    else
      Km(l,(k-1)*xg+1:k*xg) = ...
        subs(subs(K2s,[x h y], [1 hh(l) yy(k)]),z,x1);
    end
  end
  for k = 1:m
    if hh(l) <= hh(k)
      xh0 = x0(x0 <= muy(1,hh(k))/muy(1,hh(l)));
      x0h = x0(x0 > muy(1,hh(k))/muy(1,hh(l)));
      Lm(l,(k-1)*xg+1:(k-1)*xg+numel(xh0)) = ...
        subs(subs(L0s,[x h y], [1 hh(l) hh(k)]),z,xh0);
      Lm(l,k*xg-numel(x0h)+1:k*xg) = ...
        subs(subs(L1s,[x h y], [1 hh(l) hh(k)]),z,x0h);      
    else
      Lm(l,(k-1)*xg+1:k*xg) = ...
        subs(subs(L2s,[x h y], [1 hh(l) hh(k)]),z,x0);
    end
  end
end
end