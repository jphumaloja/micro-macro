function [A,B] = nmFD(lamy,muy,sigy,Wy,thy,psiy,qy,n,m,xg)
%NMFD Computes a finite difference approximation \dot z = Az+BU with xg 
%spatial grid points for an %n+m system based on given continuum parameters

yy = linspace(1/n,1,n); % grid for y
hh = linspace(1/m,1,m); % grid for eta
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % x grids for v and u, respectively
% n+m parameters (sampled from continuum on the y and eta grids)
lam = @(i,x) lamy(x,yy(i));
mu = @(j,x) muy(x,hh(j));
sig = @(i,j,x) sigy(x,yy(i),yy(j));
W = @(i,j,x) Wy(x,yy(i),hh(j));
th = @(i,j,x) thy(x,hh(i),yy(j));
psi = @(i,j,x) psiy(x,hh(i),hh(j));
q = @(i,j) qy(yy(i),hh(j));

% finite difference approximation into \dot z = Az + BU
A = zeros((n+m)*xg); B = zeros((n+m)*xg, m);
for k = 1:m
  B((n+k)*xg,k) = mu(k,1)*xg; % control operator v1(1) = U
end

% D is backward difference, -D' is forward difference
D = eye(xg) - diag(ones(1,xg-1), -1);
% fill A based on finite differences
for k = 1:n
  Ik = (k-1)*xg+1:k*xg; % index set
  A(Ik,Ik) = -xg*lam(k,x1).*D; % tranport term
  for l = 1:n
    Il = (l-1)*xg+1:l*xg; % index set
    A(Ik,Il) = A(Ik,Il) + diag(sig(k,l,x1))/n; % sigma terms (scale 1/n)
  end
  for l = n+1:n+m
    % boundary condition (scale 1/m)
    A(Ik(1),(l-1)*xg+1) = A(Ik(1),(l-1)*xg+1) + xg*lam(k,0)*q(k,l-n)/m;
    Il = (l-1)*xg+1:l*xg; % index set
    A(Ik,Il) = A(Ik,Il) + diag(W(k,l-n,x0))/m; % W terms (scale 1/m)
  end
end
for k = n+1:n+m
  Ik = (k-1)*xg+1:k*xg; % index set
  A(Ik,Ik) = -xg*mu(k-n,x0).*D'; % transport term
  for l = 1:n
    Il = (l-1)*xg+1:l*xg; % index set
    A(Ik,Il) = A(Ik,Il) + diag(th(k-n,l,x1))/n; % theta1 terms (scale 1/n)
  end
  for l = n+1:n+m
    Il = (l-1)*xg+1:l*xg;
    A(Ik,Il) = A(Ik,Il) + diag(psi(k-n,l-n,x1))/m; % psi term (scale 1/m)
  end
end
end