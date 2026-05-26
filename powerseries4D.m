function [K0s,K1s,K2s,L0s,L1s,L2s,rd] = powerseries4D(Lf,Mf,Sf,Wf,THf,PSf,Qf,N,msg)
%POWERSERIES4D Solves 4-D continuum kernel equations for given parameters
%   Output arguments:
%   K0s: kernel K restricted to P0
%   K1s: kernel K restricted to P1
%   K2s: kernel K restricted to P2
%   L0s: kernel L restricted to P0
%   L1s: kernel L restricted to P1
%   L2s: kernel L restricted to P2
%   rd: residual of the least squares fit
%   Input arguments
%   Lf: lambda (symbolic function of x and y)
%   Mf: mu (symbolic function of x and y)
%   Sf: sigma (symbolic function of x, y and h)
%   Wf: W (symbolic function of x, y and h)
%   THf: theta (symbolic function of x, y and h)
%   PSf: psi (symbolic function of x, y and h)
%   Qf: q (symbolic function of y)
%   N: order of the power series approximation
%   msg: true/false to display messages 

syms x z h y s
% compute Taylor series for parameters
L = taylor(Lf,[x y],'order',N+1);
M = taylor(Mf,[x y],'order',N+1);
W = taylor(Wf,[x y h],'order',N+1);
S = taylor(Sf,[x y h],'order',N+1);
TH = taylor(THf,[x y h],'order',N+1);
PS = taylor(PSf,[x y h],'order',N+1);
Q = taylor(Qf,[y h],'order',N+1);
% Taylor series for characteristic surface (assumes constant mu in x)
RHO = x*taylor(Mf/subs(Mf,y,h),[y h],'order',N+1);
% initialize power series for K and L
NK = N*(N+1)*(N^2+9*N+26)/24+N+1; % number of terms for K and L
KC = sym('k', [1 6*NK],'real'); % symbols for coefficients of K and L
rind = 0; % row index
K0 = 0; K1 = 0; K2 = 0; % 3 K kernels for 3 domains
L0 = 0; L1 = 0; L2 = 0; % 3 L kernels for 3 domains
% create power series with unknown coefficients
for ky=0:N
  for kh = 0:N-ky
    for kz = 0:N-ky-kh
      for kx = 0:kz
        rind = rind + 1;
        K0 = K0 + KC(rind)*x^kx*z^(kz-kx)*h^kh*y^ky;
        K1 = K1 + KC(rind+NK)*x^kx*z^(kz-kx)*h^kh*y^ky;
        K2 = K2 + KC(rind+2*NK)*x^kx*z^(kz-kx)*h^kh*y^ky;
        L0 = L0 + KC(rind+3*NK)*x^kx*z^(kz-kx)*h^kh*y^ky;
        L1 = L1 + KC(rind+4*NK)*x^kx*z^(kz-kx)*h^kh*y^ky;
        L2 = L2 + KC(rind+5*NK)*x^kx*z^(kz-kx)*h^kh*y^ky;
      end
    end
  end
end

% construct equations and boundary conditions with power series
KE0 = subs(M,y,h)*diff(K0,x) -  subs(L,x,z)*diff(K0,z) ... 
  - subs(diff(L,x),x,z)*K0 ...
  - int(subs(S,[x,y,h],[z,s,y])*subs(K0,y,s),s,0,1) ...
  - int(subs(TH,[x,y,h],[z,s,y])*subs(L0,y,s),s,0,1);
KE1 = subs(M,y,h)*diff(K1,x) -  subs(L,x,z)*diff(K1,z) ... 
  - subs(diff(L,x),x,z)*K1 ...
  - int(subs(S,[x,y,h],[z,s,y])*subs(K1,y,s),s,0,1) ...
  + int(subs(TH,[x,y,h],[z,s,y])*subs(L1,y,s),s,0,1);
KE2 = subs(M,y,h)*diff(K2,x) -  subs(L,x,z)*diff(K2,z) ... 
  - subs(diff(L,x),x,z)*K2 ...
  - int(subs(S,[x,y,h],[z,s,y])*subs(K2,y,s),s,0,1) ...
  + int(subs(TH,[x,y,h],[z,s,y])*subs(L2,y,s),s,0,1);
LE0 = subs(M,y,h)*diff(L0,x) + subs(M,x,z)*diff(L0,z) ...
  + subs(diff(M,x),x,z)*L0 ...
  - int(subs(W,[x,y,h],[z,s,y])*subs(K0,y,s),s,0,1) ...
  - int(subs(PS,[x,y,h],[z,s,y])*subs(L0,y,s),s,0,1);
LE1 = subs(M,y,h)*diff(L1,x) + subs(M,x,z)*diff(L1,z) ...
  + subs(diff(M,x),x,z)*L1 ...
  - int(subs(W,[x,y,h],[z,s,y])*subs(K1,y,s),s,0,1) ...
  - int(subs(PS,[x,y,h],[z,s,y])*subs(L1,y,s),s,0,1);
LE2 = subs(M,y,h)*diff(L2,x) + subs(M,x,z)*diff(L2,z) ...
  + subs(diff(M,x),x,z)*L2 ...
  - int(subs(W,[x,y,h],[z,s,y])*subs(K2,y,s),s,0,1) ...
  - int(subs(PS,[x,y,h],[z,s,y])*subs(L2,y,s),s,0,1);
BCK1 = (subs(M,y,h)+L)*subs(K1,z,x) + subs(TH,[y,h],[h,y]);
BCK2 = (subs(M,y,h)+L)*subs(K2,z,x) + subs(TH,[y,h],[h,y]);
BCL0 = subs(M,x,0)*subs(L0,z,0) - ...
  int(subs(K0,[z,y],[0,s])*subs(L,[x,y],[0,s])*subs(Q,[y,h],[s,y]),s,0,1);
BCL1 = (subs(M,y,h)-M)*subs(L1,z,x) + subs(PS,[y,h],[h,y]);
BCL2 = (subs(M,y,h)-M)*subs(L2,z,x) + subs(PS,[y,h],[h,y]);
BCLE = (subs(M,[x,y],[z,h]) - subs(M,x,z))*subs(L2,x,1) ...
  + subs(PS,[x,y,h],[z,h,y]);
BCKC = subs(K0,z,RHO) - subs(K1,z,RHO);
% extract coefficients for different powers of spatial variables
KE0C = coeffs(KE0,[x z h y]);
KE1C = coeffs(KE1,[x z h y]);
KE2C = coeffs(KE2,[x z h y]);
LE0C = coeffs(LE0,[x z h y]);
LE1C = coeffs(LE1,[x z h y]);
LE2C = coeffs(LE2,[x z h y]);
BK1C = coeffs(BCK1,[x h y]);
BK2C = coeffs(BCK2,[x h y]);
BL0C = coeffs(BCL0,[x h y]);
BL1C = coeffs(BCL1,[x h y]);
BL2C = coeffs(BCL2,[x h y]);
BLEC = coeffs(BCLE,[z h y]);
BKCC = coeffs(BCKC,[x h y]);
% construct set linear equations for the coefficients based on the above
neqs = [numel(KE0C), numel(KE1C), numel(KE2C), ...
  numel(LE0C), numel(LE1C), numel(LE2C), ...
  numel(BK1C), numel(BK2C), numel(BL0C), ...
  numel(BL1C), numel(BL2C), numel(BLEC), numel(BKCC)]; 
A = zeros(sum(neqs), 6*NK); % initialize A...
B = zeros(sum(neqs),1); % ...and b
EQS = {KE0C,KE1C,KE2C,LE0C,LE1C,LE2C,BK1C,BK2C,BL0C,BL1C,BL2C,BLEC,BKCC};
if msg
  disp(['Parsing data: ',num2str(sum(neqs)),' equations for ',...
    num2str(6*NK),' unknonws']);
end
rind = 0; % row index
for ii = 1:13 % loop over all equations and boundary conditions
  DAT = EQS{ii};
  for m=1:neqs(ii)
    rind = rind + 1;
    [cc, kc] = coeffs(DAT(m),KC);
    for mk=1:numel(kc) % check which K's are present in the coefficients
      str = char(kc(mk)); 
      sl = numel(str);
      if sl > 1 % look up the index and insert to A
        A(rind,str2double(str(2:sl))) = cc(mk);
      else % if no index, it's a contant; insert to B
        B(rind) = -cc(mk);
      end
    end
  end
  if msg % display number when an equation or a boundary condition is done
    disp(num2str(ii))
  end
end
% solve equations and compute residual
Kc = A\B;
rd = norm(A*Kc-B);
% insert the obtained coefficients into the power series
K0s = subs(K0,KC(1:NK),Kc(1:NK)');
K1s = subs(K1,KC(NK+1:2*NK),Kc(NK+1:2*NK)');
K2s = subs(K2,KC(2*NK+1:3*NK),Kc(2*NK+1:3*NK)');
L0s = subs(L0,KC(3*NK+1:4*NK),Kc(3*NK+1:4*NK)');
L1s = subs(L1,KC(4*NK+1:5*NK),Kc(4*NK+1:5*NK)');
L2s = subs(L2,KC(5*NK+1:6*NK),Kc(5*NK+1:6*NK)');
if msg % display number when an equation or a boundary condition is done
  disp(['Kernels solved with residual ',num2str(rd)])
end
end