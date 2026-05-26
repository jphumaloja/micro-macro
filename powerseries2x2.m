function [Ks,Ls,rd] = powerseries2x2(Lf,Mf,Sf,Wf,THf,Q,N)
%POWERSERIES2X2 Solves 2x2 kernels kernel equations for given parameters
%using power series approximation of order N (analogous to powerseries4D).

syms x z
% compute Taylor series for parameters
L = taylor(Lf,x,'order',N+1);
M = taylor(Mf,x,'order',N+1);
W = taylor(Wf,x,'order',N+1);
S = taylor(Sf,x,'order',N+1);
TH = taylor(THf,x,'order',N+1);

% initalizations
NK = (N+1)*(N+2)/2; % number of terms in power series
KC = sym('k', [1 4*NK],'real'); % symbols for coefficients
rind = 0;
Ks = 0; Ls = 0;
for tm=0:N
  for tn=0:tm
    rind = rind + 1;
    Ks = Ks + KC(rind)*x^tn*z^(tm-tn);
    Ls = Ls + KC(NK+rind)*x^tn*z^(tm-tn);
  end
end

% construct equations and boundary conditions with power series
KE1 = M*diff(Ks,x) - diff(Ks,z)*subs(L,x,z) ...
  - Ks*subs(diff(L,x),x,z) - Ks*subs(S,x,z) - Ls*subs(TH,x,z);
KE2 = M*diff(Ls,x) + diff(Ls,z)*subs(M,x,z) ...
  + Ls*subs(diff(M,x),x,z) - Ks*subs(W,x,z);
BC1 = (L+M)*subs(Ks,z,x) + TH;
BC2 = M*subs(Ls,z,0) - Q*L*subs(Ks,z,0);
% extract coefficients for different powers of spatial variables
KE1C = coeffs(KE1,[x z]);
KE2C = coeffs(KE2,[x z]);
BC1C = coeffs(BC1,x);
BC2C = coeffs(BC2,x);
% construct set linear equations for the coefficients based on the above%
neqs = [numel(KE1C), numel(KE2C), numel(BC1C), numel(BC2C)];
A = zeros(sum(neqs), 2*NK); % initialize A...
B = zeros(sum(neqs),1); % ...and b
EQS = {KE1C,KE2C,BC1C,BC2C}; % cell array for equations
rind = 0; % row index
for ii = 1:4
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
end
% solve equations and compute residual
Kc = A\B;
rd = norm(A*Kc-B);
% insert the obtained values into the power series and compute gains
Ks = subs(Ks,KC(1:NK),Kc(1:NK)');
Ls = subs(Ls,KC(NK+1:2*NK),Kc(NK+1:2*NK)');
Ks = subs(Ks,x,1);
Ls = subs(Ls,x,1);
end