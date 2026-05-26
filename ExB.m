% parameters (symbolic for solving kernel equations with power series)
syms x z y h
Lf = 1;
Mf = (2-y)/2;
Sf = x*y*(h+1/2);
Wf = x*y*(h+1/2);
THf = x*y*(h+1/2);
PSf = (y-h)/2;
Qf = 1;

% solve 4-D continuum kernel equations for the parameters
N = 10; % approximation order for power series
msg = true; % display information
[K0s,K1s,K2s,L0s,L1s,L2s,rd] = powerseries4D(Lf,Mf,Sf,Wf,THf,PSf,Qf,N,msg);

% write parameters as function handles
lamy = @(x,y) ones(size(x));
muy = @(x,y) (2-y)/2*ones(size(x));
sigy = @(x,y,h) x.*(y).*(h+0.5);
Wy = @(x,y,h) x.*(y).*(h+0.5);
thy = @(x,y,h) x.*(y).*(h+0.5);
qy = @(y,h) ones(size(y)).*ones(size(h));
psiy = @(x,y,h) (y-h).*ones(size(x))/2;

% n+m system based on continuum parameters for selected n and m
n = 10; yy = linspace(1/n,1,n); % grid for y
m = 10; hh = linspace(1/m,1,m); % grid for eta
xg = 128; % grid for x (minus one end point determined by boundary conditions)
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % spatial grids for v and u
% finite-difference approximation for the n+m system
[A,B] = nmFD(lamy,muy,sigy,Wy,thy,psiy,qy,n,m,xg);
% sample control gains for the n+m system from the 4-D continuum kernels
[Km,Lm] = controlgain(muy,n,m,xg,K0s,K1s,K2s,L0s,L1s,L2s);
Kc1 = [Km/n, Lm/m]/xg; % control gain in finite-difference approximation

% continuum system (approximation with double n,m)
nc = 20; yyc = linspace(1/nc,1,nc); % grid for y
mc = 20; hhc = linspace(1/mc,1,mc); % grid for eta
% finite-difference approximation for the continuum system
[Ac,Bc] = nmFD(lamy,muy,sigy,Wy,thy,psiy,qy,nc,mc,xg);
Kc2 = repelem(Kc1,mc/m,1); % control gain for finite differences

% cascade of n+m and continuum systems
Ae = blkdiag(A,Ac);
% auxiliary function for the continuum-based control law
zbar = @(z) [...
  repmat(mean(reshape(z((n+m)*xg+1:(n+m+nc)*xg),xg,nc),2),n,1); ...
  repmat(mean(reshape(z((n+m+nc)*xg+1:(n+m+nc+mc)*xg),xg,mc),2),m,1)];

% average paramters (symbolic for solving kernel equations with power series)
Laf = 1;
Maf = 0.75;
Saf = 0.5*x;
Waf = 0.5*x;
THaf = 0.5*x;
Qa = 1;
% solve 2x2 kernel equations for the average parameters
[Kas,Las,rda] = powerseries2x2(Laf,Maf,Saf,Waf,THaf,Qa,N);
% vectorized average parameters
la = ones(1,xg);
ma = 0.75*ones(1,xg);
sa = x1*0.5;
Wa = x0*0.5;
tha = x1*0.5;
qa = 1;
% finite difference approximation for the average part
D = eye(xg) - diag(ones(1,xg-1), -1); % backward difference
Aa = [-xg*la.*D + diag(sa), diag(Wa); diag(tha), -xg*ma.*D'];
Aa(1,xg+1) = Aa(1,xg+1) + xg*la(1)*qa;
Ba = [zeros(2*xg-1,1); xg*ma(xg)];
% control gain for finite differences
Ka = double([subs(Kas,z,x1), subs(Las,z,x0)])/xg;

% cascade of n+m and average systems
Aea = [A, B*repmat(Ka,m,1); zeros(2*xg,(n+m)*xg), Aa + Ba*Ka];

% simulate n+m system under continuum-based controls
opts = odeset('Jacobian', Ae);
T = 7; % simulation end time
z0 = [kron(trapz(hh,qy(yy',hh),2),ones(xg,1)); ones(m*xg,1); ...
  ones(nc*xg,1); ones(mc*xg,1)]; % initial condition
% solve finite-difference ODE
sol = ode45(@(t,z) Ae*z + [B*Kc1; Bc*Kc2]*zbar(z), [0, T], z0, opts);
% evaluate solution at desired time grid
tg = 513; % number of time grid points
TT = linspace(0,T,tg); % time grid
usol = deval(sol, TT); % evaluate solution at time grid
% compute control inputs
Usolc = zeros(m,tg);
for k = 1:tg
  tmp = Kc2*zbar(usol(:,k));
  Usolc(:,k) = tmp(1:mc/m:nc-1);
end
% compute norms of the solution components
nu = sqrt(sum(usol(1:n*xg,:).^2/(n*xg),1));
nv = sqrt(sum(usol(n*xg+1:(n+m)*xg,:).^2/(m*xg),1));
nmc = nu+nv; % state norm under continuum-based controls

% simulate n+m system under average-based controls
opts = odeset('Jacobian', Aea); % should expedite computations
z0 = [kron(trapz(hh,qy(yy',hh),2), ones(xg,1)); ones(m*xg,1); ...
  1*ones(xg,1); ones(xg,1)]; % initial condition
sol = ode45(@(t,z) Aea*z, [0, T], z0, opts);
usol = deval(sol, TT); % evaluate solution...
Usola = [zeros(1,(n+m)*xg), Ka]*usol; % ...and controls at tg
% compute norms of the solution components
nu = sqrt(sum(usol(1:n*xg,:).^2/(n*xg),1));
nv = sqrt(sum(usol(n*xg+1:(n+m)*xg,:).^2/(m*xg),1));
nma = nu+nv;

% simulate open-loop system
opts = odeset('Jacobian', A);
z0 = [kron(trapz(hh,qy(yy',hh),2),ones(xg,1)); ones(m*xg,1)]; % init.
% solve finite-difference ODE
sol = ode45(@(t,z) A*z, [0, T], z0, opts);
usol = deval(sol, TT); % evaluate solution at time grid
% compute norms of the solution components
nu = sqrt(sum(usol(1:n*xg,:).^2/(n*xg),1));
nv = sqrt(sum(usol(n*xg+1:(n+m)*xg,:).^2/(m*xg),1));
nmo = nu+nv; % open-loop state norm

% plot
figure(3)
RGB = orderedcolors("gem");
subplot(311)
plot(TT,nmc,'linewidth',2,'color',RGB(1,:));
hold on
plot(TT,nmc,'linestyle','none')
plot(TT,nma,'linewidth',2,'color',RGB(2,:))
plot(TT,nmo,'linewidth',2,'color',RGB(3,:))
hold off
ylabel('$\|(\mathbf{u}(t),\mathbf{v}(t))\|_E$','interpreter','latex',...
  'fontsize',13)
legend({'continuum',' ','average','autonomous',},'interpreter','latex',...
  'fontsize',13,'numcolumns', 2,'box','off','location','northeast')
set(gca,'tickdir','out','fontsize',12,'xticklabel',{},'ytick',0:2,...
  'position',get(gca,'position')+[0.01 0.01 0.06 0.02])
ylim([0,2.2])
subplot(312)
plot(TT, Usolc, 'linewidth',2,'color',RGB(1,:))
ylabel('$\mathbf{U}(t)$','interpreter','latex','rotation',0,'fontsize',13)
set(gca,'tickdir','out','fontsize',12, 'xticklabel',{},'ytick',-.5:.5:0,...
  'position',get(gca,'position')+[0.01 0.02 0.06 0.02])
ylim([-.7, 0.1])
legend('continuum','interpreter','latex',...
  'fontsize',13,'location','southeast','box','off')
subplot(313)
plot(TT, Usola, 'linewidth',2,'color',RGB(2,:))
ylabel('$U(t)$','interpreter','latex','rotation',0,'fontsize',13)
xlabel('$t$','interpreter','latex','fontsize',13)
set(gca,'tickdir','out','fontsize',12, 'ytick',-.4:.1:0,... 
  'position',get(gca,'position')+[0.01 0.03 0.06 0.02])
ylim([-0.43, 0.01])
legend('average','interpreter','latex',...
  'fontsize',13,'location','southeast','box','off')