% parameters (symbolic for solving kernel equations with power series)
syms x z y h
Lf = 1;
Mf = 2-y;
Sf = (x+1)*y*(h+1/2);
Wf = (x+1)*y*(h+1/2);
THf = (x+1)*y*(h+1/2);
PSf = y-h;
Qf = (y+1/2)*h;

% solve 4-D continuum kernel equations for the parameters
N = 10; % approximation order for power series
msg = true; % display information
[K0s,K1s,K2s,L0s,L1s,L2s,rd] = powerseries4D(Lf,Mf,Sf,Wf,THf,PSf,Qf,N,msg);

% write parameters as function handles
lamy = @(x,y) ones(size(x));
muy = @(x,y) (2-y)*ones(size(x));
sigy = @(x,y,h) (x+1).*(y).*(h+0.5);
Wy = @(x,y,h) (x+1).*(y).*(h+0.5);
thy = @(x,y,h) (x+1).*(y).*(h+0.5);
psiy = @(x,y,h) (y-h).*ones(size(x));
qy = @(y,h) (y+0.5)*h;

% Illustration of Theorem 3.1
% n+m system based on continuum parameters for selected n and m
n = 50; yy = linspace(1/n,1,n); % grid for y
m = 50; hh = linspace(1/m,1,m); % grid for eta
xg = 128; % grid for x (minus one end point determined by boundary conditions)
% finite-difference approximation for the n+m system
[A,B] = nmFD(lamy,muy,sigy,Wy,thy,psiy,qy,n,m,xg);
% sample control gains for the n+m system from the 4-D continuum kernels
[Km,Lm] = controlgain(muy,n,m,xg,K0s,K1s,K2s,L0s,L1s,L2s);
% simulate
opts = odeset('Jacobian', A+B*[Km/n Lm/m]/xg);
T = 5; % simulation end time
% initial condition
z0 = [kron(trapz(hh,qy(yy',hh),2),ones(xg,1)); ones(m*xg,1)];
% simulate, approximate integral over z in control law by endpoint rule
sol = ode45(@(t,z) (A + B*[Km/n Lm/m]/xg)*z, [0, T], z0, opts);
% evaluate solution at desired time grid
tg = 513; % number of time grid points
TT = linspace(0,T,tg); % time grid
usol = deval(sol, TT); % evaluate solution...
Usol = [Km/n Lm/m]/xg*usol; % ...and controls on the time grid
% spatial grid for plotting the PDE solution (not used)
% tmp = linspace(0,1,xg+1); % auxiliary grid
% x0 = tmp(1:xg); x1 = tmp(2:xg+1); % x spatial grids for v and u, respectively
% plot
figure(1)
surf(hh,TT,Usol'), shading interp, view([80 10]);
xlabel('$\eta$','interpreter','latex','fontsize',13)
ylabel('$t$','interpreter','latex','fontsize',13)
zlabel('$U(t,\eta)$','interpreter','latex','fontsize',13,'rotation',0)
set(gca,'xtick',[0 1],'ytick',0:5,'ztick',-1.5:.5:0,'fontsize',12)
zlim([-1.5 .2]), clim([-1.6 1.6]), colormap gray

% Illustration of Theorem 4.1
nmd = [2, 5:5:25]; % considered values of n=m
ndd = zeros(numel(nmd),tg); % for storing state norms
for k = 1:numel(nmd)
% n+m system based on continuum parameters for selected n and m
n = nmd(k); yy = linspace(1/n,1,n); % grid for y
m = nmd(k); hh = linspace(1/m,1,m); % grid for eta
% finite-difference approximation for the n+m system
[A,B] = nmFD(lamy,muy,sigy,Wy,thy,psiy,qy,n,m,xg);
% sample control gains for the n+m system from the 4-D continuum kernels
[Km,Lm] = controlgain(muy,n,m,xg,K0s,K1s,K2s,L0s,L1s,L2s);
% simulate
opts = odeset('Jacobian', A+B*[Km/n Lm/m]/xg);
T = 5; % simulation end time
% initial condition
z0 = [kron(trapz(hh,qy(yy',hh),2),ones(xg,1)); ones(m*xg,1)];
% simulate, approximate integral over z in control law by endpoint rule
sol = ode45(@(t,z) (A + B*[Km/n Lm/m]/xg)*z, [0, T], z0, opts);
usol = deval(sol, TT); % evaluate solution at time grid
% compute norms of the solution components and store data
nu = sqrt(sum(usol(1:n*xg,:).^2/(n*xg),1));
nv = sqrt(sum(usol(n*xg+1:(n+m)*xg,:).^2/(m*xg),1));
ndd(k,:) = nu+nv;
if msg
  disp(num2str(n));
end
end
% plot
figure(2)
plot(TT, ndd,'linewidth',2)
xlabel('$t$','interpreter','latex','fontsize',13)
ylabel('$\|(\mathbf{u}(t),\mathbf{v}(t))\|_E$','interpreter','latex',...
  'fontsize',13)
legend({'$n=m=2$','$n=m=5$','$n=m=10$','$n=m=15$','$n=m=20$','$n=m=25$'},...
  'interpreter','latex','fontsize',13)
set(gca,'tickdir','out','fontsize',12,'ytick',0:2,'xtick',0:5)
ylim([0 2.8])