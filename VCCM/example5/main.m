%% Example 5: require YALMIP, SDPT3
%  
%  A scalar nonlinear dynamics: dx/dt = -x+x^3+u, z=(x, Wu*u)
%
function main


%% Control synthesis via Eq.(57) with alpha=0.5 and Wu=0.1
Wu=0.1; alpha=0.5; 
A=@(rho) -(1-rho); B=1; Bw=0; C=[1;0]; D=[0;Wu]; Dw=[0;0];
nx=1; nw=size(Bw,2); nz=size(Dw,1); nH=nx+nw+nz; 
xmax=2; rho_b=[0,xmax^2]; eps=1e-10;

W=sdpvar(1); Y0=sdpvar(1); Y1=sdpvar(1); Y=@(rho) Y0+Y1*rho;
H1=@(rho) [A(rho)*W+W*A(rho)'+B*Y(rho)+Y(rho)'*B',Bw,W*C'+Y(rho)'*D';
           Bw', -alpha*eye(nw),Dw';
           C*W+D*Y(rho),Dw,-alpha*eye(nz)];
LMI=[H1(rho_b(1))<=-eps*eye(nH),H1(rho_b(2))<=-eps*eye(nH),W>=eps*eye(nx)];
optimize(LMI,[],sdpsettings('verbose',0,'solver','sdpt3'));
Y0=value(Y0); Y1=value(Y1); W=value(W);

K0=Y0/W; K1=Y1/W;

save('example5.mat','K0','K1','xmax','Wu');
% load('example5.mat');

% initial state and simulation time
x0=2; T=2;

lw=1; fz=10; lz=8;

%% set-point 1
xe=-1.9; ue=xe-xe^3;
[glpv.t,glpv.x] = ode45(@(t,x) CL_GLPV(t,x,K0,K1,xe,ue),[0 T],x0);
[vccm.t,vccm.x] = ode45(@(t,x) CL_VCCM(t,x,K0,K1,xe,ue),[0 T],x0);

subplot(131); hold on
plot(glpv.t,glpv.x(:,1),'r','linewidth',lw); 
plot(vccm.t,vccm.x(:,1),'b--','linewidth',lw);
xlabel('t');
ylabel('x');
set(gca,'fontsize',fz);
legend({'GLPV','VCCM'},'fontsize',lz);

%% set-point 2
xe=0.0; ue=xe-xe^3;
[glpv.t,glpv.x] = ode45(@(t,x) CL_GLPV(t,x,K0,K1,xe,ue),[0 T],x0);
[vccm.t,vccm.x] = ode45(@(t,x) CL_VCCM(t,x,K0,K1,xe,ue),[0 T],x0);

subplot(132); hold on
plot(glpv.t,glpv.x(:,1),'r','linewidth',lw); 
plot(vccm.t,vccm.x(:,1),'b--','linewidth',lw);
xlabel('t');
set(gca,'fontsize',fz);
legend({'GLPV','VCCM'},'fontsize',lz);

%% set-point 3
xe=1.9; ue=xe-xe^3;
[glpv.t,glpv.x] = ode45(@(t,x) CL_GLPV(t,x,K0,K1,xe,ue),[0 T],x0);
[vccm.t,vccm.x] = ode45(@(t,x) CL_VCCM(t,x,K0,K1,xe,ue),[0 T],x0);

subplot(133); hold on
plot(glpv.t,glpv.x(:,1),'r','linewidth',lw); 
plot(vccm.t,vccm.x(:,1),'b--','linewidth',lw);
ylim([1,20]);
xlabel('t');
set(gca,'fontsize',fz);
legend({'GLPV','VCCM'},'fontsize',lz);

savefig('example5.fig');
end


function dxdt = CL_GLPV(t,x,K0,K1,xe,ue)
u=ue+(K0+K1*x(1)^2)*(x(1)-xe);
dxdt=-x(1)+x(1)^3+u;
end

function dxdt = CL_VCCM(t,x,K0,K1,xe,ue)
mu=ue+xe*(xe^2-x(1)^2);
u=mu+(K0+K1*x(1)^2)*(x(1)-xe);
dxdt=-x(1)+x(1)^3+u;
end
