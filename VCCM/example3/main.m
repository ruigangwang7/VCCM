%% Example 3: require YALMIP, MOSEK, chebfun
%  dynamics: dx1/dt=-x1+x3, dx2/dt=x1^2-x2-2x1x3+x3, dx3/dt=-x2+u
%  Comparison of LQR, CCM, VCCM controllers
% 
function main

% colors
c=[0 0.4470 0.7410];
p=[0.4940 0.1840 0.5560];
r=[0.8500 0.3250 0.0980];

x0=[0.5;0.5;0.5];
dlqr=run_lqr(x0);
dccm=run_ccm(x0);
dvccm=run_vccm(x0);
savemat('small-x0.mat','dlqr','dccm','dvccm');
% load('small-x0.mat');
figure(1); 
subplot(121); hold on
plot(dvccm.tm,dvccm.x(1,:),'color',c); 
plot(dvccm.tm,dvccm.x(2,:),'color',p);
plot(dvccm.tm,dvccm.x(3,:),'color',r); 
plot(dccm.tm,dccm.x(1,:),'--','color',c);
plot(dccm.tm,dccm.x(2,:),'--','color',p);
plot(dccm.tm,dccm.x(3,:),'--','color',r);
plot(dlqr.tm,dlqr.x(1,:),'-.','color',c);
plot(dlqr.tm,dlqr.x(2,:),'-.','color',p);
plot(dlqr.tm,dlqr.x(3,:),'-.','color',r);
legend('x1 VCCM','x2 VCCM','x3 VCCM','x1 CCM','x2 CCM','x3 CCM','x1 LQR','x2 LQR','x3 LQR');
xlim([0,5]);
xlabel('t'); 

x0=[10;10;10];
dlqr=run_lqr(x0);
dccm=run_ccm(x0);
dvccm=run_vccm(x0);
savemat('large-x0.mat','dlqr','dccm','dvccm');
% load('large-x0.mat');
subplot(122); hold on
plot(dvccm.tm,dvccm.x(1,:),'color',c); 
plot(dvccm.tm,dvccm.x(2,:),'color',p);
plot(dvccm.tm,dvccm.x(3,:),'color',r); 
plot(dccm.tm,dccm.x(1,:),'--','color',c);
plot(dccm.tm,dccm.x(2,:),'--','color',p);
plot(dccm.tm,dccm.x(3,:),'--','color',r);
plot(dlqr.tm,dlqr.x(1,:),'-.','color',c);
plot(dlqr.tm,dlqr.x(2,:),'-.','color',p);
plot(dlqr.tm,dlqr.x(3,:),'-.','color',r);
legend('x1 VCCM','x2 VCCM','x3 VCCM','x1 CCM','x2 CCM','x3 CCM','x1 LQR','x2 LQR','x3 LQR');
xlim([0,8]);
xlabel('t');
ylim([-10,40]);

savefig('example3.fig')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% auxilary functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear dynamics
function dx=dynamics(t,x,u)
dx=[-x(1)+x(3);
     x(1)^2-x(2)-2*x(1)*x(3)+x(3);
    -x(2)+u];
end

% RK45 method
function x1=next_state(tau,model,x0,u0)
k1=tau*model(0,x0,u0);
k2=tau*model(tau/2,x0+k1/2,u0);
k3=tau*model(tau/2,x0+k2/2,u0);
k4=tau*model(tau,x0+k3,u0);
x1=x0+(k1+2*k2+2*k3+k4)/6;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LQR controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = run_lqr(x0)

% LQR design
A=[-1,0,1; 0,-1,1; 0,-1,0]; B=[0; 0; 1]; Q=eye(3); R=1;
[K,P,e]=lqr(A,B,Q,R);

% CL simulation
T=8; ts=0.01; tm=0:ts:T; L=length(tm); 
data.tm=tm; data.x=zeros(3,L); data.u=zeros(1,L); 

xr=[0;0;0]; ur=0; x=x0;
for t=1:L
    data.x(:,t)=x;
    u=ur-K*(x-xr);
    data.u(t)=u;
    x=next_state(ts,@dynamics,x,u);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VCCM controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = run_vccm(x0)

% Same LQR design
A=[-1,0,1; 0,-1,1; 0,-1,0]; B=[0; 0; 1]; Q=eye(3); R=1;
[K,P,e]=lqr(A,B,Q,R);

T=8; ts=0.01; tm=0:ts:T; L=length(tm); 
data.tm=tm; data.x=zeros(3,L); data.u=zeros(1,L); 

xr=[0;0;0]; ur=0; 
x=x0; y=[x; 0];

for t=1:L
    x=y(1:3);
    data.x(:,t)=x;
    u=ur-K*(x-xr);
    data.u(t)=u;
    y=next_state(ts,@CL_VCCM,y,u);
end

end

% add the VRG dynamics
function dy=CL_VCCM(t,y,u)
x=y(1:3); y2=y(4);
dy=[dynamics(t,x,u); -y2+x(1)^2-2*x(1)*x(3)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CCM controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = run_ccm(x0)
% ccm=ccm_synthesis();
load('ccm.mat');
T=8; ts=0.002;
tm=0:ts:T; L=length(tm); N=7; 
data.tm=tm; 
data.x=zeros(3,L); data.u=zeros(1,L); 
data.E=zeros(1,L); data.x_val=cell(1,L);
data.xr=zeros(3,L); data.ur=zeros(1,L);

xr=[0;0;0]; ur=0; x=x0;

%% ccm control 
pts=chebpts(N); x_val=zeros(N,3); 
for j=1:N
    x_val(j,:)=(xr+(pts(j)+1)*(x-xr)/2)';
end

for t=1:L
    data.x(:,t)=x;
    data.x_val{1,t}=x_val;
    [u,data.E(t),x_val]=ccm_ctrl(x,xr,ur,ccm,x_val);
    data.u(t)=u;
    x=next_state(ts,@dynamics,x,u);
end
end

%
% CCM control realization via path integral
% geodesic computation: require chebfun
%

function [u,E,g_val]=ccm_ctrl(x,xr,ur,ccm,x_val)

% x_val -- warm start for geodesic computation
% \delta_u=-1/2\rho(x)B'W(x)^{-1}\delta_x
% u=ur+\int_{-1}^{1}\delta_u ds
% E is the Riemannian energy functional of a geodesic

N=size(x_val,1);
[g_val,E]=geodesic(xr',x',x_val,ccm);
gf=chebfun(g_val);
dgf=diff(gf);
pts=chebpts(N);
dg_val=dgf(pts);
du_val=0*pts;
Bt=[0,0,1];
for j=1:N
    x1=g_val(j,1);
    rho=[1,x1,x1^2]*ccm.rc;
    W=ccm.W0+x1*ccm.W1+x1^2*ccm.W2;
    du_val(j)=-0.5*rho*Bt/W*dg_val(j,:)';
end
duf=chebfun(du_val);
u=ur+sum(duf);
end

%-------------------------------------------------------------------------%
% find a geodesic \gamma with \gamma(-1)=x1 and \gamma(1)=x2. 
% ccm is the Riemannian metric
% A smooth path is parameterized by x_val at Chebyshev nodes.
% E_opt is the optimal Riemannian energy functional.
%-------------------------------------------------------------------------%
function [x_val_opt,E_opt]=geodesic(x1,x2,x_val,ccm)
options = optimset('Display','off',...
                'TolFun', 1e-6,...
                'MaxIter', 2000,...
                'Algorithm', 'interior-point',...
                'AlwaysHonorConstraints', 'bounds',...
                'FinDiffType', 'forward',...
                'HessFcn', [],...
                'Hessian', 'bfgs',...
                'HessMult', [],...
                'InitBarrierParam', 0.1,...
                'InitTrustRegionRadius', sqrt(size(x_val,1)*size(x_val,2)),...
                'MaxProjCGIter', 2*size(x_val,1)*size(x_val,2),...
                'ObjectiveLimit', -1e20,...
                'ScaleProblem', 'obj-and-constr',...
                'SubproblemAlgorithm', 'cg',...
                'TolProjCG', 1e-2,...
                'TolProjCGAbs', 1e-8);
A=[]; b=[]; Aeq = []; beq = []; lb=[]; ub=[];
[x_val_opt,E_opt]=fmincon(@(x_val)rieman_energy(x_val,ccm),x_val,A,b,Aeq,beq,lb,ub,...
               @(x_val)cons(x_val,x1,x2),options);
end

%-------------------------------------------------------------------------%
% compute Riemannian energy functional of a path c with Chebyshev 
% representation of x_val.
%-------------------------------------------------------------------------%
function E=rieman_energy(x_val,ccm)
[N,M]=size(x_val); 
% N is the number of Chebyshev nodes.
% M is the dimension of a point x.
pts=chebpts(N); % Chebyshev nodes: pts(j)=-cos(j\pi/N-1),j=0,1,...,N-1
c=chebfun(x_val);
dc=diff(c);
dc_val=dc(pts); % \delta_x
v_val=0*pts; % \delta_x'M(x)\delta_x
for j=1:N
    % metric related code
    x1=x_val(j,1);
    W=ccm.W0+x1*ccm.W1+x1^2*ccm.W2;
    v_val(j)=dc_val(j,:)/W*dc_val(j,:)';
end
vf=chebfun(v_val);
E=sum(vf); % integral
end

function [c,ceq]=cons(x_val,x1,x2)
c=[]; ceq=[x_val(1,:)-x1, x_val(end,:)-x2]';
end

% CCM synthesis via sum-of-squares (SOS)
function ccm=ccm_synthesis()
sdpvar x1 x2 x3 
x=[x1; x2; x3];
f=[-x1+x3; x1^2-x2-2*x1*x3+x3; -x2];
A=jacobian(f,x); B=[0; 0; 1];
A0=[-1,0,1; 0,-1,1; 0,-1,0]; 
Q=eye(3); R=1;
[K,P,e]=lqr(A0,B,Q,R);
W0=inv(P); r0=2/R;
wc1=sdpvar(6,1);
W1=[wc1(1),wc1(2),wc1(3); 
    wc1(2),wc1(4),wc1(5);
    wc1(3),wc1(5),wc1(6)];
wc2=sdpvar(6,1);
W2=[wc2(1),wc2(2),wc2(3); 
    wc2(2),wc2(4),wc2(5);
    wc2(3),wc2(5),wc2(6)];
W=W0+W1*x1+W2*x1^2;
rc=sdpvar(2,1);
rho=r0+rc(1)*x1+rc(2)*x1^2;
lambda=0.4;
v=sdpvar(3,1);
[c1,cc1]=polynomial([x;v],4);
[c2,cc2]=polynomial([x;v],4);
[c3,cc3]=polynomial([x;v],4);
[c4,cc4]=polynomial([x;v],4);
s1=sos(c1);
s2=sos(c2);
s3=sos(c3);
s4=sos(c4);
% constraint xlb <= x1 <= xub
xub=5; xlb=-5;
alpha1=0.01;
s5=sos(v'*(W-alpha1*eye(3))*v-(xub-x1)*c1-(x1-xlb)*c2);
s6=sos(rho);
dW=(W1+2*W2*x1)*f(1);
L=dW-W*A'-A*W+rho*B*B'-2*lambda*W;
s7=sos(v'*L*v-(xub-x1)*c3-(x1-xlb)*c4);
F=[s1,s2,s3,s4,s5,s6,s7];
ops = sdpsettings('solver','mosek','verbose',1,'debug',1);
coef=[wc1; wc2; rc; cc1; cc2; cc3; cc4];
optimize(F,[],ops,coef);
ccm.lambda=lambda;
ccm.W0=W0;%value(W0);
ccm.W1=value(W1);
ccm.W2=value(W2);
ccm.rc=[r0;value(rc)];
eps=1e-4;
for i=1:3
    for j=1:3
        if abs(ccm.W0(i,j)) < eps
            ccm.W0(i,j)=0;
        end
        if abs(ccm.W1(i,j)) < eps
            ccm.W1(i,j)=0;
        end
        if abs(ccm.W2(i,j)) < eps
            ccm.W2(i,j)=0;
        end
    end
    if abs(ccm.rc(i)) < eps
        ccm.rc(i)=0;
    end
end
end
