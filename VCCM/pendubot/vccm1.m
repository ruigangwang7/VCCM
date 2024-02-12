%% virtual system (38a)
function [tm, xx] = vccm1()

% Kp=synthesis();
load('VCCM_Kp_1.mat');

l1=0.1395; lg1=0.06975; m1=0.115; mh=0.13; m2=0.0731; g=9.8;

% simulation setup
ts=0.01; T=20; tm=0:ts:T; L=length(tm); xx=zeros(4,L); uu=zeros(1,L);
x=zeros(4,1); r=[0.5,1/4*pi,0,-1/4*pi];
for kk=1:L
    xx(:,kk)=x;
    ii=floor(tm(kk)/5);
    if ii > 3
        ii=3;
    end
    xs=[r(ii+1);0;0;0]; us=-(m1*lg1+(mh+m2)*l1)*g*sin(xs(1));
    u=ctrl(x,xs,us,Kp);
    uu(kk)=u;
    % disp([tm(kk),x',u]);
    x=next_state(ts,x,u);
    if abs(u) > 5
        break
    end
end

end

function Kp=synthesis()

% parameters
l1=0.1395; lg1=0.06975; l2=0.078; 
m1=0.115; mh=0.13; m2=0.0731;
fv1=1.3E-3; fv2=2.2E-5; g=9.8;

% control design
M=@(x) [(m1+m2)*l1^2, m2*l1*l2*cos(x);
    m2*l1*l2*cos(x), m2*l2^2];
syms y1 y2 y3 y4
f=[fv1*y3+m2*l1*l2*sin(y1-y2)*y4^2; -m2*l1*l2*sin(y1-y2)*y3^2+fv2*y4]+[-(m1*lg1+(mh+m2)*l1)*g*sin(y1); -m2*l2*g*sin(y2)];
Df=matlabFunction(jacobian(f,[y1; y2; y3; y4]));

x1m=-1/3*pi; x1M=1/3*pi; x1N=11; 
x1g=linspace(x1m,x1M,x1N);
y1m=-1/3*pi; y1M=1/3*pi; y1N=11;
y1g=linspace(y1m,y1M,y1N);
y2m=-1/18*pi; y2M=1/18*pi; y2N=5;
y2g=linspace(y2m,y2M,y2N);
y3m=-1; y3M=1; y3N=5;
y3g=linspace(y3m,y3M,y3N);
y4m=-1; y4M=1; y4N=5;
y4g=linspace(y4m,y4M,y4N);

ny=4; nv=1;
W=sdpvar(ny);
sdpvar t
Yp=sdpvar(nv,ny,x1N,y1N,y2N,y3N,y4N,'full');
LMI = [W>=0.01*eye(ny); W <= t*eye(ny); t >= 0];
A=zeros(4,4); A(1,3)=1; A(2,4)=1; B=zeros(4,1);
lbd=0.5;
for x1i=1:x1N
    Mi=inv(M(x1g(x1i))); B(3:4)=Mi(:,1);
    for y1i=1:y1N
        for y2i=1:y2N
            for y3i=1:y3N
                for y4i=1:y4N
                    A(3:4,:)=-Mi*Df(y1g(y1i),y2g(y2i),y3g(y3i),y4g(y4i));
                    Y=Yp(:,:,x1i,y1i,y2i,y3i,y4i);
                    H=(A*W+B*Y)+(A*W+B*Y)'+2*lbd*W;
                    LMI=[LMI; H <= 0];
                end
            end
        end
    end
end
optimize(LMI,t)
W=double(W);
Yp=double(Yp);
% Yp(:,:,1,1,1,1,1)/W
Kp.dom=[x1m,x1M,x1N;
        y1m,y1M,y1N;
        y2m,y2M,y2N;
        y3m,y3M,y3N;
        y4m,y3M,y4N];
Kp.K=zeros(nv,ny,x1N,y1N,y2N,y3N,y4N);
for x1i=1:x1N
    for y1i=1:y1N
        for y2i=1:y2N
            for y3i=1:y3N
                for y4i=1:y4N
                    Y=Yp(:,:,x1i,y1i,y2i,y3i,y4i);
                    Kp.K(:,:,x1i,y1i,y2i,y3i,y4i)=Y/W;
                end
            end
        end
    end
end
save('VCCM_Kp_1.mat','Kp');
end

%% realization
function u=ctrl(x,xs,us,Kp)
L=20; K=zeros(1,4);
x1=x(1)-x(2); 
for jj=0:L
    y=(xs*(L-jj)+x*jj)/L;
    Ks=lookup_table(Kp,[x1;y]);
    K=K+Ks;
end
K=K/(L+1);
u=us+K*(x-xs);
end