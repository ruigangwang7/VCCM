%% Example 4:
%  dynamics: dx1/dt=-x1-x2+r, dx2/dt=1-exp(-x2)+u
%
function main
lw=1; 
ts=0.05; T=30; tm=0:ts:T; r=0*tm; L=length(tm); x=zeros(2,L);
for k=1:L
    if tm(k) <= 10 || tm(k) >= 20
            r(k)=3;
        else
            r(k)=-3;
    end
end
figure(1); 
subplot(121); hold on
plot(tm,r,'k--','linewidth',lw);
for i=1:3
    x(:,1)=[0;5];
    for k=1:L-1
        x(:,k+1)=rk45(ts,@(t,x)dyn(t,x,r(k),0,i),x(:,k));
    end
    switch i
        case 1
            c=[1, 0, 0];
        case 2
            c=[0,0.5,0];
        case 3
            c=[0, 0, 1];
    end
    plot(tm,x(2,:),'color',c,'linewidth',lw);
end
xlim([0,T]); ylim([-5,8]); grid on
xlabel('t'); ylabel('x_2');
legend('Ref','GSC 1', 'GSC 2', 'VCCM');

q=2;
ts=0.05; T=12; tm=0:ts:T; r=3*sin(q*tm); dr=q*3*cos(q*tm); L=length(tm); x=zeros(2,L);
subplot(122); hold on
plot(tm,r,'k--','linewidth',lw);
for i=1:3
    x(:,1)=[0;5];
    for k=1:L-1
        x(:,k+1)=rk45(ts,@(t,x)dyn(t,x,r(k),dr(k),i),x(:,k));
    end
    switch i
        case 1
            c=[1, 0, 0];
        case 2
            c=[0,0.5,0];
        case 3
            c=[0, 0, 1];
    end
    plot(tm,x(2,:),'color',c,'linewidth',lw);
end
xlim([0,T]); ylim([-5,8]); grid on
xlabel('t'); ylabel('x_2');
legend('Ref','GSC 1', 'GSC 2', 'VCCM');

savefig('example4.fig');
end

function dx=dyn(t,x,r,dr,i)
x1=x(1); x2=x(2); 
switch i
    case 1
        u=exp(-r)-1+x1-(3+exp(-r))*(x2-r);
    case 2
        u=x1+exp(-x2)-1;
    case 3
        u=dr+x1-3*(x2-r)+exp(-x2)-1;
end
dx=[-x1-x2+r; 1-exp(-x2)+u];
end

function x1=rk45(tau,model,x0)
k1=tau*model(0,x0);
k2=tau*model(tau/2,x0+k1/2);
k3=tau*model(tau/2,x0+k2/2);
k4=tau*model(tau,x0+k3);
x1=x0+(k1+2*k2+2*k3+k4)/6;
end