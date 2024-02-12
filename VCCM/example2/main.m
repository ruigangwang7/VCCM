%% Example 2
%  dynamical system: dx1/dt=2x1+x2+u, dx2/dt=x1+x1^2-x1^3-x2
%
function main
%% virtual linear system and LQR control design
% A=[2,1; 1,-1]; B=[1; 0]; [K,S,e]=lqr(A,B,eye(2),1);
K=[5, 1.5]; T=6;
% set-point
x1r=0; x2r=x1r+x1r^2-x1r^3; ur=-(2*x1r+x2r); xr=[x1r; x2r];
% simulate linear controller with VRG
x0=[7; 7; 0];
[t1,x1]=ode45(@(t,x)CL_VRG(t,x,xr,ur,K),[0,T],x0);
% simulate linear controller
x0=[7; 7]; 
[t2,x2]=ode45(@(t,x)CL_LINEAR(t,x,xr,ur,K),[0,T],x0);
save('example2.mat','t1','x1','t2','x2');

% load('example2.mat');
figure(1); hold on
plot(t1,x1(:,1),'r'); 
plot([0,6],[0,0],'r--'); 
plot(t1,x1(:,2),'b');
plot(t1,x1(:,3),'b--');
plot(t2,x2(:,1),'r-.'); 
plot(t2,x2(:,2),'b-.'); 
xlim([0,6]);
ylim([-40,20]);
legend('x_1 VRG-D','y_1^* VRG-D','x_2 VRG-D','y_2^* VRG-D','x_1 Linear','x_2 Linear');
xlabel('t'); 
set(gca,'fontsize',10);
grid on

savefig('example2.fig')
end

% Linear controller
function dx=CL_LINEAR(t,x,xr,ur,K)
x1=x(1); x2=x(2);
u=ur-K*(x-xr);
dx=[2*x1+x2+u; x1+x1^2-x1^3-x2];
end

% Linear controller with VRG
function dx=CL_VRG(t,x,xr,ur,K)
x1=x(1); x2=x(2); x2s=x(3); x1r=xr(1); x2r=xr(2);
mu=ur+x2r-x2s; xs=[x1r; x2s];
u=mu-K*(x(1:2)-xs);
dx=[2*x1+x2+u; x1+x1^2-x1^3-x2; x1r+x1^2-x1^3-x2s];
end

