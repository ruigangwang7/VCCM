%% Pendubot Example: require YALMIP, SEDUMI, SYMBOLIC
%  control of pendubot via different choices of virtual systems
%
function main

lw=0.5;
c=[0,0,0.9];
p=[0,0.5,0];
o=[0.9,0,0];

r=[0.5,1/4*pi,0,-1/4*pi];
T=20;

% load('vccm-1.mat'); 
[tm, xx]=vccm1();
save('vccm-1.mat','tm','xx');
subplot(121);
hold on
plot(tm,xx(1,:),'color',c,'linewidth',lw); 
subplot(122);
plot(tm,xx(2,:),'color',c,'linewidth',lw); hold on

% load('vccm-2.mat');
[tm, xx]=vccm2();
save('vccm-2.mat','tm','xx');
subplot(121);
plot(tm,xx(1,:),'color',p,'linewidth',lw);
subplot(122);
plot(tm,xx(2,:),'color',p,'linewidth',lw);

% load('vccm-3.mat');
[tm, xx]=vccm3();
save('vccm-3.mat','tm','xx');
subplot(121);
plot(tm,xx(1,:),'color',o,'linewidth',lw);
plot([0 5 5 10 10 15 15 20], [r(1),r(1),r(2),r(2),r(3),r(3),r(4),r(4)],'k--');
xlim([0,T]);
ylim([-pi/3,pi/3]);
ylabel('$\theta_1$','FontSize',14,'interpreter','latex');
grid on
xlabel('t');
set(gca,'fontsize',12);

subplot(122);
plot(tm,xx(2,:),'color',o,'linewidth',lw);
xlim([0,T]);
ylim([-pi/18,pi/18]);
grid on
ylabel('$\theta_2$','FontSize',14,'interpreter','latex');
xlabel('t');
set(gca,'fontsize',12);
% legend({'virtsys-1','virtsys-2','virtsys-3'},'fontsize',10);

savefig('pendubot.fig')
end
