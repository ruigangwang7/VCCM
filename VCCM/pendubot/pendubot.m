%% pendubot nonlinear dynamics
function dx=pendubot(x,u)
% parameters
l1=0.1395; lg1=0.06975; l2=0.078; 
m1=0.115; mh=0.13; m2=0.0731;
fv1=1.3E-3; fv2=2.2E-5; g=9.8;
% dynamics
theta1=x(1); theta2=x(2); dtheta1=x(3); dtheta2=x(4); 
M=[(m1+m2)*l1^2, m2*l1*l2*cos(theta1-theta2);
    m2*l1*l2*cos(theta1-theta2), m2*l2^2];
C=[fv1,m2*l1*l2*sin(theta1-theta2)*dtheta2;
   -m2*l1*l2*sin(theta1-theta2)*dtheta1,fv2];
G=[-(m1*lg1+(mh+m2)*l1)*g*sin(theta1); -m2*l2*g*sin(theta2)];
T=[u; 0];
dx=0*x;
dx(1:2)=x(3:4);
dx(3:4)=M\(T-C*x(3:4)-G);
end