%% RK45 method
function x1=next_state(ts,x0,u0)
k1=ts*pendubot(x0,u0);
k2=ts*pendubot(x0+k1/2,u0);
k3=ts*pendubot(x0+k2/2,u0);
k4=ts*pendubot(x0+k3,u0);
x1=x0+(k1+2*k2+2*k3+k4)/6;
end