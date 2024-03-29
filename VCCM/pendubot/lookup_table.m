%% interpolation
function K=lookup_table(Kp,rho)
N=length(rho); id= [ones(N,1),Kp.dom(:,3)]; w=0.5*ones(N,2); 
step=Kp.dom(:,2)-Kp.dom(:,1);
for k=1:N
    x=rho(k);
    if x<=Kp.dom(k,1)
        id(k,2)=id(k,1);
    else
        if x >= Kp.dom(k,2)
            id(k,1)=id(k,2);
        else
            id(k,1)=1+floor((x-Kp.dom(k,1))/step(k));
            id(k,2)=id(k,1)+1;
            w(k,1)=id(k,1)-(x-Kp.dom(k,1))/step(k);
            w(k,2)=1-w(k,1);
        end
    end
end
K=zeros(size(Kp.K,1),size(Kp.K,2));
for k=0:2^N-1
    c=1; b=ones(N,1);
    for j=1:N
        if bitget(k,j)
            b(j)=b(j)+1;
        end
        c=c*w(j,b(j));
    end
    K=K+c*Kp.K(:,:,id(1,b(1)),id(2,b(2)),id(3,b(3)),id(4,b(4)),id(5,b(5)));
end
end