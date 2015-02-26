clear all
close all
tic
colors = spring(30);
for S=1:100;
    n=10;
    for ind=1:100
    for runs=1:S
        J1 = zeros(n);
        J = triu(normrnd(0,1,n),0);
   %J = J1 + (triu(J1,1))';
       % J = 0.5*(J+J');
        sig=randi([0:1],1,n);
        sig=2*sig-1;        
        E = sig*J*sig';
        E_store(runs) = E;
    end
    E_average(ind,S) = mean(E_store(:));
    end
    E_variance(S) = var(E_average(:,S));
end
t =1:0.1:S;
figure(1);
plot((1:S), E_variance(:),'k-'), hold on;
plot((1:0.1:S), 0.5*n*(n+1)./t(:) ,'r-');
toc