clear all
close all
tic
colors = spring(10);
for ind = 1:10;
    n = ind;
    for runs=1:5000
        J = normrnd(0,1/sqrt(n),n);
        J = 0.5*(J+J');
        E = trace(J);
        D = 0:2^n-1;
        B = dec2bin(D);
        B=2*B-97;
        for i=1:2^n
            sig = B(i,:);
            E(i) = sig*J*sig';
        end
        E_store(runs,ind) = min(E);
    end
    E_average(ind) = mean(E_store(:,ind));
    E_variance(ind) = var(E_store(:,ind));
    figure(1);
    plot((1:runs),E_store(:,ind),'o','MarkerEdgeColor', colors(ind,:),'MarkerFaceColor',colors(ind,:)),hold all
end
    t=1:ind;
    figure(2);
    plot((1:ind),E_average(:),'k-');
    figure(3);
    plot((1:ind),E_variance(:),'k-');
toc