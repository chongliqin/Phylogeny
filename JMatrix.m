clear all 
close all
tic
colors = spring(30);
ind=1;
for S=1:1;
for ind = 1:1;
    n = ind;
    for runs=1:1
        J = normrnd(0,1,n);
        J = 0.5*(J+J');
        D = 0:2^n-1;
        B = dec2bin(D);
        B2=2*B-97;
        for i=1:2^n
         sig=randi([0:1],1,ind);
         sig=2*sig-1;
           sig= B2(i,:);
            E = sig*J*sig';
        end
        E_store(runs,ind) = E;
    end
    E_average(ind) = mean(E_store(:,ind));
    E_variance(ind) = var(E_store(:,ind));
   figure(1);
   plot((1:runs),E_store(:,ind),'o','MarkerEdgeColor', colors(ind,:),'MarkerFaceColor',colors(ind,:)),hold all
   figure; 
end
    Variance(S) = var(E_average);
end
   %  t = 1:S;
   %  figure(2);
   %  plot((1:S),Variance(t),'k-');
 %   figure(2);
  %  plot((1:ind),E_average(:)), hold on;
  %  plot((1:ind),sqrt(t)/sqrt(30), 'r-');
   % plot((1:ind),-sqrt(t)/sqrt(30), 'r-');
%%
%  t = 1:ind;
%     figure(3);
%     plot((1:ind),E_variance(:)), hold on;
%     plot((1:ind),2*t.^2-t, 'r-');
    %figure(4);
    
    toc
%end