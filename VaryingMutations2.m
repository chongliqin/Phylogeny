%Keeping branching events and steps per branch constant - varying the
%number of positions.
clear all
close all
clear MutationVariance;
clear foo;
resolution=1;
% steps_per_branch=10;
branches=3;
N_seq = 1;
N_pos = 100;
alpha =zeros(700,1);
max_run=100;
for i=1:size(alpha,1)
    steps_per_branch=resolution*i;
   % alpha = zeros(size(alpha));
    for runs=1:max_run
        clear foo;
        proteins = randi(0:1,1,N_pos);
        proteins = 2*proteins-ones(size(proteins));
        % evolution and duplication....
        for bct = 1:branches
            for pct = 1:size(proteins,1)                
                % evolution step
                for t = 1:steps_per_branch
                    hot = randi(N_pos);
                    proteins(pct,hot) = -1*proteins(pct,hot);
                end              
            end            
            proteins = repmat(proteins,2,1);
            % duplication step
            for count = 1:size(proteins,1)/2
                foo(2*count-1:2*count,:) = proteins([count,count+size(proteins,1)/2],:);
            end
            proteins=foo;
        end       
        % N_seq=size(proteins,1);
        C=cov(proteins');
        C2=cov(proteins);
        alpha(i) = alpha(i)+C(1,5);
        beta(i,runs)=C(1,5); %+(i-1)*max_run
        %C = C-diag(diag(C));
%         [V,D] = eig(C);
    end
        alpha(i)=alpha(i)/max_run;
end

MSE = zeros(size(alpha));

for i=1:size(beta,1)
    for j=1:size(beta,2)
        MSE(i) = (beta(i,j)-alpha(i))*(beta(i,j)-alpha(i))+MSE(i);
    end
end
%%
figure(1)
plot((1:i)*resolution,alpha(:),'k');

figure(2)
plot((1:i)*resolution,beta(:,1:max_run),'k');

figure(3)
plot((1:i)*resolution, MSE(:),'k'),hold on;
X = (1:i)*resolution;
X =X(:);
Y = MSE;
myfittype=fittype('a*tanh(b*X)','dependent',{'Y'},'independent',{'X'},'coefficients',{'a','b'});
Trend=fit(X,Y,myfittype);
figure(4)
plot(Trend,X,Y);
Trend
%%
% figure(2)
% imagesc(C)
% colorbar
% 
% %C2=cov(proteins');
% figure(3)
% imagesc(C2)
% colorbar


