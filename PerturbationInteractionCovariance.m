%Keeping mutations per branch and number of positions constant - varying
%the number of branching events.
close all;
clear all;
clear MutationVariance;
N_pos=160;
steps_per_branch=200;
resolution=1;
branches = 10;
proteins = randi(0:1,1,N_pos);
proteins = (proteins*2)-ones(size(proteins));
J=eye(N_pos);
J(100,2) =0;
J(2,100) =0;
E=proteins(:)'*(J*proteins(:));
E_new=E;
clear foo;

% evolution and duplication....
for bct = 1:branches
    for pct = 1:size(proteins,1)
        
        % evolution step
        for t = 1:steps_per_branch
            hot = randi(N_pos);
            proteins(pct,hot) = - proteins(pct,hot);
            E_new =proteins(pct,:)*(J*proteins(pct,:)');
            DeltaE=E_new-E;
            if(DeltaE>0)
                ProbAccept=exp(-DeltaE);
                decider=rand;
                if ProbAccept<rand
                    E_new=E;
                    proteins(pct,hot) = - proteins(pct,hot);
                end
            end
            E=E_new;
        end
        
    end
    
    proteins = repmat(proteins,2,1);
    % duplication step
    for count = 1:size(proteins,1)/2
        foo(2*count-1:2*count,:) = proteins([count,count+size(proteins,1)/2],:);
    end
    proteins=foo;
end

%Position Covariance
C=cov(proteins);

%Sequence Covariance
C2=cov(proteins');

[V,D]=eig(C2);
[V1,D1]=eig(C);
% [V,D]=eig(C);
% This sorts the eigenvalues and corresponding eigenvectors in descending
% order.
[D,order] = sort(diag(D), 'descend');
[D1,order1] = sort(diag(D1), 'descend');
V=V(:,order);
V1=V(:,order1);
%% 
figure(1)
imagesc(C)
colorbar

figure(2)
imagesc(C2)
colorbar

figure(3)
hist(D1,200);

for i=1:4
    figure(i+3)
    imagesc(V(:,i)*V(:,i)');
    colorbar;
end

% figure(3)
% imagesc(V(:,1)*V(:,1)');
% colorbar
% 
% figure(4)
% imagesc(V(:,2)*V(:,2)');
% colorbar
% 
% figure(5)
% imagesc(V(:,3)*V(:,3)');
% colorbar
% 
% figure(6)
% imagesc(V(:,4)*V(:,4)');
% colorbar
% var(C(:))
% MutationVariance(s)=var(C(:));
% 
% figure(1)
% plot((1:s)*resolution,MutationVariance(:));
% xlabel('Number of Branch Events')
% ylabel('Variance')