%Keeping mutations per branch and number of positions constant - varying
%the number of branching events.
close all;
clear all;
clear MutationVariance;
N_pos=160;
steps_per_branch=200;
resolution=1;
branches = 5;
proteins = randi(0:1,1,N_pos);
proteins = (proteins*2)-ones(size(proteins));
clear foo;

% evolution and duplication....
for bct = 1:branches
    for pct = 1:size(proteins,1)
        
        % evolution step
        for t = 1:steps_per_branch
            hot = randi(N_pos);
            proteins(pct,hot) = - proteins(pct,hot);
        end
        
    end
    
    proteins = repmat(proteins,2,1);
    % duplication step
    for count = 1:size(proteins,1)/2
        foo(2*count-1:2*count,:) = proteins([count,count+size(proteins,1)/2],:);
    end
    proteins=foo;
end
C=cov(proteins);
C2=cov(proteins');

[V,D]=eig(C2);
% This sorts the eigenvalues and corresponding eigenvectors in descending
% order.
[D,order] = sort(diag(D), 'descend');
V=V(:,order);
% 
% figure(1)
% imagesc(C)
% colorbar

figure(2)
imagesc(C2)
colorbar

for i=1:4
    figure(i+2)
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