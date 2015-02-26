%Keeping branching events and steps per branch constant - varying the
%number of positions.
clear MutationVariance;
steps_per_branch=20;
branches=10;
resolution=10;
%steps_per_branch=1;
N_seq = 1;
for s=1:2 % s changes the number of position. 
    N_pos = resolution*s;

    proteins = ones(N_seq,N_pos);
    
    clear foo;
    
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
  %  C2=cov(proteins');
    C = C-diag(diag(C));
    [V,D] = eig(C);
    C2 = C(1:size(C,1)/2,size(C,2)/2+1:size(C,2));
    C3 = C(1:size(C,1)/4,size(C,2)/4+1:size(C,2)/2);
    C4 = C(1:size(C,1)/8,size(C,2)/8+1:size(C,2)/4);
    C5 = C(1:size(C,1)/16,size(C,2)/16+1:size(C,2)/8);
    C6 = C(1:size(C,1)/32,size(C,2)/32+1:size(C,2)/16);
    C7 = C(1:size(C,1)/64,size(C,2)/64+1:size(C,2)/32);
    C8 = C(1:size(C,1)/128,size(C,2)/128+1:size(C,2)/64);
    C9 = C(1:size(C,1)/256,size(C,2)/256+1:size(C,2)/128);
    C10 = C(1:size(C,1)/512,size(C,2)/512+1:size(C,2)/1024);
    N_pos
    var(C(:))
    MutationVariance(s)=var(C(:));
    MutationVariance2(s)=var(C2(:));
    MutationVariance3(s)=var(C3(:));
    MutationVariance4(s)=var(C4(:));
    MutationVariance5(s)=var(C5(:));
    MutationVariance6(s)=var(C6(:));
    MutationVariance7(s)=var(C7(:));
    MutationVariance8(s)=var(C8(:));
    MutationVariance9(s)=var(C9(:));
   % MutationVariance10(s)=var(C10(:));
   
    MutationMean(s)=mean(C(:));
    MutationMean2(s)=mean(C2(:));
    MutationMean3(s)=mean(C3(:));
    MutationMean4(s)=mean(C4(:));
    MutationMean5(s)=mean(C5(:));
    MutationMean6(s)=mean(C6(:));
    MutationMean7(s)=mean(C7(:));
    MutationMean8(s)=mean(C8(:));
    MutationMean9(s)=mean(C9(:));  
    
    MatrixMean(s,:)=[mean(C2(:)),mean(C3(:)),mean(C4(:)),mean(C5(:)),mean(C6(:)),mean(C7(:)),mean(C8(:)),mean(C9(:))];
end
figure(1)
plot((10:s)*resolution, MutationVariance(10:s));%, hold on;
% plot((1:s)*resolution,MutationVariance2(:)), hold on;
% plot((1:s)*resolution,MutationVariance3(:)), hold on;
% plot((1:s)*resolution,MutationVariance4(:)), hold on;
% plot((1:s)*resolution,MutationVariance5(:)), hold on;
% plot((1:s)*resolution,MutationVariance6(:)), hold on;
% plot((1:s)*resolution,MutationVariance7(:)), hold on;
% plot((1:s)*resolution,MutationVariance8(:)), hold on;
% plot((1:s)*resolution,MutationVariance9(:)), hold on;
% plot((1:s)*resolution,MutationVariance10(:)), hold on;
xlabel('Number of positions')
ylabel('variance')
%%
figure(2)
imagesc(C)
colorbar

%C2=cov(proteins');
figure(3)
imagesc(C2)
colorbar


%%
figure(4)
imagesc(C3)
colorbar

%C2=cov(proteins');
figure(5)
imagesc(C4)
colorbar

%%
figure(6)
imagesc(C5)
colorbar

%C2=cov(proteins');
figure(7)
imagesc(C6)
colorbar

%%
figure(8)
imagesc(C7)
colorbar

%C2=cov(proteins');
figure(9)
imagesc(C8)
colorbar

%%
figure(10)
imagesc(C9)
colorbar

figure(11)
hist(C(:),100);
xlabel('Covariance Values');
ylabel('Frequency');

figure(12)
hist(C2(:),100);
xlabel('Covariance Values');
ylabel('Frequency');
%C2=cov(proteins');
% figure(11)
% imagesc(C10)
% colorbar
