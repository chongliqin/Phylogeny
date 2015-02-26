%Keeping branching events and positions constant - varying the mutations
%per branch.
clear all;
close all;
clear MutationVariance;
N_pos=160;
branches=10;
%resolution=1;
%steps_per_branch=1;

J = eye(N_pos);
%interactions C is similar to J.
J(100,2)=0;
J(2,100)=0;
J(100,3)=0;
J(3,100)=0;
J(102,4)=0;
J(4,102)=0;
tic
%for s=1:400   
    C =zeros(2^branches);
    for j=1:1
    steps_per_branch = 10;
    proteins = randi(0:1,1,N_pos);
    proteins = (proteins*2) - ones(size(proteins));
    E_ij = -1*J.*(proteins'*proteins);
    E_0 = 0.5*sum(E_ij(:));
    clear foo;
    % evolution and duplication....
    for bct = 1:branches
        for pct = 1:size(proteins,1)
            
            % evolution step
            for t = 1:steps_per_branch
                hot = randi(N_pos);
                proteins(pct,hot) = - proteins(pct,hot);
                E_ij=-1*J.*(proteins(pct,:)'*proteins(pct,:));
                E_new=0.5*sum(E_ij(:));
                Delta_E=E_new-E_0;
                if E_new>= E_0
                    prob_accept=exp(-Delta_E);
                    decider=rand;
                    if prob_accept<decider
                        E_new=E_0;
                        proteins(pct,hot)=-proteins(pct,hot);
                    end
                end
            end
            E_0=E_new;
        end
        
        proteins = repmat(proteins,2,1);
        % duplication step
        for count = 1:size(proteins,1)/2
            foo(2*count-1:2*count,:) = proteins([count,count+size(proteins,1)/2],:);
        end
        proteins=foo;
    end
     C=cov(proteins');
     C2=cov(proteins);
    end
    
    figure(1)
    imagesc(C2);
    colorbar
    
    figure(2)
    imagesc(C);
    colorbar
   %%
   [V_C,D_C] =eig(C);
   D_C = eig(C);
   [V_C2,D_C2] = eig(C2);
   D_C2=eig(C2);
   
   figure(3)
   hist(C2(:),50);
%%   
   figure(4)
   hist(D_C2(:),50);
   
%%
n = N_pos;
[LAMBDA, lambda_idx] = sort(D_C2,'descend');
EVT = V_C2(:,lambda_idx);

greyColor = [120, 120, 120]/255;
evtA = EVT(:,1);
t1 = max(abs(evtA));
evtB = EVT(:,2);
t2 = max(abs(evtB));
t = 2/sqrt(n);
figure
set(gca,'fontsize',15,'fontname','Helvetica','fontweight','b')
otherHandle = plot(evtA, evtB, 'o', 'markerFaceColor', greyColor, 'markerEdgeColor', greyColor, 'markerSize', 4);
for i = 1:n
    text(evtA(i), evtB(i), num2str(i));
end
hold on
circle([0,0],t,1000,'k--'); 
plot([-t, t], [0, 0], 'k--')
plot([0, 0], [-t, t], 'k--')

% legend([h1, h2, h3, h4, otherHandle], 'sector 1', 'sector 2', 'sector 3', 'sector 4', 'other component pairs')
%title('1st Eigenvector against 2nd Eigenvector')

set(gca, 'xTick', [-t1, 0, t1]);
set(gca, 'yTick', [-t2, 0, t2]);
xlim([-t1, t1]);
ylim([-t2, t2]);
xlabel('EVT 1')
ylabel('EVT 2')

greyColor = [120, 120, 120]/255;
evtA = EVT(:,length(LAMBDA)-1);
t1 = max(abs(evtA));
evtB = EVT(:,length(LAMBDA));
t2 = max(abs(evtB));
figure
set(gca,'fontsize',15,'fontname','Helvetica','fontweight','b')
otherHandle = plot(evtA, evtB, 'o', 'markerFaceColor', greyColor, 'markerEdgeColor', greyColor, 'markerSize', 4);
for i = 1:n
    text(evtA(i), evtB(i), num2str(i));
end
hold on
circle([0,0],t,1000,'k--'); 
plot([-t2, t2], [0, 0], 'k--')
plot([0, 0], [-t2, t2], 'k--')

% legend([h1, h2, h3, h4, otherHandle], 'sector 1', 'sector 2', 'sector 3', 'sector 4', 'other component pairs')
%title('1st Eigenvector against 2nd Eigenvector')

set(gca, 'xTick', [-t1, 0, t1]);
set(gca, 'yTick', [-t2, 0, t2]);
xlim([-t1, t1]);
ylim([-t2, t2]);
xlabel(horzcat('EVT ',num2str(length(LAMBDA)-1)))

ylabel(horzcat('EVT ',num2str(length(LAMBDA))))

toc    
   % N_seq=size(proteins,1);
   
   % C2(s)=C(3,7)/100;
    %clear C;
  %  C2=cov(proteins');
    %C = C-diag(diag(C));
   % steps_per_branch
    %var(C(:))
    %MutationVariance(s)=var(C(:));
%end
%  t=[1:s]'.*resolution;
%  C2=C2';
%  BestFit=fit(t,C2,'exp1','Normalize','on');
%  figure(1)
%  plot(BestFit,t,C2);
%  xlabel('Number of Mutations per Branch')
%  ylabel('covariance')
