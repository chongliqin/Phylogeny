close all;
clear all;
%--------------------------------------------------------------------------
%Set initial conditions.
%--------------------------------------------------------------------------

for test = 1:1
Npos = 160;
NSeqs = 1;
N_samples_per_leaf = 1000;
steps_per_branch=[400,10];
%No interactions
J=eye(Npos);
%Turning on interactions.
J(100,2)=0;
J(2,100)=0;
J(100,3)=0;
J(3,100)=0;
J(102,4)=0;
J(4,102)=0;
proteins = randi(0:1,1,Npos);
proteins = (proteins*2) - ones(size(proteins));
E_ij = -1*J.*(proteins'*proteins);
E_0 = 0.5*sum(E_ij(:));

%--------------------------------------------------------------------------
%1st step duplication.
%--------------------------------------------------------------------------

proteins =repmat(proteins,NSeqs*N_samples_per_leaf,1);

%--------------------------------------------------------------------------
%The star phylogeny evolution Step
%--------------------------------------------------------------------------

for pct = 1:size(proteins,1) 
    E_ij=-1*J.*(proteins(pct,:)'*proteins(pct,:));
    E_0=0.5*sum(E_ij(:));
    for t = 1:steps_per_branch(1)
        hot = randi(Npos);
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
        E_0=E_new;
    end
end

%--------------------------------------------------------------------------
%2nd Duplication and reordering step
%--------------------------------------------------------------------------
% 
% proteins=repmat(proteins,2,1);
% for count = 1:size(proteins,1)/2
%     foo(2*count-1:2*count,:) = proteins([count,count+size(proteins,1)/2],:);
% end
% proteins=foo;
% 
% %--------------------------------------------------------------------------
% %2nd Step of phylogeny
% %--------------------------------------------------------------------------
% 
% for pct = 1:size(proteins,1) 
%     E_ij=-1*J.*(proteins(pct,:)'*proteins(pct,:));
%     E_0=0.5*sum(E_ij(:));
%     for t = 1:steps_per_branch(2)
%         hot = randi(Npos);
%         proteins(pct,hot) = - proteins(pct,hot);
%         E_ij=-1*J.*(proteins(pct,:)'*proteins(pct,:));
%         E_new=0.5*sum(E_ij(:));
%         Delta_E=E_new-E_0;
%         if E_new>= E_0
%             prob_accept=exp(-Delta_E);
%             decider=rand;
%             if prob_accept<decider
%                 E_new=E_0;
%                 proteins(pct,hot)=-proteins(pct,hot);
%             end
%         end
%         E_0=E_new;
%     end
% end

%--------------------------------------------------------------------------
%Covariance matrices of the position (C) and seque'nces(C2).
%--------------------------------------------------------------------------

C=proteins*proteins';
C2=proteins'*proteins;

%--------------------------------------------------------------------------
%Eigenvalue Calculation
%--------------------------------------------------------------------------

D=eig(C);
D2=eig(C2);

D_store(:,test) = D(D>1);
D2_store(:,test) = D2(D2>1);

end
%%
figure(1)
imagesc(C);

figure(2)
imagesc(C2);

figure(3)
hist(D_store(:),100),hold on;
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')

figure(4)
hist(D2_store(:),100), hold on;
% h = findobj(gca,'Type','patch');
% set(h,'FaceColor',[0 .5 .5],'EdgeColor','w')

















