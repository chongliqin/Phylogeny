close all;
clear all;
%--------------------------------------------------------------------------
%Set initial conditions.
%--------------------------------------------------------------------------

Npos = 200;
NSeqs = 1;
N_samples_per_leaf = 1000;
steps_per_branch=[400,400];
%No interactions
J=eye(Npos);
%Turning on interactions.
J(100,2)=1;
J(2,100)=1;
J(100,30)=1;
J(30,100)=1;
% J(10,4)=0;
% J(10,4)=0;
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
        E_ij=-1*J.*(proteins(pct,:)'*proteins(pct,:))
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
proteins=repmat(proteins,2,1);
for count = 1:size(proteins,1)/2
    foo(2*count-1:2*count,:) = proteins([count,count+size(proteins,1)/2],:);
end
proteins=foo;

%--------------------------------------------------------------------------
%2nd Step of phylogeny
%--------------------------------------------------------------------------
for pct = 1:size(proteins,1) 
    E_ij=-1*J.*(proteins(pct,:)'*proteins(pct,:));
    E_0=0.5*sum(E_ij(:));
    for t = 1:steps_per_branch(2)
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
%Covariance matrices of the position (C) and sequences(C2).
%--------------------------------------------------------------------------
C=cov(proteins);
% C=proteins'*proteins;
% C = C./2000;
C2=cov(proteins');

%%
[V1,D1] = eig(C);
[D1,order]=sort(diag(D1),'descend');
V1=V1(:,order);
D1=D1;

[V,D] = eig(C2);
[D,order] = sort(diag(D),'descend');
V=V(:,order);
%D = D(1:size(D1,1));

%--------------------------------------------------------------------------
%Eigenvalue Calculation
%--------------------------------------------------------------------------
% D=eig(C);
% D2=eig(C2);


figure(1)
imagesc(C);
colorbar

figure(2)
imagesc(C2);
colorbar

figure(3)
hist(D,size(D,1));
%%
figure(4)
spacing = (max(D1)-min(D1))/100;
[height,bin] = hist(D1, min(D1):spacing:max(D1));
bar(bin,height./(spacing*Npos)),hold on;
% hist(D1,size(D1,1)),hold on;
scale = ceil(D1)/10;
X=min(D1-scale):0.01:max(D1+scale);
Y=size(D1,1)/size(D,1);
v=(max(D1)-min(D1))/((1+sqrt(Y))^2-(1-sqrt(Y))^2);
v=1;
a=v*(1-sqrt(Y))^2;
b=v*(1+sqrt(Y))^2;
f = real(sqrt((b-X).*(X-a))./(2.0*pi*X*Y*v));
plot(X,f,'r');

no_of_eig=4;

% for i=1:no_of_eig
%     figure(i+2)
%     imagesc(V(:,i)*V(:,i)');
%     colorbar
% end
% 
% for i=no_of_eig+1:2*no_of_eig
%     figure(i+2)
%     imagesc(V1(:,i-no_of_eig)*V1(:,i-no_of_eig)');
%     colorbar
% end
%figure(3)
% hist(D,100);
% 
% figure(4)
% hist(D2,100);

















