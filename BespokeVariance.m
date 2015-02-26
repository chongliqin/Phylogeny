close all;
clear all;
% This generates S number of sequences with N number of positions
%-------------------------------------------------------------------
N=50;
S=1000;
m=200;
resolution=0.01;
s=100.0;
J=eye(N);
T=100;

for alpha=41:99
al=0.01*alpha;
TaylorCov=eye(N);
TaylorCov(1,1)=sqrt(1+al);
TaylorCov(2,2)=sqrt(1-al);

% TaylorCov=sqrt(3)*eye(N);

% J(1,2)=s*resolution;
% J(2,1)=J(1,2);
% This is the interaction matrix
%----------------------------------
for t=1:T
%for s=1:200
%   m=s;
%for s=1:2000
% J(3,4)=0.5;
% J(4,3)=J(3,4);
% J(1,3)=1.0;
% J(3,1)=J(1,3);
% J(2,3)=1;
% J(3,2)=J(2,3);
% J(1,3)=1;
% J(3,1)=J(1,3);
proteins = zeros(S,N);
%Generating the proteins to be a sequence of -1 and 1s.
%-------------------------------------------------------
for i=1:S
    proteins(i,:) = randi(0:1,1,N);
    proteins(i,:) = (proteins(i,:)*2) - ones(size(proteins(i,:)));
end
proteins=proteins*TaylorCov;
for i=1:S
    %Calculating the corresponding energy.
    %-------------------------------------
    E_o(i) = 0.5*proteins(i,:)*(J*proteins(i,:)');
end

%Monte Carlo Algorithm
%-----------------------
for i=1:S
    for M=1:m
        hot = randi(N);
        proteins(i,hot)=-proteins(i,hot);
        E_new(i) = 0.5*proteins(i,:)*(J*proteins(i,:)');
        DeltaE= E_new(i)-E_o(i);
        if(DeltaE>0)
            ProbAccept=exp(-DeltaE);
            decider=rand;
            if rand>ProbAccept
                E_new(i)=E_o(i);
                proteins(i,hot)=-proteins(i,hot);
            end
        end
        E_o(i) =E_new(i);
    end
end

C = cov(proteins);
[V,D] = eig(C);
EV(t,:) = diag(D);
end
D1=EV(:);
D1=sort(D1);

FarBulk=D1(size(D1,1)-T+1:size(D1,1));
[FBheight,FBbin]=hist(FarBulk,100);
spacing=(max(FarBulk)-min(FarBulk))/100;
FBheight=double(FBheight(:))./(spacing*size(FarBulk,1));
FBbin=FBbin(:);
fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1]);
GaussianFit = fittype('1.0/(sqrt(2.0*pi)*Sig)*exp(-(FBbin-mean)^2/(2.0*Sig^2))', 'dependent', {'FBheight'}, 'independent', {'FBbin'}, 'coefficients',{'Sig','mean'}, 'options', fo);
GFit=fit(FBbin,FBheight,GaussianFit);
FBparam(alpha,:) =coeffvalues(GFit);


LeftBulk=D1(1:T);
[LBheight,LBbin]=hist(LeftBulk,100);
spacing=(max(LeftBulk)-min(LeftBulk))/100;
LBheight=LBheight(:)./(spacing*size(LeftBulk,1));
LBbin=LBbin(:);
GaussianFit = fittype('1.0/(sqrt(2.0*pi)*Sig)*exp(-(FBbin-mean)^2/(2.0*Sig^2))', 'dependent', {'FBheight'}, 'independent', {'FBbin'}, 'coefficients',{'Sig','mean'}, 'options', fo);
Git=fit(LBbin,LBheight,GaussianFit);
LBparam(alpha,:)=coeffvalues(Git);

end
% Data(s) =C(1,2);
% Diag(s) = C(1,1);
% Data2(s) = C(1,3);
% Data3(s) = C(2,3);
% end;
%%
close all
% X = (1:s);
% X =X(:);
% Y = Data(:);
% figure(2)
% scatter((1:s), Data(:));
% myfittype=fittype('b*tanh(a*X)','dependent',{'Y'},'independent',{'X'},'coefficients',{'a','b'});
% Trend=fit(X,Y,myfittype);
% figure(4)
% plot(Trend,X,Y);
% Trend
% Z=Data3(:);
% myfittype=fittype('b*tanh(a*X)','dependent',{'Z'},'independent',{'X'},'coefficients',{'a','b'});
% Trend=fit(X,Z,myfittype);
% figure(3)
% plot(Trend,X,Z);
% figure(3)
% scatter((1:s)*resolution,Data2(:));
% figure(5)
% scatter((1:s)*resolution,Data3(:));
% %Figures of outputs.
%--------------------
%%
% % % % % % % % close all;
% % % % % % % % figure(1)
% % % % % % % % imagesc(C);
% % % % % % % % colorbar;
% % % % % % % % figure(2)
% % % % % % % % D1=EV(:);
% % % % % % % % hist(EV(:),1000);
% % % % % % % % 
% % % % % % % % % figure(3);
% % % % % % % % % D2=D1(D1>2);
% % % % % % % % % hist(D2(:),100);
% % % % % % % % % 
% % % % % % % % % figure(4)
% % % % % % % % % D3=D1(D1<0.25);
% % % % % % % % % hist(D3(:),100);
% % % % % % % % D1=sort(D1);
% % % % % % % % 
% % % % % % % % spacing = (max(D1)-min(D1))/1000;
% % % % % % % % [height,bin] = hist(D1, min(D1):spacing:max(D1));
% % % % % % % % bar(bin,height./(spacing*size(D1,1))),hold on;
% % % % % % % % % hist(D1,size(D1,1)),hold on;
% % % % % % % % scale = ceil(D1)/10;
% % % % % % % % X=min(D1-scale):0.01:max(D1+scale);
% % % % % % % % Y=size(proteins,2)/size(proteins,1);
% % % % % % % % v=(max(D1)-min(D1))/((1+sqrt(Y))^2-(1-sqrt(Y))^2);
% % % % % % % % v=1;
% % % % % % % % a=v*(1-sqrt(Y))^2;
% % % % % % % % b=v*(1+sqrt(Y))^2;
% % % % % % % % f =real(sqrt((b-X).*(X-a))./(2.0*pi*X*Y*v));
% % % % % % % % plot(X,f,'r');
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % 
% % % % % % % % FarBulk=D1(size(D1,1)-T+1:size(D1,1));
% % % % % % % % [FBheight,FBbin]=hist(FarBulk,100);
% % % % % % % % spacing=(max(FarBulk)-min(FarBulk))/100;
% % % % % % % % %sum=sum(FBheight,2);
% % % % % % % % FBheight=double(FBheight(:))./(spacing*size(FarBulk,1));
% % % % % % % % FBbin=FBbin(:);
% % % % % % % % %fo = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1,1]);
% % % % % % % % GaussianFit = fittype('1.0/(sqrt(2.0*pi)*S)*exp(-(FBbin-m)^2/(2.0*S^2))', 'dependent', {'FBheight'}, 'independent', {'FBbin'}, 'coefficients',{'S','m'});
% % % % % % % % GFit=fit(FBbin,FBheight,GaussianFit);
% % % % % % % % figure(3)
% % % % % % % % bar(FBbin,FBheight),hold on;
% % % % % % % % plot(GFit);
% % % % % % % % % 
% % % % % % % % LeftBulk=D1(1:T);
% % % % % % % % [LBheight,LBbin]=hist(LeftBulk,100);
% % % % % % % % spacing=(max(LeftBulk)-min(LeftBulk))/100;
% % % % % % % % LBheight=LBheight(:)./(spacing*size(LeftBulk,1));
% % % % % % % % LBbin=LBbin(:);
% % % % % % % % 
% % % % % % % % Git=fit(LBbin,LBheight,GaussianFit);
% % % % % % % % figure(4)
% % % % % % % % bar(LBbin,LBheight),hold on;
% % % % % % % % plot(Git);


% fo2=fitoptions('Method','NonlinearLeastSquares','StartPoint', [1,1,0.2]);
% GaussianFit = fittype('a*exp(-b*(FBbin-c)^2)', 'dependent', {'FBheight'}, 'independent', {'FBbin'}, 'coefficients',{'a','b','c'},'options', fo2);
% axis([0.14 0.24 0 12])
% figure(2)
% plot((1:s)*resolution, Data(:), 'k'),hold on;
% J
% Y=-inv(J)

% X =(1:s)*resolution;
% Y =1./(2.*X);
% plot(X, Y, 'k');
