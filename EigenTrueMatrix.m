%This plots the exact eigenvalues of the true covariance matrix
clear all;
close all;
m=400;
n=160;
b=7;
a=exp(-4*m/n);
% 
% for i=1:2^(b-1)
%     EigVal(i)=1-a;
% end
% for i=2^(b-1)+1:2^(b-1)+2^(b-2)
%     EigVal(i)=1+a-2*a^2;
% end

s=0;
EigVal=zeros(2^b,1);
double(EigVal);
for j=1:b-1
    for i=s+1:s+2^(b-j)
        if j==1
            EigVal(i) =1-a;
        else
            EigVal(i) =1-2^(j-1)*a^(j);
            for k=1:j-1
            EigVal(i)=a^k*2^(k-1)+EigVal(i);   
            end
        end
    end
    s=i;
end
EigVal(2^b)=1;

for k=1:b
    EigVal(2^b)= EigVal(2^b)+a^k*2^(k-1);
end

EigVal(2^b-1)=EigVal(2^b)-2*a^b*2^(b-1);
figure(1)
hist(EigVal,2^b);

% 
% for i=1:b
%     for i=1:2^(b-1)
%         EigVal(i)=1-a;
%     end
%     for j=1:b
%     if  sum(2^s,s=0..j-1)<i<=sum(2^s,s=..j)
%         EigVal(i)=1+sum(2^(k-1)*a^k,k=1..j-1);
%     end
%     end
% end