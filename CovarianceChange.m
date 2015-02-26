%The change in covariance as a sequence changes from 1 to -1
clear all;
close all;
N_pos=100;
Initial_Sequence = randi(1,1,N_pos);
Secondary_Sequence = Initial_Sequence;
Comparison_Matrix =[Initial_Sequence; Secondary_Sequence]';
C = cov(Comparison_Matrix);
C2(1) = C(2,2);
for i=1:N_pos
    Secondary_Sequence(i) = 0;
    Comparison_Matrix =[Secondary_Sequence; Secondary_Sequence]';
    C = cov(Comparison_Matrix)
    C2(i+1) = C(1,2);
end;
figure(1);
plot((1:N_pos+1),C2(:));
    