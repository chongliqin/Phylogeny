close all
clear all
%--------------------------------------------------------------------------
%The Initial Data
%--------------------------------------------------------------------------
NumOfSequences =10;
NumOfNodes =160;
NumOfMutations =400;
J=eye(NumOfNodes);
% J=randi(0:1,NumOfNodes);
% J=J+J-ones(size(J));
E=zeros(NumOfSequences,1);
E_new=zeros(NumOfSequences,1);
proteins = zeros(NumOfSequences,NumOfNodes);

for i=1:NumOfSequences
    proteins(i,:) = randi(0:1,1,NumOfNodes);
    proteins(i,:) = proteins(i,:)+proteins(i,:)-ones(size(proteins(i,:)));
    E(i)= proteins(i,:)*(J*proteins(i,:)');
    for j=1:NumOfMutations
        change = randi(NumOfNodes);
        proteins(i,change)=-proteins(i,change);
        E_new(i) = proteins(i,:)*(J*proteins(i,:)');
        delta_E = E_new(i)-E(i);
        if(delta_E>0)
            ProbAccept=exp(-delta_E);
            decider=rand;
            if ProbAccept<decider
                E_new(i) = E(i);
                proteins(i,change)=-proteins(i,change);
            end
        end
        E(i)=E_new(i);
    end   
end
C = cov(proteins);
C2 = cov(proteins');
figure(1)
imagesc(proteins);
figure(2)
imagesc(C);
figure(3)
imagesc(C2);