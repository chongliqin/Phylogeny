clear all
close all
tic
colors = spring(30);
for ind = 1:30;
    n = ind;
    for runs=1:500
        J = normrnd(0,1/sqrt(n),n);
        J = 0.5*(J+J');
        E = trace(J);
        for row=1:n
            for column=row+1:n
                if E+2*J(row,column)<E-2*J(row,column)
                    E=E+2*J(row,column);
                else
                    E=E-2*J(row,column);
                end
            end
        end
        E_store(runs,ind) = E;
    end
    E_average(ind) = mean(E_store(:,ind));
    E_variance(ind) = var(E_store(:,ind));
    figure(1);
    plot((1:runs),E_store(:,ind),'o','MarkerEdgeColor', colors(ind,:),'MarkerFaceColor',colors(ind,:)),hold all
end
    t=1:ind;
    figure(2);
    plot((1:ind),E_average(:),'k-');
    figure(3);
    plot((1:ind),E_variance(:),'k-');
toc