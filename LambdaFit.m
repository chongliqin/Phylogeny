X= [50;100;200;300;400;500;800];
Y= [18.84/115.6;9.33/115.6;4.643/115.6;6.208/231.1;4.651/231.1;3.709/231.1;2.319/231.1];
Trend=fit(X,Y,'power1');
figure(2)
plot(Trend,X,Y);


X=[10;20;40;50;80;100];
Y=[0.2415;0.1053;0.05257;0.04371;0.0256;0.0212];
Trend2=fit(X,Y,'power1');
figure(3)
plot(Trend2,X,Y);

X=(41:90)*0.01;
Y=FBparam(2,41:90);
X=X(:);
Y=Y(:);
Trend3=fit(X,Y,'poly1');
figure(4)
plot(Trend3,X,Y);

X=(41:90)*0.01;
Y=LBparam(2,41:90);
X=X(:);
Y=Y(:);
Trend4=fit(X,Y,'poly1');
figure(5)
plot(Trend4,X,Y);

