



## Chemotaxis
x=linspace(0, 1);
C0 = 1;
alpha = -5; 
kinvD = 0.1;

lambda1 = 0.5*(alpha + sqrt(alpha^2 - 4*kinvD));
lambda2 = 0.5*(alpha - sqrt(alpha^2 - 4*kinvD));

a = (2-exp(lambda1))/(exp(lambda2)-exp(lambda1));
b = exp(lambda2.*x)-exp(lambda1.*x);
c = C0.*(exp(lambda1.*x) + a.*b);

plot(x, c);


