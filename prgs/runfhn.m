clear

function xdot = fhn(x,t)
	global ALPHA BETA;
	
	xdot(1) = x(1) - x(1)^3 - x(2) + ALPHA;
	xdot(2) = BETA * (x(1) - x(2));
	
endfunction	

#############
#
D1 = 1; D2 = 10*D1; 
alpha = (1/10)^(3/2); beta  = 1;

L = 7.3*10;
ngrid = 100;
dt = 1e-4;
nloops = 10000000;

c1 = alpha^(1/3).*ones(1,ngrid)+ 0.1.*randn(1,ngrid);
c2 = alpha^(1/3).*ones(1,ngrid);
x = linspace(0, L, ngrid);

for n=1:nloops
	r1 = c1 - c1.^3 - c2 + alpha;
	r2 = beta.*(c1 - c2);

	[dum diffc1] = der1d(c1, x, 'n', 'c');
 	[dum diffc2] = der1d(c2, x, 'n', 'c');
 
	c1 = c1 + (r1 + D1.*diffc1).*dt;
	c2 = c2 + (r2 + D2.*diffc2).*dt; 	
	
	c1(1)=c1(2); c1(end)=c1(end-1);
	c2(1)=c2(2); c2(end)=c2(end-1);

	if rem(n,1000)==0
		plot(x, c1, x, c2);
		 pause(0.01);
	endif
endfor

A = [x', c1'];
save -ascii tmp.dat A
 
return

#############
D1 = 1;
D2 = 10.1*D1; 

k = linspace(0, 1, 200); 

alpha = (1/10)^(3/2);
beta  = 1;

a = -7/10+k.^2*D1;
b = beta + k.^2*D2;


lambda1 = 0.5*(-(a+b) + sqrt((a+b).^2-4.*(beta+a.*b)));
lambda2 = 0.5*(-(a+b) - sqrt((a+b).^2-4.*(beta+a.*b)));

return
##############
#
global ALPHA BETA;

ALPHA = (0.1)^(1.5); BETA = 7/10 + 0.1; 
fpoint = ALPHA^(2/3);

t = linspace(0,200,10000);
x = lsode("fhn", [fpoint, fpoint+0.1], t);

plot(x(:,1), x(:,2))
hold on;
plot(ALPHA^(1/3), ALPHA^(1/3), 'o');
hold off;

return;

###
alpha=1e-3; 
beta = linspace(0, 2);
b = beta - 1 + 3*alpha.^(2/3); 
c = 3*beta.*alpha.^(2/3); 

lambda_p =  0.5*(-b + sqrt(b.^2-4.*c)); 
lambda_m =  0.5*(-b - sqrt(b.^2-4.*c)); 
#plot(beta, lambda_p, beta, real(lambda_m)) 


###
alpha=1e-2; 
beta = 2;
D1 = 1;
D2 = 1;
L = 100;
a = alpha^(1/3);

for n=1:100

	k=n*pi/L;

	b =  beta + k^2*D2 -1 +3*a^2 + k^2*D1;
	c = (beta + k^2*D2)*(-1+3*a^2+k^2*D1) + beta;

	coeff=[1, b, c];

	ret(n,:)=roots(coeff)';
	
	karray(n) =k;
end

plot(karray, real(ret(:,1)), karray, real(ret(:,2)));

find(real(ret(:,1))>0.0)
find(real(ret(:,2))>0.0)
