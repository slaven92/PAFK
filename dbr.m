clear all
close all
clc
set (0, "defaultaxesfontname", "/usr/share/fonts/truetype/msttcorefonts/arial.ttf")
set (0, "defaultaxesfontsize", 10)
set (0, "defaulttextfontname", "arial")
set (0, "defaulttextfontsize", 5) 


lamda0=1550; % nanometri
n1=3.43;
n2=3.17;
r=abs((n1-n2)/(n1+n2));
t=2*sqrt(n1*n2)/(n1+n2);

delta=300;
N=1000;
m=50;
L1=lamda0./(4.*n1);
L2=lamda0./(4.*n2);
lmbda=linspace(lamda0-delta,lamda0+delta,N);

beta1=2.*pi.*n1./lmbda;
beta2=2.*pi.*n2./lmbda;

phiplus=beta1.*L1+beta2.*L2;
phiminus=beta1.*L1-beta2.*L2;

matrix=zeros(2,2,N);
matrix(1,1,:)=(exp(j.*phiplus)-r^2.*exp(-j.*phiminus))./t.^2;
matrix(1,2,:)=r.*(exp(j.*phiplus)-exp(-j.*phiminus))./t.^2;
matrix(2,1,:)=r.*(exp(-j.*phiplus)-exp(j.*phiminus))./t.^2;
matrix(2,2,:)=(exp(-j.*phiplus)-r^2.*exp(j.*phiminus))./t.^2;

for i=1:N
matrix(:,:,i)=matrix(:,:,i)^m;
end



reflection=abs(matrix(2,1,:)./matrix(1,1,:));
phase=unwrap(angle(matrix(2,1,:)./matrix(1,1,:)));

reflection=reshape(reflection,1,[]);
phase=reshape(phase,1,[]);

plot(lmbda,reflection.^2)
xlabel('aaaa');
ylabel('aaaa');
legend(['lala'])
figure(2)
plot(lmbda,phase)

