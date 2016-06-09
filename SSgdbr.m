clear all
close all
clc

lamda0=1550; % nanometri
n1=3.2255;
n2=3.26;
r=abs((n1-n2)/(n1+n2));
t=2*sqrt(n1*n2)/(n1+n2);

delta=200;
N=2000;
m=5;
lmbda=linspace(lamda0-delta,lamda0+delta,N);

L1=lamda0./(4.*n1);
L2=lamda0./(4.*n2);
N1=100;
ma=100;
duzina=65000;
dn=(duzina-N1.*(L1+L2))./((N1-1)*N1);
%dn=(ma.*lamda0./2 - (L1.*n1+n2.*L2).*N1./2)./((N1.*(N1-1).*(n1+n2))./2);
Lovi=L1;
Lovi(2)=L2;

for i=3:2:N1
  Lovi(i)=Lovi(i-2)+dn;
  Lovi(i+1)=Lovi(i-1)+dn;
end

sum(Lovi)
%Lovi=linspace(2000,3000,100);
%Lovi=[L1 L2];
beta1=2.*pi.*n1./lmbda;
beta2=2.*pi.*n2./lmbda;

%L1=5000;
%N1=20;
%dd=500;
%dn=(L1.*(n1+n2).*N1-dd.*lamda0./2)./(n1.*(N1-1).*(N1-1)+n2.*N1.^2);
%Lovi=linspace(L1,L1-(2.*N1-1).*dn,2.*N1)




for i=1:2:length(Lovi)
  phiplus=beta1.*Lovi(i)+beta2.*Lovi(i+1);
  phiminus=beta1.*Lovi(i)-beta2.*Lovi(i+1);
  matrix=zeros(2,2,N);
  matrix(1,1,:)=(exp(j.*phiplus)-r^2.*exp(-j.*phiminus))./t.^2;
  matrix(1,2,:)=r.*(exp(j.*phiplus)-exp(-j.*phiminus))./t.^2;
  matrix(2,1,:)=r.*(exp(-j.*phiplus)-exp(j.*phiminus))./t.^2;
  matrix(2,2,:)=(exp(-j.*phiplus)-r^2.*exp(j.*phiminus))./t.^2;
  
  if(i==1)
    out=matrix;
  else
  
    for k=1:length(lmbda)
      out(:,:,k)=out(:,:,k)*matrix(:,:,k);
    end
  
  end
end

matrix=out;

for i=1:N
matrix(:,:,i)=matrix(:,:,i)^m;
end

reflection=abs(matrix(2,1,:)./matrix(1,1,:));
reflection=reshape(reflection,1,[]);
plot(lmbda,reflection)