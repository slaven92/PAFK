clear all
close all
clc

lamda0=1540; % nanometri
ni2=3.2609;
ni1=3.2255;
%deltan=ni2.*1./lamda0; %radi po formuliiiiiiii
deltan=0; %radi po formuliiiiiiii
n2=3.2609+deltan;
n1=3.2255+deltan;
r=abs((n1-n2)/(n1+n2));
t=2*sqrt(n1*n2)/(n1+n2);

delta=20;
N=5000;
bursts=10;
z0=50000;
z1=3000;
%L1=lamda0./(4.*n1);
L1=lamda0./(4.*ni1);
%L2=lamda0./(4.*n2);
L2=lamda0./(4.*ni2);
lmbda=linspace(lamda0-delta,lamda0+delta,N);
m=round(z1/(L1+L2));
z1=m.*(L1+L2);
L3=round((z0-z1)./(lamda0./(2*ni2)));
L3=L3.*lamda0./(2*ni2);

beta1=2.*pi.*n1./lmbda;
beta2=2.*pi.*n2./lmbda;
beta3=beta2;

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
matrix2=matrix; %treba nam da bi dodali na kraj jos jedan ovakav

%% pomnoziti matrix sa refl(nije potreban) i transl
%refl=ones(2,2,N);
%refl(1,2,:)=r;
%refl(2,1,:)=r;
%refl=refl./t;

transm=zeros(2,2,N);
transm(1,1,:)=exp(j.*beta3.*L3);
%transm(1,2,:)=exp(-j.*beta3.*L3);
%transm(2,1,:)=exp(j.*beta3.*L3);
transm(2,2,:)=exp(-j.*beta3.*L3);

for i=1:N
%matrix(:,:,i)=(matrix(:,:,i)*refl(:,:,i))*transm(:,:,i);
matrix(:,:,i)=matrix(:,:,i)*transm(:,:,i);
end
%%----------------------------------------


%% podici matrix na burst ekvivalentno sa N burstova
for i=1:N
matrix(:,:,i)=matrix(:,:,i)^bursts;
end
%%-------------------------



%% pomnoziti sa jos jednom matrix na m ekv sa dodavanjem jos m resetaka na kraj;
for i=1:N
matix(:,:,i)=matrix(:,:,i)*matrix2(:,:,i);
end
%%====================================

reflection=abs((matrix(2,1,:)./matrix(1,1,:)));
phase=unwrap(angle(matrix(2,1,:)./matrix(1,1,:)));

reflection=reshape(reflection,1,[]);
phase=reshape(phase,1,[]);

plot(lmbda,reflection.^2) 
%figure(2)
%plot(lmbda,phase)
