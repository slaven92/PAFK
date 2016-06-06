% clear all
% close all
% clc

ng = 4;
lambda0 = 1315;
Lb= 4500;

n1 = 3.43; n2 = 3.17; 
d1 = 200; d2 = 290;
V01 = 2*pi/lambda0*d1*sqrt(n1^2-n2^2);
b01 = 1 - log(1+V01^2/2)/(V01^2/2);
neff01 = sqrt(b01*n1^2+(1-b01)*n2^2);
V02 = 2*pi/lambda0*d2*sqrt(n1^2-n2^2);
b02 = 1 - log(1+V02^2/2)/(V02^2/2);
neff02 = sqrt(b02*n1^2+(1-b02)*n2^2);

L1 = round(lambda0/(4*neff01)); L2 = round(lambda0/(4*neff02));
m = round(Lb/(L1+L2));
L3 = round(43000-Lb); N = 10;

lambda = 1280:0.001:1350;
V1 = 2.*pi./lambda.*d1.*sqrt(n1.^2-n2.^2);
b1 = 1 - log(1+V1.^2./2)./(V1.^2./2);
neff1 = sqrt(b1.*n1.^2+(1-b1).*n2.^2);
V2= 2.*pi./lambda.*d2.*sqrt(n1.^2-n2.^2);
b2 = 1 - log(1+V2.^2./2)./(V2.^2./2);
neff2 = sqrt(b2.*n1.^2+(1-b2).*n2.^2);

alpha = 5*(1/10^7);
beta1 = 2.*pi./lambda.*neff1-1i.*alpha./2; 
beta2 = 2.*pi./lambda.*neff2-1i.*alpha./2;
r = (neff2 - neff1)./(neff2 + neff1);
t = sqrt(1-r.^2);

delta_lambda_B =0;
delta_neff1 = delta_lambda_B/(4*L1); delta_neff2 = delta_lambda_B/(4*L2);

phi_plus = beta1.*(lambda0+delta_lambda_B)./(4*neff01+delta_neff1) + ...
    beta2.*(lambda0+delta_lambda_B)./(4*neff02+delta_neff2);
phi_minus = beta1.*(lambda0+delta_lambda_B)./(4*neff01+delta_neff1) - ...
    beta2.*(lambda0+delta_lambda_B)./(4*neff02+delta_neff2);

T11 = 1./(t.^2).*(exp(1i.*phi_plus)-r.^2.*exp(-1i.*phi_minus));
T21 = r./(t.^2).*(exp(1i.*phi_plus)-exp(-1i.*phi_minus));
T12 = T21;
T22 = 1./t.^2.*(exp(-1i.*phi_plus)-r.^2.*exp(1i.*phi_minus));

ksi = log(0.5.*(T11+T22)+sqrt(0.25.*(T11+T22).^2-1));
Delta = 1i.*(T22-T11)./(T22+T11);
meff = tanh(m.*ksi)./tanh(ksi);

Tg11 = (1+1i.*meff.*Delta).*cosh(m.*ksi);
Tg21 = T21./T11.*meff.*(1+1i.*Delta).*cosh(m.*ksi);
Tg12 = T12./T22.*meff.*(1-1i.*Delta).*cosh(m.*ksi);
Tg22 = (1-1i.*meff.*Delta).*cosh(m.*ksi);

for i=1:length(lambda)
    Tg = [Tg11(i) Tg12(i); Tg21(i) Tg22(i)];
    Tb = [exp(1i*beta2(i)*L3) 0; 0 exp(-1i*beta2(i)*L3)];
    Tsgdbr1 = (Tb*Tg)^(N-1); %N ili N-1
    Tsgdbr = Tg*Tsgdbr1;
    Reflectivity(i) = (abs(Tsgdbr(2,1)/Tsgdbr(1,1)))^2;
end;

r_bar = (neff02 - neff01)/(neff02 + neff01);
kappa = 2*r_bar/(L1+L2);

LF = (4-1)*45000 + 4000 ;
LB = (10-1)*43000 + 4500;
kappa_avg_F = 4*4000*kappa/LF;
kappa_avg_B = 10*4500*kappa/LB;
LeffF = (1/(2*kappa_avg_F))*tanh(kappa_avg_F*LF);
LeffB = (1/(2*kappa_avg_B))*tanh(kappa_avg_B*LB);
La = 500000; Lp = 75000;
Lc = La+Lp+LeffF+LeffB;
delta_lambda_mode = lambda0^2/(2*ng*Lc);

limit = 1280+round((1350-1280)/delta_lambda_mode)*delta_lambda_mode;
mode_comb = 1280:delta_lambda_mode:limit;
for i = 1:length(mode_comb)
    peak(i) = 0.05;
end;

figure(1)
hold all
plot(lambda,Reflectivity);
plot(lambda,RF_45_4_4);
stem(mode_comb,peak,'.','MarkerSize',0.1);
ax = gca;
ax.XTick = 1280:5:1350;
xlim([1280 1350]);
ylim([0 1]);
grid on

