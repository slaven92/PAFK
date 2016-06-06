clear all
clc
close all

pkg load signal

Fs = 1000;
T = 1/Fs;
L = 10000;
t = (0:L-1)*T;

%S=(chirp(t,0.5,20,10));
S=sign(chirp(t,0.5,20,10));
S=[S S S S S S S S S S];
L=10*L;
%X = S + 2*randn(size(t));
X = S;

Y = fft(X);
u=ifft(Y);
P2 = (Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;


plot(f,abs((P1)));
figure(2)
plot(X)