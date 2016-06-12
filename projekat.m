clear all
clc
close all


ncl=3.17;
nco=3.43;
lambda=1550;
d=[200 290];

v=2.*pi./lambda.*d.*sqrt(nco^2-ncl^2);

temp=v.^2/2;
b=1-log(1+temp)./temp;

neff=sqrt(nco^2.*b + ncl^2.*(1-b));

gama=2.*pi./lambda.*sqrt(neff.^2-ncl.^2);

deff=d+2./gama;

kapa=4./pi.*2.*pi./lambda.*45.*(nco^2-neff.^2)./(2.*neff.*deff).*1000;

r=(neff(1)-neff(2))./(neff(1)+neff(2));
r=abs(r);

rg=tanh(log((1+r)./(1-r)));


%% drugi deo


la=500; % aktivna oblast u [um]
lp=75; % fazna oblast u [um]

lam=1540; % lambda za slucaj projekta
deltaLambda=[4 4.5];% 1=front 2=back ovo je spejsing izmedju pikova ogledala treba odabrati pravi [nm]
ng=3.8; % grupni indeks
envelopa=100; % propusni opseg za oba ogledala [nm]

z1=lam.^2./(2.*ng.*envelopa)./1000; % ovako se racuna z1 u [um]
z0=lam.^2./(2.*ng.*deltaLambda)./1000;% 1=front 2=back ovako se racuna z0 za oba ogledala [um]

%u slucjau da zelimo da zaokruzimo vrednosti za duzinu z0
%z0=round(z0);
%deltaLambda=lam.^2./(2.*3.8.*z0)./1000;
%--------------------------------------------------------
kapa1=kapa(2); % jacina sprege, pitanje je da li se koriti ovo [um^-1]

nbf=5;  % broj burstova prednjeg ogledala
nbb=10;  % broj burstova zadnjeg ogledala
%nbb=8:10;  % broj burstova zadnjeg ogledala

%back ogledalo
kuperb=nbb.*z1.*kapa1./((nbb-1).*z0(2)+z1);% efektivna jacina sprezanja [um-1]

leffb=1./2./kuperb.*tanh(kuperb.*((nbb-1).*z0(2)+z1)); % efekivna duzina [um]

fwhmb=lam^2./(2.*ng.*leffb).*3.76./pi./1000; %fwhm [nm]
%=======================================================================

%front ogledalo
kuperf=nbf.*z1.*kapa1./((nbf-1).*z0(1)+z1);% efektivna jacina sprezanja [um-1]

lefff=1./2./kuperf.*tanh(kuperf.*((nbf-1).*z0(1)+z1));% efekivna duzina [um]

fwhmf=lam^2./(2.*ng.*lefff).*3.76./pi./1000; %fwhm [nm]
%======================================================

%figure(1)

%plot(nbb,fwhmb,nbb,deltaLambda(2).*ones(1,length(nbb)),nbf,fwhmf,nbf,deltaLambda(1).*ones(1,length(nbf)))

%% racunanje razmak modova cele supljine


lc=leffb'*ones(1,length(nbf)) + ones(length(nbb),1)*lefff + la + lp;%duzina supljine od broja burstova

deltamode=lam.^2./(2.*ng.*lc)./1000; % spejsing modova cele supljine [nm]

smaler=zeros(length(nbb),length(nbf)); % ovde se stavlja manji fwhm od oba ogledala [nm]

for i=1:length(nbb)
  smaler(i,:)=min(fwhmb(i),fwhmf);
end


Rbprim=tanh(kapa1.*nbb.*z1.^2);
Rfprim=tanh(kapa1.*nbf.*z1).^2;
Rb=(Rbprim).*exp(-2.*5.*leffb.*1e-4);
Rf=(Rfprim).*exp(-2.*5.*lefff.*1e-4);


t=smaler<(3.*deltamode); % uslov da fwhm mora da bude manji od 2*rastojanje modova cele supljine
a=find(t!=0); % sadrzi mesta na kojima se nalaze slu;ajevi koji zadovoljavaju uslov

%figure(2)
%hold all
%surf(2.*deltamode-smaler)
%surf(2.*deltamode)
%surf(smaler)
%hold off



%%struje ogledala potrebne za tjunovanje
nbb1=10;
nbf1=5;
lc1=750; % uneti izabranu ukupnu duzinu na osnovu broja burstova iz prethodnog dela [um]
leff1=lefff; % popuniti odozgo za prednje ogledalo [um]
leff2=leffb; % popuniti odozgo za zadnje ogledalo [um]
gama1=la./lc1.*0.043; % faktor konfiniranja 
ni=0.7; % gubici ubacivanja struje
q=1.6e-19; % nalektrisanje elektrona [C]
V=0.007.*la.*7.*3.*1e-12; % zapremina aktivne oblasti [cm3] 
B=0.5e-10; %[cm3/s]


dl=[1 2 3 0 3 2 1; 3 3.5 4 0 1.5 1 0.5]; % pomeraji za odg tal duz; gore prednje, dole zadnje[nm]
neff1=3.2609; % efektivni indeks prelamanja upisati ili ovaj ili onaj iz prvog dela zadatka
lam1=lam; % upisati lambda na kojoj ce se raditi
deltan=-neff1.*dl./lam1; % pomeraj indeska za koji se ostvaruje pomeraj dl

konc=deltan./(0.01.*gama1).*1e18; %koncentracija koja je potrebna da se dobije dn[cm-3]

C=(-7.*neff1./lam1+1)./(log10(90));
%current=q.*V.*B.*konc.^2./ni; % struja koja je potrebna za dl[A]
current=10.^((deltan+1)./C)-10; % struja koja je potrebna za dl[mA]

%%treshhold za postizanje laserovanja
%rf=tanh(kapa1.*nbf1.*z1);
%rb=tanh(kapa1.*nbb1.*z1);

r=[0.1 0.2 0.3 0.4 0.5 0.6 0.7; 0.1 0.2 0.3 0.4 0.5 0.6 0.7]; % front gore back dole


%alphap=5+5.*konc.*1e-18; %gubici u ogledalima i fazi[cm-1]
alphap=5; %gubici u ogledalima i fazi[cm-1]
alphaa=15;% gubici u akt [cm-1]
alphai=(alphap.*(lp+leff1+leff2)+alphaa.*la)./lc1; %srednji gubici [cm-1]

alpham=1./lc1.*log(1./(r(1,:).*r(2,:))).*10000; % gubici ogledala [cm-1]



th=7.*3.*la.*1e-8./ni.*60.*exp((alphai+alpham)./(gama1.*650)); %struja praga [A]

%%
nd=ni.*alpham./(alphai+alpham);

hnq=1.98644568e-25./(lam1.*1e-9.*1.60217662e-19);
struja=150E-3;%[A]

pow=nd.*hnq.*(struja-th).*1000; %[mW]