clear all;
close all;

%% Set path to ako's and roza's functions

rozita=genpath('C:\work\fuggvenyek_algoritmusok\ako_functions');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\saisir');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\matzler_Mie');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\roza_functions\');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\cmu\+cmu\');
addpath(rozita,'C:\work\PMMA_model\Measured_data\');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\geom2d\geom2d\');
addpath(rozita,'C:\work\');

%% load the measured pollen data, take the ripples from it, calculate the Mie

a=0.000012;
aa=a; % the space arround the sphere
G=(a*a*pi*4.0);
NMAX=76;
c=3e08;
hplanck=6.62607004081e-34;
mu=1.0;

n_res=3;

wn=(200:0.1:1500);% choose the theoretical wavenumber region  
%wn=(4219:1.0:4222);% choose the theoretical wavenumber region  
% wn=linspace(4207,4208,10); % electric mode 40 1st WGM (TM mode - a coeff
% wn=linspace(4221,4222,10); % magnetic mode 36 2nd WGM (TE mode - b coeff)
% wn=linspace(4268,4269,10); % electric mode 36 2nd WGM (TM mode - a coeff)
% wn=linspace(4216,4217,10); % magnetic mode 32 3nd WGM (TE mode - b coeff)
% wn=linspace(4247,4248,10); % electric mode 32 3nd WGM (TM mode - a coeff

k=2*pi*wn*100;          % is given in meter-1
sf=a.*k;                % size parameter

KVar=size(wn,2);

n0=1.46.*ones(1,KVar);
kk=0;

for n_i=[0 0.00004 0.0004 0.004]

    kk=kk+1;
    
    an=zeros(KVar,NMAX);
    bn=zeros(KVar,NMAX);
    cn=zeros(KVar,NMAX);
    dn=zeros(KVar,NMAX);
    
    Zn.d=n0(1,1:KVar).*ones(1,KVar);
    Zn.v=num2str(wn);
    Zn.i='Real part of ref. ind.';
    
    Znp.d=n_i.*ones(1,KVar);
    Znp.v=num2str(wn);
    Znp.i='Imaginary part of ref. ind.';
    
    ZRefIndexComplex.d=complex(Zn.d,Znp.d); % calculate the complex ref. ind
    ZRefIndexComplex.v=num2str(wn);
    ZRefIndexComplex.i='Complex refractive index';
    
    % calculate cn, dn, an, bn --- the Mie scattering coeficients
    
    for ii=1:KVar
        
        eff=Mie_abcd_vec(ZRefIndexComplex.d(1,ii),sf(1,ii));
       %  eff=mie_abcd(ZRefIndexComplex.d(1,ii),sf(1,ii));
        
        an(ii,1:size(eff,2))=eff(1,:);
        bn(ii,1:size(eff,2))=eff(2,:);
        cn(ii,1:size(eff,2))=eff(3,:);
        dn(ii,1:size(eff,2))=eff(4,:);
        
    end
    
    n_front=(2.*(1:NMAX)+1);
    A=n_front.*(real(an(:,1:NMAX)+bn(:,1:NMAX)));
    C_ext(kk,:)=((2*pi)./(k.*k)).*(sum(A,2))';

    clear A;
end

%% plot the exctinction

%legend('n_i=0','n_i=0.0000016','n_i=0.000016','n_i=0.00016','n_i=0.0016');

figure;
orient(gcf,'landscape');
plot(sf,C_ext(1,:)./(pi*a*a),'b','Linewidth',1.25);
%axis([4147 4268 2.25 2.6]);
set(gca,'XDir','Reverse','LineWidth', 1.25,'FontSize', 18);
xlabel('Size parameter [a.u.]','FontSize',18);
ylabel('Q_{ext}','FontSize',18);
legend({'n_i=0'},'Box','on','Linewidth',0.5);

figure;
orient(gcf,'landscape');
subplot(2,2,1);
plot(wn,C_ext(1,:)./(pi*a*a),'b','Linewidth',1.25);
%axis([4147 4270 2.25 2.6]);
set(gca,'XDir','Reverse','LineWidth', 1.25,'FontSize', 18);
xlabel('Wavenumber [cm^{-1}]','FontSize',18);
ylabel('Q_{ext}','FontSize',18);
legend({'n_i=0'},'Box','on','Linewidth',0.5);

subplot(2,2,2);
plot(wn,C_ext(2,:)./(pi*a*a),'b','Linewidth',1.25);
%axis([4147 4268 2.25 2.6]);
set(gca,'XDir','Reverse','LineWidth', 1.25,'FontSize', 18);
xlabel('Wavenumber [cm^{-1}]','FontSize',18);
ylabel('Q_{ext}','FontSize',18);
legend({'n_i=0.00004'},'Box','on','Linewidth',0.5);

subplot(2,2,3);
plot(wn,C_ext(3,:)./(pi*a*a),'b','Linewidth',1.25);
%axis([4147 4268 2.25 2.6]);
set(gca,'XDir','Reverse','LineWidth', 1.25,'FontSize', 18);
xlabel('Wavenumber [cm^{-1}]','FontSize',18);
ylabel('Q_{ext}','FontSize',18);
legend({'n_i=0.0004'},'Box','on','Linewidth',0.5);

subplot(2,2,4);
plot(wn,C_ext(4,:)./(pi*a*a),'b','Linewidth',1.25);
axis([4147 4268 2.25 2.6]);
set(gca,'XDir','Reverse','LineWidth', 1.25,'FontSize', 18);
xlabel('Wavenumber [cm^{-1}]','FontSize',18);
ylabel('Q_{ext}','FontSize',18);
legend({'n_i=0.004'},'Box','on','Linewidth',0.5);