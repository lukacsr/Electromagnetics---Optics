% ***********************************************************************
%     3-D code for calculation of the electric field inside and outside of
%     a sphere
% ***********************************************************************
%
%     Program author: Rozalia Lukacs
%
%     First version: 23.01.2018
%
%     Version: 09. 07. 2018
%
%     This MATLAB M-file implements the analytical solutions of the 
%     scattering of a plane wave on a dielectric field.
%     It gives both the three dimensional eletric field and magnetic 
%     field inside and outside of the sphere.
%     It uses the spherical harmonics representation of the fields
%     following Bohren and Huffman.
%
%     The program can plot in three dimensions the electric field.
%     It also plots the 2 dimensional plots of the electric and
%     magnetic fields.
%
% ***********************************************************************

clear all;
close all;

%% Set path to ako's and roza's functions

rozita=genpath('C:\work\fuggvenyek_algoritmusok\ako_functions');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\roza_functions\');
addpath(rozita,'C:\work\fuggvenyek_algoritmusok\cmu\+cmu\');
addpath(rozita,'C:\work\PMMA_model\Measured_data\');
addpath(rozita,'C:\work\');

%% load the measured pollen data, take the ripples from it, calculate the Mie

a=0.000012; % radius of the sphere
aa=a; % the space arround the sphere for 3rd WGM 32
%aa=a/2; % the space arround the sphere for 2nd WGM 36
%aa=a/4; % the space arround the sphere for 1st WGM 40
G=(a*a*pi*4.0);
NMAX=66;
c=3e08; % speed of light
hplanck=6.62607004081e-34; % plank constant
mu=1.0;

Nx=240;   % grid points on Ox axis
Ny=240;   % grid points on Oy axis
Nz=240;   % grid points on Oy axis

n_res=5;

kor_a=Nx*(1/4); % the red disk arround 3rd WGM 32
%kor_a=Nx*(2/6); % the red disk arround 2nd WGM 36
%kor_a=Nx*(4/10); % the red disk arround 1st WGM 40

% wn=linspace(512,518,10);% choose the theoretical wavenumber region
% wn=linspace(4147,4148,10);% choose the theoretical wavenumber region                 % wavenumber in cm-1
wn=linspace(4207,4208,10); % electric mode 40 1st WGM
% wn=linspace(4242,4243,10); % electric mode 41 1st WGM
% wn=linspace(4167,4168,10); % electric mode 35 2nd WGM
% wn=linspace(4221,4222,10); % magnetic mode 36 2nd WGM
% wn=linspace(4268,4269,10); % electric mode 36 2nd WGM
% wn=linspace(4251,4252,10); % electric mode 32 3nd WGM
% wn=linspace(4221,4222,10); % magnetic mode 32 3nd WGM

k=2*pi*wn*100;          % is given in meter-1
sf=a.*k;                % size parameter

KVar=size(wn,2);

n0=1.46.*ones(1,KVar);

mode_n=40; % the order of the Whsipering Gallery Mode

an=zeros(KVar,NMAX);
bn=zeros(KVar,NMAX);
cn=zeros(KVar,NMAX);
dn=zeros(KVar,NMAX);

Zn.d=n0(1,1:KVar).*ones(1,KVar);
Zn.v=num2str(wn);
Zn.i='Real part of ref. ind.';

Znp.d=0.0*ones(1,KVar);
Znp.v=num2str(wn);
Znp.i='Imaginary part of ref. ind.';

ZRefIndexComplex.d=complex(Zn.d,Znp.d); % calculate the complex ref. ind
ZRefIndexComplex.v=num2str(wn);
ZRefIndexComplex.i='Complex refractive index';

% calculate cn, dn, an, bn --- the Mie scattering coeficients

for ii=1:KVar
    
    eff=Mie_abcd_vec(ZRefIndexComplex.d(1,ii),sf(1,ii));
    
    an(ii,1:size(eff,2))=eff(1,:);
    bn(ii,1:size(eff,2))=eff(2,:);
    cn(ii,1:size(eff,2))=eff(3,:);
    dn(ii,1:size(eff,2))=eff(4,:);
    
end

%% calculate the electric field inside the sphere for a given n and Phi

% define the grid

NN=Nx*Ny*Nz;

Lx=2*a+2*aa;
Ly=2*a+2*aa;
Lz=2*a+2*aa;

dx=Lx/Nx; % spacing in x-direction
dy=Ly/Ny; % spacing in y-direction
dz=Ly/Nz; % spacing in z-direction

x=zeros(1,NN);
y=zeros(1,NN);
z=zeros(1,NN);

% in order to understand the indexes, check my notes from Wesleyan stay
% in the topic, internal/scattered fields

indexx=repmat(repelem(Nx:-1:1,Ny),1,Nz);

if mod(Ny,2)==0
    indexy=repmat(repmat([1:Ny,Ny:-1:1],1,floor(Nx/2)),1,Nz);
else
    indexy=repmat([repmat([1:Ny,Ny:-1:1],1,floor(Nx/2)),1:Ny],1,Nz);
end

indexz=repelem(1:Nz,Nx*Ny);
indexl=1:Nx*Ny*Nz;

x(indexl)=((indexx-0.5).*dx)-a-aa; % the Descartes coordinate x
y(indexl)=((indexy-0.5).*dy)-a-aa; % the Descartes coordinate y
z(indexl)=((indexz-0.5).*dz)-a-aa; % the Descartes coordinate z

% transfer the cartesian coordiantes to polar coordianates

rhop=sqrt(x.*x + y.*y + z.*z);      % the r coordinate

for kk=1:NN                         % the phi coordinate
    
    if (y(kk)<0)
        phi(kk)=atan2(y(kk),x(kk))+2*pi;
    else
        phi(kk)=atan2(y(kk),x(kk));
    end
    
    % the theta coordinate
    thetap(kk)=atan2(sqrt(x(kk).*x(kk)+y(kk).*y(kk)),z(kk));
    
end

dtheta=abs(thetap(NN-1)-thetap(NN));
dphi=abs(phi(NN-Nx)-phi(NN));
dro=abs(rhop(NN-1)-rhop(NN));

% calculate the cos and sin of theta and phi

koszphi=cos(phi);
szinphi=sin(phi);

koszthetap=cos(thetap);
szinthetap=sin(thetap);

rhopkor = rhop < a; % the point is inside the circle rhopkor takes value 1
% if is outside then then takes value 0

rhopout = rhop >= a; % the point is outside the circle rhopout takes value 1
% if is inside then takes value 0

Area_section=pi*a*a;

E1=zeros(KVar,Nx*Ny*Nz);
H1=zeros(KVar,Nx*Ny*Nz);
Es=zeros(KVar,Nx*Ny*Nz);
Hs=zeros(KVar,Nx*Ny*Nz);

AA1=zeros(KVar,4);
AA2=zeros(KVar,4);
AA3=zeros(KVar,4);

for n=mode_n:mode_n
    
    % calculate i^2*n
    
    if mod(n,2)==0
        in2=1;
    end
    
    if mod(n,2)==1
        in2=-1;
    end
    
    % calculate prefactor from n
    
    nu=(n+(1/2)); %  Bohren-Huffman p.86
    numinus=(nu-1);
    nuplus=(nu+1);
    n2=(n^2.*(n+1)^2);
    nn=2*n+1;
    n22=(nn)^2./(n2);
    
    for ii=1:KVar
        
        % here I need k1, which is the wavenumber inside the
        % sphere and not the k, which is the wavenumber outside of the sphere!
        % it is assumed that the outside medium is air. The program must be
        % modified for a different outside medium
        
        % give E0
        
        %E0=(hplanck/(2*pi))*c*k(ii);
        E0=1.0;
        
        [pe, te]=Mie_pt_vec(koszthetap,n);
        
        tau_n(1,:)=te(:,n)';
        pi_n(1,:)=pe(:,n)';
        tau_n2=tau_n.*tau_n;
        pi_n2=pi_n.*pi_n;
        
        % calculation of the spherical bessel functions
        
        k1=k(ii).*ZRefIndexComplex.d(1,ii);
        Z=k1.*rhop;
        
        sqZ= sqrt((0.5*pi)./Z);
        
        Jn=besselj(nu, Z.*rhopkor).*sqZ;
        Jn_minus=besselj(numinus, Z.*rhopkor).*sqZ;
        Jn_plus=besselj(nuplus, Z.*rhopkor).*sqZ;
        Jn_per_Z=Jn./Z;
        Dj=((Jn_per_Z)+((nu/(2*nu+1)).*Jn_minus)-(((nu+1)/(2*nu+1)).*Jn_plus));
        
        % calculate the vector spherical harmonics for generating function
        % Jnu - inside the sphere
        
        Mo1n=[koszphi.*pi_n.*Jn; -szinphi.*tau_n.*Jn;];
        
        Ne1n=[koszphi.*n.*(n+1).*szinthetap.*pi_n.*Jn_per_Z;  koszphi.*tau_n.*Dj;  -szinphi.*pi_n.*Dj;];
        
        Me1n=[-szinphi.*pi_n.*Jn;  -koszphi.*tau_n.*Jn;];
        
        No1n=[szinphi.*n.*(n+1).*szinthetap.*pi_n.*Jn_per_Z;  szinphi.*tau_n.*Dj;  koszphi.*pi_n.*Dj;];
        
        Mo1n_cc=Mo1n;
        Me1n_cc=Me1n;
        No1n_cc=No1n;
        Ne1n_cc=Ne1n;
        
        % check the ortognality of the vector sperical harmonics
        
        AA1(ii,1)=sum((Mo1n_cc(1,:).*Me1n(1,:)+Mo1n_cc(2,:).*Me1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA1(ii,2)=sum((Me1n_cc(1,:).*Mo1n(1,:)+Me1n_cc(2,:).*Mo1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA1(ii,3)=sum((Me1n_cc(1,:).*Me1n(1,:)+Me1n_cc(2,:).*Me1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA1(ii,4)=sum((Mo1n_cc(1,:).*Mo1n(1,:)+Mo1n_cc(2,:).*Mo1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA2(ii,1)=sum((Ne1n_cc(1,:).*No1n(1,:)+Ne1n_cc(2,:).*No1n(2,:)+Ne1n_cc(3,:).*No1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA2(ii,2)=sum((No1n_cc(1,:).*Ne1n(1,:)+No1n_cc(2,:).*Ne1n(2,:)+No1n_cc(3,:).*Ne1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA2(ii,3)=sum((No1n_cc(1,:).*No1n(1,:)+No1n_cc(2,:).*No1n(2,:)+No1n_cc(3,:).*No1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA2(ii,4)=sum((Ne1n_cc(1,:).*Ne1n(1,:)+Ne1n_cc(2,:).*Ne1n(2,:)+Ne1n_cc(3,:).*Ne1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA3(ii,1)=sum((0.*No1n(1,:)+Me1n_cc(1,:).*No1n(2,:)+Me1n_cc(2,:).*No1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA3(ii,2)=sum((0.*No1n_cc(1,:)+Me1n(1,:).*No1n_cc(2,:)+Me1n(2,:).*No1n_cc(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA3(ii,3)=sum((0.*Ne1n(1,:)+Mo1n_cc(1,:).*Ne1n(2,:)+Mo1n_cc(2,:).*Ne1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA3(ii,4)=sum((0.*Ne1n_cc(1,:)+Mo1n(1,:).*Ne1n_cc(2,:)+Mo1n(2,:).*Ne1n_cc(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        
        Vcell=(Lx*Ly*Lz);
        
        AA1(ii,:)=AA1(ii,:)./Vcell;
        AA2(ii,:)=AA2(ii,:)./Vcell;
        AA3(ii,:)=AA3(ii,:)./Vcell;
        
        % I calculate the absolute square of the electric field at k1 wavenumber
        % and given n and Phi inside the sphere
        
        E1r=    0               - 1i.*dn(ii,n).*Ne1n(1,:);
        
        E1t=cn(ii,n).*Mo1n(1,:) - 1i.*dn(ii,n).*Ne1n(2,:);
        
        E1p=cn(ii,n).*Mo1n(2,:) - 1i.*dn(ii,n).*Ne1n(3,:);
        
        E1(ii,:)=E0.*E0.*in2.*n22.*(abs(E1r).*abs(E1r)+abs(E1t).*abs(E1t)+abs(E1p).*abs(E1p));
        
        % I calculate the absolute square of the magnetic field at k1 wavenumber
        % and given n an Phi
        
        H1r=    0               + 1i.*cn(ii,n).*No1n(1,:);
        
        H1t=dn(ii,n).*Me1n(1,:) + 1i.*cn(ii,n).*No1n(2,:);
        
        H1p=dn(ii,n).*Me1n(2,:) + 1i.*cn(ii,n).*No1n(3,:);
        
        scaling=(k1./(c*k(ii)*mu)).*(k1./(c*k(ii)*mu));
        
        H1(ii,:)=scaling.*E0.*E0.*in2.*n22.*(abs(H1r).*abs(H1r) + abs(H1t).*abs(H1t) + abs(H1p).*abs(H1p));
        
        % calculation of the Hankel functions of first kind
        
        clear Me1n Mo1n Ne1n No1n;
        
        Ro=k(ii).*rhop;
        sqRo= sqrt((0.5*pi)./Ro);
        
        Hn1=besselh(nu,1,Ro).*sqRo.*rhopout;
        Hn1_minus=besselh(numinus,1,Ro).*sqRo.*rhopout;
        Hn1_plus=besselh(nuplus,1,Ro).*sqRo.*rhopout;
        Hn1_per_Ro=Hn1./Ro;
        Dhn1=((Hn1./(Ro))+((nu/(2*nu+1)).*Hn1_minus)-(((nu+1)/(2*nu+1)).*Hn1_plus));
        
        % calculate the vector spherical harmonics for generating function
        % Hnu1 - outside the sphere
        
        Mo1n=[koszphi.*pi_n.*Hn1; -szinphi.*tau_n.*Hn1;];
        
        Ne1n=[koszphi.*n.*(n+1).*szinthetap.*pi_n.*Hn1_per_Ro;  koszphi.*tau_n.*Dhn1;  -szinphi.*pi_n.*Dhn1;];
        
        Me1n=[-szinphi.*pi_n.*Hn1;  -koszphi.*tau_n.*Hn1;];
        
        No1n=[szinphi.*n.*(n+1).*szinthetap.*pi_n.*Hn1_per_Ro;  szinphi.*tau_n.*Dhn1;  koszphi.*pi_n.*Dhn1;];
        
        Mo1n_cc=[rot90(Mo1n(1,:)'); rot90(Mo1n(2,:)');];
        Me1n_cc=[rot90(Me1n(1,:)'); rot90(Me1n(2,:)');];
        No1n_cc=[rot90(No1n(1,:)'); rot90(No1n(2,:)'); rot90(No1n(3,:)');];
        Ne1n_cc=[rot90(Ne1n(1,:)'); rot90(Ne1n(2,:)'); rot90(Ne1n(3,:)');];
        
        % check the ortognality of the vector sperical harmonics
        
        AA4(ii,1)=sum((Mo1n_cc(1,:).*Me1n(1,:)+Mo1n_cc(2,:).*Me1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA4(ii,2)=sum((Me1n_cc(1,:).*Mo1n(1,:)+Me1n_cc(2,:).*Mo1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA4(ii,3)=sum((Me1n_cc(1,:).*Me1n(1,:)+Me1n_cc(2,:).*Me1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA4(ii,4)=sum((Mo1n_cc(1,:).*Mo1n(1,:)+Mo1n_cc(2,:).*Mo1n(2,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA5(ii,1)=sum((Ne1n_cc(1,:).*No1n(1,:)+Ne1n_cc(2,:).*No1n(2,:)+Ne1n_cc(3,:).*No1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA5(ii,2)=sum((No1n_cc(1,:).*Ne1n(1,:)+No1n_cc(2,:).*Ne1n(2,:)+No1n_cc(3,:).*Ne1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA5(ii,3)=sum((No1n_cc(1,:).*No1n(1,:)+No1n_cc(2,:).*No1n(2,:)+No1n_cc(3,:).*No1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA5(ii,4)=sum((Ne1n_cc(1,:).*Ne1n(1,:)+Ne1n_cc(2,:).*Ne1n(2,:)+Ne1n_cc(3,:).*Ne1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA6(ii,1)=sum((0.*No1n(1,:)+Me1n_cc(1,:).*No1n(2,:)+Me1n_cc(2,:).*No1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA6(ii,2)=sum((0.*No1n_cc(1,:)+Me1n(1,:).*No1n_cc(2,:)+Me1n(2,:).*No1n_cc(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA6(ii,3)=sum((0.*Ne1n(1,:)+Mo1n_cc(1,:).*Ne1n(2,:)+Mo1n_cc(2,:).*Ne1n(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        AA6(ii,4)=sum((0.*Ne1n_cc(1,:)+Mo1n(1,:).*Ne1n_cc(2,:)+Mo1n(2,:).*Ne1n_cc(3,:)).*rhop.*rhop.*szinthetap.*dtheta.*dphi.*dro);
        
        AA4(ii,:)=AA4(ii,:)./Vcell;
        AA5(ii,:)=AA5(ii,:)./Vcell;
        AA6(ii,:)=AA6(ii,:)./Vcell;
        
        % I calculate the absolute square of the electric field at k wavenumber
        % and given n and Phi outside the sphere
        
        Esr=1i.*an(ii,n).*Ne1n(1,:)  - 0 ;
        
        Est=1i.*an(ii,n).*Ne1n(2,:) - bn(ii,n).*Mo1n(1,:);
        
        Esp=1i.*an(ii,n).*Ne1n(3,:) - bn(ii,n).*Mo1n(2,:);
        
        Es(ii,:)=E0.*E0.*in2.*n22.*(abs(Esr).*abs(Esr)+abs(Est).*abs(Est)+abs(Esp).*abs(Esp));
        
        % I calculate the absolute square of the magnetic field at k wavenumber
        % and given n an Phi outside the sphere
        
        Hsr=1i.*bn(ii,n).*No1n(1,:) + 0 ;
        
        Hst=1i.*bn(ii,n).*No1n(2,:) + an(ii,n).*Me1n(1,:);
        
        Hsp=1i.*bn(ii,n).*No1n(3,:) + an(ii,n).*Me1n(2,:);
        
        scaling=(k(ii)./(c*k(ii)*mu)).*(k(ii)./(c*k(ii)*mu));
        
        Hs(ii,:)=scaling.*E0.*E0.*in2.*n22.*(abs(Hsr).*abs(Hsr) + abs(Hst).*abs(Hst) + abs(Hsp).*abs(Hsp));
        
    end
    
end

%%
% electric field at the magnetic mode 40 1st WGM
%
% Electric=ceil((E1(n_res,:))./max(E1(n_res,:)));
% Color=ceil((E1(n_res,:))./max(E1(n_res,:)));
%
% Electric=ceil((E1(n_res,:))./1);
% %Color=ceil((E1(n_res,:))./1);
%
% figure;
% %scatter3(x,y,z,Electric+1,Color,'filled');
% scatter3(x,y,z,Electric+1,'filled');
%
% Magnetic=ceil((H1(n_res,:))./max(H1(n_res,:)));
% Color=ceil((H1(n_res,:))./max(H1(n_res,:)));
%
% figure;
% scatter3(x,y,z,Magnetic+1,Color,'filled');

% electric field at the electric mode 40 1st WGM

% Electric=ceil(abs(E1(3,:)).*1e-3);
% %Color=ceil(abs(E1(3,:)).*(1e-4));
%
% figure;
% scatter3(x,y,z,Electric+1,'filled');
% %scatter3(x,y,z,Electric+1,Color,'filled');
%
% Magnetic=ceil(abs(H1(3,:)).*(1e11));
% %Color=ceil(abs(H1(3,:)).*(1e2));
%
% figure;
% scatter3(x,y,z,Magnetic+1,'filled');
% scatter3(x,y,z,Magnetic+1,Color,'filled');

% electric field at the magnetic mode 36 2nd WGM

% Electric=ceil(abs(E_1(5,:)).*(1e27));
% Color=ceil(abs(E_1(5,:)).*(1e27));
%
% figure;
% scatter3(x,y,z,Electric+1,Color,'filled');
%
% % magnetic field at the magnetic mode 36 2nd WGM
%
% Magnetic=ceil(abs(H_1(5,:)).*(1e27));
% Color=ceil(abs(H_1(5,:)).*(1e27));
%
% figure;
% scatter3(x,y,z,Magnetic+1,Color,'filled');

% electric field at the electric mode 36 2nd WGM

% Electric=ceil(abs(E1(5,:)).*(1e-11));
% Color=ceil(abs(E1(5,:)).*(1e-11));
%
% figure;
% scatter3(x,y,z,Electric+1,Color,'filled');
%
% Magnetic=ceil(abs(H1(5,:)).*(1e8));
% Color=ceil(abs(H1(5,:)).*(1e8));
%
% figure;
% scatter3(x,y,z,Magnetic+1,Color,'filled');

% electric field at the magnetic mode 32 3rd WGM

% Electric=ceil(abs(E_1(4,:)).*(1e27));
% Color=ceil(abs(E_1(4,:)).*(1e27));
%
% figure;
% scatter3(x,y,z,Electric+1,Color,'filled');
%
% Magnetic=ceil(abs(H_1(4,:)).*(1e27));
% Color=ceil(abs(H_1(4,:)).*(1e27));
%
% figure;
% scatter3(x,y,z,Magnetic+1,Color,'filled');

% electric field at the electric mode 32 3rd WGM

% Electric=ceil(abs(E1(3,:)).*(1e1));
% Color=ceil(abs(E1(3,:)).*(1e1));
%
% figure;
% scatter3(x,y,z,Electric+1,Color,'filled');
%
% Magnetic=ceil(abs(H1(3,:)).*(1e1));
% Color=ceil(abs(H1(3,:)).*(1e1));
%
% figure;
% scatter3(x,y,z,Magnetic+1,Color,'filled');

%% cut out a plane in Xy at a given z position

%Eresonance3=sum(Es(n_res,:))
% Hresonance4=reshape(Hs(n_res,:),Nx,Ny,Nz);
Eresonance=reshape(E1(n_res,:)+Es(n_res,:),Nx,Ny,Nz);
Hresonance=reshape(H1(n_res,:)+Hs(n_res,:),Nx,Ny,Nz);

ang=0:0.01:2*pi;
xp=kor_a*cos(ang);
yp=kor_a*sin(ang);

EEE(:,:)=abs(Eresonance(:,:,Nz/2));
FFF(:,:)=abs(Hresonance(:,:,Nz/2));

figure;
pcolor(EEE');
shading interp
hold on;
plot((Nx/2)+xp,(Ny/2)+yp,'r')
title('Electric field in XOY plane')
set(gcf,'Color',[1 1 1]);

figure;
pcolor(FFF');
shading interp
hold on;
plot((Nx/2)+xp,(Ny/2)+yp,'r')
title('Magnetic field in XOY plane')
set(gcf,'Color',[1 1 1]);

% figure
% plot(str2num(wn),real(an(:,n)),'g--',str2num(wn),real(bn(:,n)),'m--',str2num(wn),real(cn(:,n)),'r--',str2num(wn),real(dn(:,n)),'b--');
% legend('a_n','b_n','c_n','d_n');

%% cut out a plane in XZ at a given y position


Eresonance=reshape(E1(n_res,:)+Es(n_res,:),Nx,Ny,Nz);
Hresonance=reshape(H1(n_res,:)+Hs(n_res,:),Nx,Ny,Nz);

ang=0:0.01:2*pi;
xp=kor_a*cos(ang);
yp=kor_a*sin(ang);

EEE(:,:)=abs(Eresonance(:,Ny/2,:));
FFF(:,:)=abs(Hresonance(:,Ny/2,:));

figure;
pcolor(EEE');
shading interp
hold on;
plot((Nx/2)+xp,(Nz/2)+yp,'r')
title('Electric field in XOZ plane')
set(gcf,'Color',[1 1 1]);

figure;
pcolor(FFF');
shading interp
hold on;
plot((Nx/2)+xp,(Nz/2)+yp,'r');
title('Magnetic field in XOZ plane');
set(gcf,'Color',[1 1 1]);

%% cut out a plane in YZ at a given x position

EEE(:,:)=abs(Eresonance(Nx/2,:,:));
FFF(:,:)=abs(Hresonance(Nx/2,:,:));

figure;
pcolor(EEE');
shading interp
hold on;
plot((Ny/2)+xp,(Nz/2)+yp,'r')
title('Electric field in YOZ plane')
set(gcf,'Color',[1 1 1]);

figure;
pcolor(FFF');
shading interp
hold on;
plot((Ny/2)+xp,(Nz/2)+yp,'r');
title('Magnetic field in YOZ plane')
set(gcf,'Color',[1 1 1]);

