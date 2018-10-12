% ***********************************************************************
%     3-D FDTD code with Mur's second order boundaries
% ***********************************************************************
%
%     Program author: Rozalia Lukacs
%
%     First version: 26.11.2014
%
%     Version: 05. 01. 2015
%
%     This MATLAB M-file implements the finite-difference time-domain
%     solution of Maxwell's curl equations over a three-dimensional
%     Cartesian space lattice comprised of uniform cubic grid cells. Each
%     unit cell within the Yee lattice contains electric field components
%     sampled along the edges of the cube and magnetic field components
%     normal to the faces of the cube.
%
%     A homogeneus sphere is modelled in air. The radius of the sphere is 10
%     micrometer.
%     The space step is choosen to be in micro meter.
%
%     At all the boundaries Mur's second order absorbing boundary
%     conditions are used.
%
%
%     The problem of plane wave source is resolved via the total
%     field/scattered field approach.
%
%     For the problem of leakage the Mathed Numerical Dispersion technic is
%     used.
%
% ***********************************************************************

clear all;

%***********************************************************************
%     Fundamental constants
%***********************************************************************

cc=2.99792458e8;            % speed of light in free space
mu0=4.0*pi*1.0e-7;          % permeability of free space
eps0=1.0/(cc*cc*mu0);       % permittivity of free space
eta0=sqrt(mu0/eps0);

nutilde=4707.6;              % wavenumber in cm
lambda=1/(100*nutilde);     % lambda must be in meter
%lambda=6.0606060606e-6;    % center frequency of source excitation
freq=cc/lambda;             % center wavelength of source excitation
omega=2.0*pi*freq;


%***********************************************************************
%     Grid parameters
%***********************************************************************

ie=120;     %number of grid cells and Ex samples along x-direction
je=120;     %number of grid cells and Ey samples along y-direction
ke=120;     %number of grid cells and Ez samples along z-direction

i0=20;      %coordinate origins for calculation of the incident field
i1=100;      %the incident field is a plane wave
j0=20;
j1=100;
k0=20;
k1=100;

kobs=60;

ic=60;      %the center of the sphere
jc=60;
kc=60;

radius=5.4e-6;           %the radius of the sphere

dx=0.5e-6;               %space increment of cubic lattice
dt=dx/(sqrt(3.0)*cc);    %time step

nmax=200;                %total number of time steps

theta=0;  %the orientation of the Kinc to the +z axis (0-180)
phi=90;    %the orientation of the Kinc to the +x axis (0-360)
psi=0;     %orientation angle for Einc to the reference direction (Kinc X z) (0-360)
%in this way we can define the polarization of the incident-wave
E0=1.0;

S=(cc*dt)/dx;            % Courant stability factor
Nlambda=lambda/dx;       % the grid sampling resolution in space cells
% per free-space wavelength

% constants for the calculation of the incident field

c3=(sind(psi)*sind(phi)+cosd(psi)*cosd(theta)*cosd(phi));
c4=(-sind(psi)*cosd(phi)+cosd(psi)*cosd(theta)*sind(phi));
c5=(-cosd(psi)*sind(theta));
c6=(cosd(psi)*sind(phi)-sind(psi)*cosd(theta)*cosd(phi));
c7=(-cosd(psi)*cosd(phi)-sind(psi)*cosd(theta)*sind(phi));
c8=(sind(psi)*sind(theta));

%c9=((i1-i0)+1)*((j1-j0)+1)*((k1-k0)+1);

%***********************************************************************
%    Numerical Phase velocities
%***********************************************************************

% calculate numerical phase velocity for phi=0

vp0=pi*cc/(Nlambda*asin((1/S)*sin((pi*S)/Nlambda)));

% calculate numerical phase velocity for arbitrary phi

A=(dx*cosd(phi))/2;
B=(dx*sind(phi))/2;
C=(1/(S^2))*(sin((pi*S)/Nlambda))^2;

ktilda=(2*pi)/lambda;

for m=1:3
    
    ktilda=ktilda-(sin(A*ktilda)*sin(A*ktilda)+sin(B*ktilda)*sin(B*ktilda)-C)/...
        (A*sin(2*A*ktilda)+B*sin(2*B*ktilda));
end

vp=omega/ktilda;

%***********************************************************************
%    Auxiliary grid for incident field
%***********************************************************5************

%approximation

%delta1=sqrt((sind(theta)^4)*((cosd(phi)^4)+(sind(phi)^4))+(cosd(theta)^4))*dx;

%rigorous calculation

%solve Eq 5.65

A=(pi*cosd(phi)*sind(theta))/Nlambda;
B=(pi*sind(phi)*sind(theta))/Nlambda;
C=(pi*cosd(theta))/Nlambda;
D=(sin((pi*S)/Nlambda))^2;

f=@(xi) S^2*(sin(A*xi)^2+sin(B*xi)^2+sin(C*xi)^2)-D;
df=@(xi) S^2*(2*A*sin(A*xi)*cos(A*xi)+2*B*sin(B*xi)*cos(B*xi)+2*C*sin(C*xi)*cos(C*xi));

[xi, erxi]=newton(f,df,1,0.1e-7,20);

A=size(xi);

%solve Eq 5.67

B=(pi*xi(A(1,2))*S)/Nlambda;
C=sin((pi*S)/Nlambda);

g=@(s) s*sin(B/s)-C;
dg=@(s) sin(B/s)-cos(B/s)*(B/s);

[s, eres]=newton(g,dg,0.5,0.5e-5,10);

A=size(s);

delta1=(cc*dt)/s(A(1,2));

%***********************************************************************
%     Material parameters
%***********************************************************************

media=2; % air and PMMA

eps=[1.0 2.22904];
sig=[0.0 0.0];
mur=[1.0 0.999992];
sim=[0.0 0.0];

%***********************************************************************
%     Updating coefficients
%***********************************************************************

for i=1:media
    eaf  =dt*sig(i)/(2.0*eps0*eps(i));
    ca(i)=(1.0-eaf)/(1.0+eaf);
    cb(i)=(dt/(eps0*eps(i)*dx))/(1+eaf);
    haf  =dt*sim(i)/(2.0*mu0*mur(i));
    da(i)=(1.0-haf)/(1.0+haf);
    db(i)=(dt/(mu0*mur(i)*dx))/(1+haf);
end

%***********************************************************************
%     Field arrays
%***********************************************************************

ex=zeros(ie+1,je+2,ke+2);   % ex in time step n+1
ey=zeros(ie+2,je+1,ke+2);   % ey in time step n+1
ez=zeros(ie+2,je+2,ke+1);   % ez in time step n+1

hx=zeros(ie+2,je+1,ke+1);   % hx in time step n+1
hy=zeros(ie+1,je+2,ke+1);   % hy in time step n+1
hz=zeros(ie+1,je+1,ke+2);   % hz in time step n+1

Exinc=zeros(ie+2,je+2,ke+2); % Exinc in time step n+1
Eyinc=zeros(ie+2,je+2,ke+2); % Eyinc in time step n+1
Ezinc=zeros(ie+2,je+2,ke+2); % Ezinc in time step n+1

Hxinc=zeros(ie+2,je+2,ke+2); % Hxinc in time step n+1
Hyinc=zeros(ie+2,je+2,ke+2); % Hyinc in time step n+1
Hzinc=zeros(ie+2,je+2,ke+2); % Hzinc in time step n+1

ex0=zeros(ie+2,je+2,ke+2);  % ex in time step n-1
ex1=zeros(ie+2,je+2,ke+2);  % ex in time step n

ey0=zeros(ie+2,je+2,ke+2);  % ey in time step n-1
ey1=zeros(ie+2,je+2,ke+2);  % ey in time step n

ez0=zeros(ie+2,je+2,ke+2);  % ez in time step n-1
ez1=zeros(ie+2,je+2,ke+2);  % ez in time step n

Einc=zeros(1,5*ie*je);
Hinc=zeros(1,5*ie*je);

Einc0=zeros(1,5*ie*je);       % Einc at time step n-1
Einc1=zeros(1,5*ie*je);       % Einc at time step n

Eincd=zeros(ie+2,je+2,ke+2);
Hincd=zeros(ie+2,je+2,ke+2);

%***********************************************************************
%     Geometry specification (main grid)
%***********************************************************************

%  Initialize entire main grid to free space (vacuum)

caex(1:ie+1,1:je+2,1:ke+2)=ca(1);
cbex(1:ie+1,1:je+2,1:ke+2)=cb(1);

caey(1:ie+2,1:je+1,1:ke+2)=ca(1);
cbey(1:ie+2,1:je+1,1:ke+2)=cb(1);

caez(1:ie+2,1:je+2,1:ke+1)=ca(1);
cbez(1:ie+2,1:je+2,1:ke+1)=cb(1);

dahx(1:ie+2,1:je+1,1:ke+1)=da(1);
dbhx(1:ie+2,1:je+1,1:ke+1)=db(1);

dahy(1:ie+1,1:je+2,1:ke+1)=da(1);
dbhy(1:ie+1,1:je+2,1:ke+1)=db(1);

dahz(1:ie+1,1:je+1,1:ke+2)=da(1);
dbhz(1:ie+1,1:je+1,1:ke+2)=db(1);

%  Add object

for j=j0:j1
    for k=k0:k1
        for i=i0:i1
            d=sqrt((i-ic)^2+(j-jc)^2+(k-kc)^2);
            if  ((d*dx)<=radius)
                caex(i,j,k)=ca(2);
                cbex(i,j,k)=cb(2);
                caey(i,j,k)=ca(2);
                cbey(i,j,k)=cb(2);
                caez(i,j,k)=ca(2);
                cbez(i,j,k)=cb(2);
            end
        end
    end
end

%***********************************************************************
%     Coordinate origins for the incident field
%***********************************************************************

if (0 <= theta)&&(theta <= 90)
    if (0<=phi)&&(phi<=90)
        O=[i0 j0 k0];
    end
    
    if (90<phi)&&(phi<=180)
        O=[i1 j0 k0];
    end
    
    if (180<phi)&&(phi<=270)
        O=[i1 j1 k0];
    end
    
    if (270<phi)&&(phi<=360)
        O=[i0 j1 k0];
    end
end

if (90 < theta)&&(theta <= 180)
    if (0<=phi)&&(phi<=90)
        O=[i0 j0 k1];
    end
    
    if (90<phi)&&(phi<=180)
        O=[i1 j0 k1];
    end
    
    if (180<phi)&&(phi<=270)
        O=[i1 j1 k1];
    end
    
    if (270<phi)&&(phi<=360)
        O=[i0 j1 k1];
    end
end

% the normal unit vector of the incident plane wave

Kinc=[sind(theta)*cosd(phi) sind(theta)*sind(phi) cosd(theta)];

figure
set(gcf,'DoubleBuffer','on')

%***********************************************************************
%     BEGIN TIME-STEPPING LOOP
%***********************************************************************
for n=1:nmax
    
    %***********************************************************************
    %     Calculation of the incident field (introduction of a plane wave)
    %***********************************************************************
    
    % Generation of Look-Up Table
    
    %Einc(1)=E0*exp(omega*n*dt-dx/(vp*dt));
    %Einc(1)=E0*exp(omega*n*dt);
    %Einc(1)=E0*sin(omega*n*dt);
    %Einc(1)=E0*sin(omega*n*dt-delta1/(vp*dt));
    %Einc(1)=E0*cos(omega*n*dt-delta1/(vp*dt));
    Einc(1)=E0*cos(omega*n*dt);
    
    % Einc(3)=origin of the coordinate system 
    % do I need to chnge this
    
    for m0=1:5*ie*je-1
        Hinc(1,m0)=Hinc(1,m0)+(dt/((vp0/vp)*mu0*delta1))*(Einc(1,m0+1)-Einc(1,m0));
    end
    
    for m0=2:5*ie*je
        Einc(1,m0)=Einc(1,m0)+(dt/((vp0/vp)*eps0*delta1))*(Hinc(1,m0)-Hinc(1,m0-1));
    end
    
    for i=i0-1:i1
        for j=j0-1:j1
            for k=k0-1:k1
                Rcomp=[i-O(1) j-O(2) k-O(3)];
                d=dot(Kinc,Rcomp);
                d=(dx/delta1)*d;
                dp=d-floor(d);
                Eincd(i,j,k)=(1-dp)*Einc(6+floor(d))+dp.*Einc(6+floor(d)+1);
                dp=d+0.5-floor(d+0.5);
                Hincd(i,j,k)=(1-dp).*Hinc(6+floor(d+0.5)-1)+dp.*Hinc(6+floor(d+0.5));
            end
        end
    end
    
    % Incident-Field Components
    
    Hxinc=c3.*Hincd;
    Hyinc=c4.*Hincd;
    Hzinc=c5.*Hincd;
    Exinc=c6.*Eincd;
    Eyinc=c7.*Eincd;
    Ezinc=c8.*Eincd;
    
    %***********************************************************************
    %     Store the electric fields Ex, Ey and Ez from previous time steps
    %***********************************************************************
    
    if n==2
        ex1=ex(1:ie+1,1:je+2,1:ke+2);
        ey1=ey(1:ie+2,1:je+1,1:ke+2);
        ez1=ez(1:ie+2,1:je+2,1:ke+1);
        
        Einc1=Einc;
    end
    
    if n>=3
        ex0=ex1;
        ey0=ey1;
        ez0=ez1;
        ex1=ex(1:ie+1,1:je+2,1:ke+2);
        ey1=ey(1:ie+2,1:je+1,1:ke+2);
        ez1=ez(1:ie+2,1:je+2,1:ke+1);
        
        Einc0=Einc1;
        Einc1=Einc;
    end
    
    %***********************************************************************
    %     Update magnetic fields
    %***********************************************************************
    % following Schneider
    
    hx(1:ie+2,1:je+1,1:ke+1)=hx(1:ie+2,1:je+1,1:ke+1)+dbhx(1:ie+2,1:je+1,1:ke+1).*...
        ((ey(1:ie+2,1:je+1,2:ke+2)-ey(1:ie+2,1:je+1,1:ke+1))-...
        (ez(1:ie+2,2:je+2,1:ke+1)-ez(1:ie+2,1:je+1,1:ke+1)));
    
    hy(1:ie+1,1:je+2,1:ke+1)=hy(1:ie+1,1:je+2,1:ke+1)+dbhy(1:ie+1,1:je+2,1:ke+1).*...
        ((ez(2:ie+2,1:je+2,1:ke+1)-ez(1:ie+1,1:je+2,1:ke+1))-...
        (ex(1:ie+1,1:je+2,2:ke+2)-ex(1:ie+1,1:je+2,1:ke+1)));
    
    hz(1:ie+1,1:je+1,1:ke+2)=hz(1:ie+1,1:je+1,1:ke+2)+dbhz(1:ie+1,1:je+1,1:ke+2).*...
        ((ex(1:ie+1,2:je+2,1:ke+2)-ex(1:ie+1,1:je+1,1:ke+2))-...
        (ey(2:ie+2,1:je+1,1:ke+2)-ey(1:ie+1,1:je+1,1:ke+2)));
    
    %***********************************************************************
    %     Add correction terms to the magnetic field at TF/SF interface
    %***********************************************************************
 
    % after schneider's program --- I changed the signs
    
    %i=i0 face
    
    hz(i0-1,j0:j1-1,k0:k1)=hz(i0-1,j0:j1-1,k0:k1)-dbhz(i0-1,j0:j1-1,k0:k1).*Eyinc(i0,j0:j1-1,k0:k1);
    hy(i0-1,j0:j1,k0:k1-1)=hy(i0-1,j0:j1,k0:k1-1)+dbhy(i0-1,j0:j1,k0:k1-1).*Ezinc(i0,j0:j1,k0:k1-1);
    
    %i=i1 face
    
    hz(i1,j0:j1-1,k0:k1)=hz(i1,j0:j1-1,k0:k1)+dbhz(i1,j0:j1-1,k0:k1).*Eyinc(i1,j0:j1-1,k0:k1);
    hy(i1,j0:j1,k0:k1-1)=hy(i1,j0:j1,k0:k1-1)-dbhy(i1,j0:j1,k0:k1-1).*Ezinc(i1,j0:j1,k0:k1-1);
    
    %j=j0 face
    
    hz(i0:i1-1,j0-1,k0:k1)=hz(i0:i1-1,j0-1,k0:k1)+dbhz(i0:i1-1,j0-1,k0:k1).*Exinc(i0:i1-1,j0,k0:k1);
    hx(i0:i1,j0-1,k0:k1-1)=hx(i0:i1,j0-1,k0:k1-1)-dbhx(i0:i1,j0-1,k0:k1-1).*Ezinc(i0:i1,j0,k0:k1-1);
    
    %j=j1 face
    
    hz(i0:i1-1,j1,k0:k1)=hz(i0:i1-1,j1,k0:k1)-dbhz(i0:i1-1,j1,k0:k1).*Exinc(i0:i1-1,j1,k0:k1);
    hx(i0:i1,j1,k0:k1-1)=hx(i0:i1,j1,k0:k1-1)+dbhx(i0:i1,j1,k0:k1-1).*Ezinc(i0:i1,j1,k0:k1-1);
    
    %k=k0 face
    
    hy(i0:i1-1,j0:j1,k0-1)=hy(i0:i1-1,j0:j1,k0-1)-dbhy(i0:i1-1,j0:j1,k0-1).*Exinc(i0:i1-1,j0:j1,k0);
    hx(i0:i1,j0:j1-1,k0-1)=hx(i0:i1,j0:j1-1,k0-1)+dbhx(i0:i1,j0:j1-1,k0-1).*Eyinc(i0:i1,j0:j1-1,k0);
    
    %k=k1 face
    
    hy(i0:i1-1,j0:j1,k1)=hy(i0:i1-1,j0:j1,k1)+dbhy(i0:i1-1,j0:j1,k1).*Exinc(i0:i1-1,j0:j1,k1);
    hx(i0:i1,j0:j1-1,k1)=hx(i0:i1,j0:j1-1,k1)-dbhx(i0:i1,j0:j1-1,k1).*Eyinc(i0:i1,j0:j1-1,k1);    
     
    %***********************************************************************
    %     Update electric fields (Ex, Ey and Ez) in main grid
    %***********************************************************************
    %
    % following Schneider
    
    ex(1:ie+1,2:je+1,2:ke+1)=caex(1:ie+1,2:je+1,2:ke+1).*ex(1:ie+1,2:je+1,2:ke+1)+...
        cbex(1:ie+1,2:je+1,2:ke+1).*...
        ((hz(1:ie+1,2:je+1,2:ke+1)-hz(1:ie+1,1:je,2:ke+1))-...
        (hy(1:ie+1,2:je+1,2:ke+1)-hy(1:ie+1,2:je+1,1:ke)));
    
    ey(2:ie+1,1:je+1,2:ke+1)=caey(2:ie+1,1:je+1,2:ke+1).*ey(2:ie+1,1:je+1,2:ke+1)+...
        cbey(2:ie+1,1:je+1,2:ke+1).*...
        ((hx(2:ie+1,1:je+1,2:ke+1)-hx(2:ie+1,1:je+1,1:ke))-...
        (hz(2:ie+1,1:je+1,2:ke+1)-hz(1:ie,1:je+1,2:ke+1)));
    
    ez(2:ie+1,2:je+1,1:ke+1)=caez(2:ie+1,2:je+1,1:ke+1).*ez(2:ie+1,2:je+1,1:ke+1)+...
        cbez(2:ie+1,2:je+1,1:ke+1).*...
        ((hy(2:ie+1,2:je+1,1:ke+1)-hy(1:ie,2:je+1,1:ke+1))-...
        (hx(2:ie+1,2:je+1,1:ke+1)-hx(2:ie+1,1:je,1:ke+1)));
    
    %***********************************************************************
    %     Add correction terms to the electric field at TF/SF interface
    %***********************************************************************

    % i=i0 face
    
    ey(i0,j0:j1-1,k0:k1)=ey(i0,j0:j1-1,k0:k1)+cbey(i0,j0:j1-1,k0:k1).*Hzinc(i0-1,j0:j1-1,k0:k1);
    ez(i0,j0:j1,k0:k1-1)=ez(i0,j0:j1,k0:k1-1)-cbez(i0,j0:j1,k0:k1-1).*Hyinc(i0-1,j0:j1,k0:k1-1);
    
    % i=i1 face
    
    ey(i1,j0:j1-1,k0:k1)=ey(i1,j0:j1-1,k0:k1)-cbey(i1,j0:j1-1,k0:k1).*Hzinc(i1,j0:j1-1,k0:k1);
    ez(i1,j0:j1,k0:k1-1)=ez(i1,j0:j1,k0:k1-1)+cbez(i1,j0:j1,k0:k1-1).*Hyinc(i1,j0:j1,k0:k1-1);
    
    % j=j0 face
    
    ex(i0:i1-1,j0,k0:k1)=ex(i0:i1-1,j0,k0:k1)-cbex(i0:i1-1,j0,k0:k1).*Hzinc(i0:i1-1,j0-1,k0:k1);
    ez(i0:i1,j0,k0:k1-1)=ez(i0:i1,j0,k0:k1-1)+cbez(i0:i1,j0,k0:k1-1).*Hxinc(i0:i1,j0-1,k0:k1-1);
    
    %j=j1 face
    
    ex(i0:i1-1,j1,k0:k1)=ex(i0:i1-1,j1,k0:k1)+cbex(i0:i1-1,j1,k0:k1).*Hzinc(i0:i1-1,j1,k0:k1);
    ez(i0:i1,j1,k0:k1-1)=ez(i0:i1,j1,k0:k1-1)-cbez(i0:i1,j1,k0:k1-1).*Hxinc(i0:i1,j1,k0:k1-1);
    
    %k=k0 face
    
    ex(i0:i1-1,j0:j1,k0)=ex(i0:i1-1,j0:j1,k0)+cbex(i0:i1-1,j0:j1,k0).*Hyinc(i0:i1-1,j0:j1,k0-1);
    ey(i0:i1,j0:j1-1,k0)=ey(i0:i1,j0:j1-1,k0)-cbey(i0:i1,j0:j1-1,k0).*Hxinc(i0:i1,j0:j1-1,k0-1);
    
    %k=k1 face
    
    ex(i0:i1-1,j0:j1,k1)=ex(i0:i1-1,j0:j1,k1)-cbex(i0:i1-1,j0:j1,k1).*Hyinc(i0:i1-1,j0:j1,k1);
    ey(i0:i1,j0:j1-1,k1)=ey(i0:i1,j0:j1-1,k1)+cbey(i0:i1,j0:j1-1,k1).*Hxinc(i0:i1,j0:j1-1,k1);
    
    %***********************************************************************
    %     Second order Mur Absorbing Boundary Conditions
    %***********************************************************************
    
    % Boundary conditions are calculated from the second time step
    
    % Boundary conditions for the electric field
    
    if n>=3
        % x=0 boundary
        
        ez(1,2:je+1,2:ke)=-ez0(2,2:je+1,2:ke)+((cc*dt-dx)/(cc*dt+dx)).*...
            (ez(2,2:je+1,2:ke)+ez0(1,2:je+1,2:ke))+((2.0*dx)/(cc*dt+dx)).*...
            (ez1(1,2:je+1,2:ke)+ez1(2,2:je+1,2:ke))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ez1(1,3:je+2,2:ke)-2.0.*ez1(1,2:je+1,2:ke)+ez1(1,1:je,2:ke)+...
            ez1(2,3:je+2,2:ke)-2.0.*ez1(2,2:je+1,2:ke)+ez1(2,1:je,2:ke)+...
            ez1(1,2:je+1,3:ke+1)-2.0.*ez1(1,2:je+1,2:ke)+ez1(1,2:je+1,1:ke-1)+...
            ez1(2,2:je+1,3:ke+1)-2.0.*ez1(2,2:je+1,2:ke)+ez1(2,2:je+1,1:ke-1));
        
        ey(1,2:je,2:ke+1)=-ey0(2,2:je,2:ke+1)+((cc*dt-dx)/(cc*dt+dx)).*...
            (ey(2,2:je,2:ke+1)+ey0(1,2:je,2:ke+1))+((2*dx)/(cc*dt+dx)).*...
            (ey1(1,2:je,2:ke+1)+ey1(2,2:je,2:ke+1))+...
            ((cc*dt)^2/(2*dx*(cc*dt+dx))).*...
            (ey1(1,2:je,3:ke+2)-2.*ey1(1,2:je,2:ke+1)+ey1(1,2:je,1:ke)+...
            ey1(2,2:je,3:ke+2)-2.*ey1(2,2:je,2:ke+1)+ey1(2,2:je,1:ke)+...
            ey1(1,3:je+1,2:ke+1)-2.*ey1(1,2:je,2:ke+1)+ey1(1,1:je-1,2:ke+1)+...
            ey1(2,3:je+1,2:ke+1)-2.*ey1(2,2:je,2:ke+1)+ey1(2,1:je-1,2:ke+1));
        
        % x=ie+2 boundary
        
        ez(ie+2,2:je+1,2:ke)=-ez0(ie+1,2:je+1,2:ke)+...
            ((cc*dt-dx)/(cc*dt+dx)).*...
            (ez(ie+1,2:je+1,2:ke)+ez0(ie+2,2:je+1,2:ke))+...
            ((2.0*dx)/(cc*dt+dx)).*...
            (ez1(ie+1,2:je+1,2:ke)+ez1(ie+2,2:je+1,2:ke))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ez1(ie+1,3:je+2,2:ke)-2.0.*ez1(ie+1,2:je+1,2:ke)+ez1(ie+1,1:je,2:ke)+...
            ez1(ie+2,3:je+2,2:ke)-2.0.*ez1(ie+2,2:je+1,2:ke)+ez1(ie+2,1:je,2:ke)+...
            ez1(ie+1,2:je+1,3:ke+1)-2.0.*ez1(ie+1,2:je+1,2:ke)+ez1(ie+1,2:je+1,1:ke-1)+...
            ez1(ie+2,2:je+1,3:ke+1)-2.0.*ez1(ie+2,2:je+1,2:ke)+ez1(ie+2,2:je+1,1:ke-1));
        
        ey(ie+2,2:je,2:ke+1)=-ey0(ie+1,2:je,2:ke+1)+...
            ((cc*dt-dx)/(cc*dt+dx)).*...
            (ey(ie+1,2:je,2:ke+1)+ey0(ie+2,2:je,2:ke+1))+...
            ((2*dx)/(cc*dt+dx)).*...
            (ey1(ie+2,2:je,2:ke+1)+ey1(ie+1,2:je,2:ke+1))+...
            ((cc*dt)^2/(2*dx*(cc*dt+dx))).*...
            (ey1(ie+1,3:je+1,2:ke+1)-2.*ey1(ie+1,2:je,2:ke+1)+ey1(ie+1,1:je-1,2:ke+1)+...
            ey1(ie+2,3:je+1,2:ke+1)-2.*ey1(ie+2,2:je,2:ke+1)+ey1(ie+2,1:je-1,2:ke+1)+...
            ey1(ie+1,2:je,3:ke+2)-2.*ey1(ie+1,2:je,2:ke+1)+ey1(ie+1,2:je,1:ke)+...
            ey1(ie+2,2:je,3:ke+2)-2.*ey1(ie+2,2:je,2:ke+1)+ey1(ie+2,2:je,1:ke));
        
        % y=0 boundary
        
        ez(2:ie+1,1,2:ke)=-ez0(2:ie+1,2,2:ke)+((cc*dt-dx)/(cc*dt+dx)).*...
            (ez(2:ie+1,2,2:ke)+ez0(2:ie+1,1,2:ke))+((2.0*dx)/(cc*dt+dx)).*...
            (ez1(2:ie+1,2,2:ke)+ez1(2:ie+1,1,2:ke))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ez1(3:ie+2,1,2:ke)-2.0.*ez1(2:ie+1,1,2:ke)+ez1(1:ie,1,2:ke)+...
            ez1(3:ie+2,2,2:ke)-2.0.*ez1(2:ie+1,2,2:ke)+ez1(1:ie,2,2:ke)+...
            ez1(2:ie+1,1,3:ke+1)-2.0.*ez1(2:ie+1,1,2:ke)+ez1(2:ie+1,1,1:ke-1)+...
            ez1(2:ie+1,2,3:ke+1)-2.0.*ez1(2:ie+1,2,2:ke)+ez1(2:ie+1,2,1:ke-1));
        
        ex(2:ie,1,2:ke+1)=-ex0(2:ie,2,2:ke+1)+((cc*dt-dx)/(cc*dt+dx)).*...
            (ex(2:ie,2,2:ke+1)+ex0(2:ie,1,2:ke+1))+((2.0*dx)/(cc*dt+dx)).*...
            (ex1(2:ie,2,2:ke+1)+ex1(2:ie,1,2:ke+1))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ex1(3:ie+1,1,2:ke+1)-2.*ex1(2:ie,1,2:ke+1)+ex1(1:ie-1,1,2:ke+1)+...
            ex1(3:ie+1,2,2:ke+1)-2.*ex1(2:ie,2,2:ke+1)+ex1(1:ie-1,2,2:ke+1)+...
            ex1(2:ie,1,3:ke+2)-2.*ex1(2:ie,1,2:ke+1)+ex1(2:ie,1,1:ke)+...
            ex1(2:ie,2,3:ke+2)-2.*ex1(2:ie,2,2:ke+1)+ex1(2:ie,2,1:ke));
        
        % y=je+2 boundary
        
        ez(2:ie+1,je+2,2:ke)=-ez0(2:ie+1,je+1,2:ke)+...
            ((cc*dt-dx)/(cc*dt+dx)).*...
            (ez(2:ie+1,je+1,2:ke)+ez0(2:ie+1,je+2,2:ke))+...
            ((2.0*dx)/(cc*dt+dx)).*...
            (ez1(2:ie+1,je+1,2:ke)+ez1(2:ie+1,je+2,2:ke))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ez1(3:ie+2,je+1,2:ke)-2.0.*ez1(2:ie+1,je+1,2:ke)+ez1(1:ie,je+1,2:ke)+...
            ez1(3:ie+2,je+2,2:ke)-2.0.*ez1(2:ie+1,je+2,2:ke)+ez1(1:ie,je+2,2:ke)+...
            ez1(2:ie+1,je+1,3:ke+1)-2.0.*ez1(2:ie+1,je+1,2:ke)+ez1(2:ie+1,je+1,1:ke-1)+...
            ez1(2:ie+1,je+2,3:ke+1)-2.0.*ez1(2:ie+1,je+2,2:ke)+ez1(2:ie+1,je+2,1:ke-1));
        
        ex(2:ie,je+2,2:ke+1)=-ex0(2:ie,je+1,2:ke+1)+...
            ((cc*dt-dx)/(cc*dt+dx)).*...
            (ex(2:ie,je+1,2:ke+1)+ex0(2:ie,je+2,2:ke+1))+...
            ((2.0*dx)/(cc*dt+dx)).*...
            (ex1(2:ie,je+1,2:ke+1)+ex1(2:ie,je+2,2:ke+1))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ex1(3:ie+1,je+1,2:ke+1)-2.*ex1(2:ie,je+1,2:ke+1)+ex1(1:ie-1,je+1,2:ke+1)+...
            ex1(3:ie+1,je+2,2:ke+1)-2.*ex1(2:ie,je+2,2:ke+1)+ex1(1:ie-1,je+2,2:ke+1)+...
            ex1(2:ie,je+1,3:ke+2)-2.*ex1(2:ie,je+1,2:ke+1)+ex1(2:ie,je+1,1:ke)+...
            ex1(2:ie,je+2,3:ke+2)-2.*ex1(2:ie,je+2,2:ke+1)+ex1(2:ie,je+2,1:ke));
        
        % z=0 boundary
        
        ex(2:ie,2:je+1,1)=-ex0(2:ie,2:je+1,2)+((cc*dt-dx)/(cc*dt+dx)).*...
            (ex(2:ie,2:je+1,2)+ex0(2:ie,2:je+1,1))+((2.0*dx)/(cc*dt+dx)).*...
            (ex1(2:ie,2:je+1,2)+ex1(2:ie,2:je+1,1))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ex1(3:ie+1,2:je+1,1)-2.0.*ex1(2:ie,2:je+1,1)+ex1(1:ie-1,2:je+1,1)+...
            ex1(3:ie+1,2:je+1,2)-2.0.*ex1(2:ie,2:je+1,2)+ex1(1:ie-1,2:je+1,2)+...
            ex1(2:ie,3:je+2,1)-2.0.*ex1(2:ie,2:je+1,1)+ex1(2:ie,1:je,1)+...
            ex1(2:ie,3:je+2,2)-2.0.*ex1(2:ie,2:je+1,2)+ex1(2:ie,1:je,2));
        
        ey(2:ie+1,2:je,1)=-ey0(2:ie+1,2:je,2)+((cc*dt-dx)/(cc*dt+dx)).*...
            (ey(2:ie+1,2:je,2)+ey0(2:ie+1,2:je,1))+((2.0*dx)/(cc*dt+dx)).*...
            (ey1(2:ie+1,2:je,2)+ey1(2:ie+1,2:je,1))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ey1(3:ie+2,2:je,1)-2.*ey1(2:ie+1,2:je,1)+ey1(1:ie,2:je,1)+...
            ey1(3:ie+2,2:je,2)-2.*ey1(2:ie+1,2:je,2)+ey1(1:ie,2:je,2)+...
            ey1(2:ie+1,3:je+1,1)-2.*ey1(2:ie+1,2:je,1)+ey1(2:ie+1,1:je-1,1)+...
            ey1(2:ie+1,3:je+1,2)-2.*ey1(2:ie+1,2:je,2)+ey1(2:ie+1,1:je-1,2));
        
        % z=ke+2 boundary
        
        ex(2:ie,2:je+1,ke+2)=-ex0(2:ie,2:je+1,ke+1)+...
            ((cc*dt-dx)/(cc*dt+dx)).*...
            (ex(2:ie,2:je+1,ke+1)+ex0(2:ie,2:je+1,ke+2))+...
            ((2.0*dx)/(cc*dt+dx)).*...
            (ex1(2:ie,2:je+1,ke+1)+ex1(2:ie,2:je+1,ke+2))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ex1(3:ie+1,2:je+1,ke+1)-2.0.*ex1(2:ie,2:je+1,ke+1)+ex1(1:ie-1,2:je+1,ke+1)+...
            ex1(3:ie+1,2:je+1,ke+2)-2.0.*ex1(2:ie,2:je+1,ke+2)+ex1(1:ie-1,2:je+1,ke+2)+...
            ex1(2:ie,3:je+2,ke+1)-2.0.*ex1(2:ie,2:je+1,ke+1)+ex1(2:ie,1:je,ke+1)+...
            ex1(2:ie,3:je+2,ke+2)-2.0.*ex1(2:ie,2:je+1,ke+2)+ex1(2:ie,1:je,ke+2));
        
        ey(2:ie+1,2:je,ke+2)=-ey0(2:ie+1,2:je,ke+1)+...
            ((cc*dt-dx)/(cc*dt+dx)).*...
            (ey(2:ie+1,2:je,ke+1)+ey0(2:ie+1,2:je,ke+2))+...
            ((2.0*dx)/(cc*dt+dx)).*...
            (ey1(2:ie+1,2:je,ke+1)+ey1(2:ie+1,2:je,ke+2))+...
            ((cc*dt)^2/(2.0*dx*(cc*dt+dx))).*...
            (ey1(3:ie+2,2:je,ke+1)-2.*ey1(2:ie+1,2:je,ke+1)+ey1(1:ie,2:je,ke+1)+...
            ey1(3:ie+2,2:je,ke+2)-2.*ey1(2:ie+1,2:je,ke+2)+ey1(1:ie,2:je,ke+2)+...
            ey1(2:ie+1,3:je+1,ke+1)-2.*ey1(2:ie+1,2:je,ke+1)+ey1(2:ie+1,1:je-1,ke+1)+...
            ey1(2:ie+1,3:je+1,ke+2)-2.*ey1(2:ie+1,2:je,ke+2)+ey1(2:ie+1,1:je-1,ke+2));
        
    end
    
    %}
    
    % Boundary conditions for the incident electric field
    %
%     if n>=3
%         
% %         % If delta1=dx
% %     
% %         % x=0 boundary
% %         
% %         Einc(1)=-Einc0(2)+((cc*dt-dx)/(cc*dt+dx)).*(Einc(2)+Einc0(1))+...
% %             ((2.0*dx)/(cc*dt+dx)).*(Einc1(1)+Einc1(2));
% %         
% %         % x=h boundary
% %         
% %         Einc(n+3)=-Einc0(n+2)+((cc*dt-dx)/(cc*dt+dx)).*(Einc(n+2)+Einc0(n+3))+...
% %             ((2.0*dx)/(cc*dt+dx)).*(Einc1(n+3)+Einc1(n+2));
%     
%         % If delta1=~dx
%         
%         % x=0 boundary
%         
%         Einc(1)=-Einc0(2)+((((vp0/vp)*eps0)*dt-delta1)/(((vp0/vp)*eps0)*dt+delta1)).*(Einc(2)+Einc0(1))+...
%             ((2.0*delta1)/(((vp0/vp)*eps0)*dt+delta1)).*(Einc1(1)+Einc1(2));
%         
%         % x=h boundary
%         
%         Einc(n+3)=-Einc0(n+2)+((((vp0/vp)*eps0)*dt-delta1)/(((vp0/vp)*eps0)*dt+delta1)).*(Einc(n+2)+Einc0(n+3))+...
%             ((2.0*delta1)/(((vp0/vp)*eps0)*dt+delta1)).*(Einc1(n+3)+Einc1(n+2));
%         
%     end
    
    
    %***********************************************************************
    %     Visualize fields
    %***********************************************************************
    
    % visualize the whole simulation box
    
    timestep=int2str(n);
    tview(:,:)=ex(kobs,:,:);
    zview(:,:)=ex(:,:,kobs);
    sview(:,:)=ex(:,kobs,:);
    yview(:,:)=ey(:,kobs,:);
    zfild(:,:)=ez(:,:,kobs);
    xmagn(:,:)=hx(kobs,:,:);
    ymagn(:,:)=hy(kobs,:,:);
    zmagn(:,:)=hz(kobs,:,:);
    tviewinc(:,:)=Ezinc(:,:,kobs);
    hviewinc(:,:)=Hyinc(kobs,:,:);
    %exTF1(:,:)=ex(kobs,j0:j1,k0:k1);
    %exALL1(:,:)=ex(kobs,:,:);

    % visualize the sphere
    
    exSPxnull(:,:)=ex(kobs,40:80,40:80);
    exSPznull(:,:)=ex(40:80,40:80,kobs);
    exSPynull(:,:)=ex(40:80,kobs,40:80);
    xmagnSP(:,:)=hx(kobs,40:80,40:80);
    ymagnSP(:,:)=hy(kobs,40:80,40:80);
    zmagnSP(:,:)=hz(kobs,40:80,40:80);
    
    %{
    imagesc(tviewinc');
    shading flat;
    caxis auto;
    %caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Exinc(i=55,j,k), time step = ',timestep]);
    xlabel('j coordinate'); ylabel('k coordinate');
%
    imagesc(tviewinc');
    shading flat;
    caxis auto;
    %caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Eyinc(i,j=55,k), time step = ',timestep]);
    xlabel('i coordinate'); ylabel('k coordinate');

    imagesc(tviewinc');
    shading flat;
    caxis auto;
    %caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ezinc(i,j,k=55), time step = ',timestep]);
    xlabel('i coordinate'); ylabel('j coordinate');
    %
    imagesc(hviewinc');
    shading flat;
    caxis auto;
    %caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Hxinc(i=55,j,k), time step = ',timestep]);
    xlabel('j coordinate'); ylabel('k coordinate');
    %}
    %{
    %subplot('position',[0.15 0.45 0.7 0.45]),
    subplot('position',[0.05 0.45 0.45 0.45])
    imagesc(tview');
    shading flat;
    caxis auto;
    caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ex(i=60,j,k), time step = ',timestep]);
    xlabel('j coordinate'); ylabel('k coordinate');
    
    %F(n)=getframe;
    
    %subplot('position',[0.15 0.10 0.7 0.25]),
    subplot('position',[0.05 0.05 0.45 0.30])
    imagesc(zview');
    shading flat;
    caxis auto;
    caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ex(i,j,k=60), time step = ',timestep]);
    xlabel('i coordinate'); ylabel('j coordinate');

    %F(n)=getframe;
    %}
    
%     %subplot('position',[0.15 0.45 0.7 0.45]),
%     imagesc(sview');
%     shading flat;
%     caxis([-1.0 1.0]);
%     caxis auto;
%     colorbar;
%     axis image; axis xy;
%     title(['Ex(i,j=60,k), time step = ',timestep]);
%     xlabel('i coordinate'); ylabel('k coordinate');
  
%{
    %subplot('position',[0.55 0.45 0.35 0.30]) % 4 subplots
    subplot('position',[0.03 0.15 0.30 0.45]) % 3 subplots
    imagesc(exSPxnull');
    shading flat;
    caxis auto;
    caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ex(i=60,j,k) around the sphere , time step = ',timestep]);
    xlabel('j coordinate'); ylabel('k coordinate');
    
    F(n)=getframe;
    
    %}
    %subplot('position',[0.15 0.10 0.7 0.25]), % 2 subplots
    %subplot('position',[0.55 0.05 0.35 0.30])  % 4 subplots
    %subplot('position',[0.37 0.15 0.30 0.45])  % 3 subplots
    imagesc(exSPznull');
    shading flat;
    caxis auto;
    caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ex(i,j,k=50) around the sphere, time step = ',timestep]);
    xlabel('i coordinate'); ylabel('j coordinate');

    %F(n)=getframe;
    
    %{
    %subplot('position',[0.15 0.45 0.7 0.45]), % 4 subplots
    subplot('position',[0.70 0.15 0.30 0.45]) % 3 subplots
    imagesc(exSPynull');
    shading flat;
    caxis([-1.0 1.0]);
    caxis auto;
    colorbar;
    axis image; axis xy;
    title(['Ex(i,j=60,k) around the sphere, time step = ',timestep]);
    xlabel('i coordinate'); ylabel('k coordinate');
    %}
    %{
    subplot('position',[0.15 0.10 0.7 0.25]),
    imagesc(yview');
    shading flat;
    caxis auto;
    caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ey(i,j=55,k), time step = ',timestep]);
    xlabel('i coordinate'); ylabel('k coordinate');

    %F1(n)=getframe;
    %
    %subplot('position',[0.15 0.10 0.7 0.25]),
    imagesc(zfild');
    shading flat;
    caxis auto;
    %caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Ez(i,j,k=55), time step = ',timestep]);
    xlabel('i coordinate'); ylabel('j coordinate');
    %
    subplot('position',[0.15 0.10 0.7 0.25]),imagesc(zmagn');
    shading flat;
    caxis auto;
    %caxis([-1.0 1.0]);
    colorbar;
    axis image; axis xy;
    title(['Hz(i=55,j,k), time step = ',timestep]);
    xlabel('j coordinate'); ylabel('k coordinate');
    %}
    
    pause(0.05)
    
    %***********************************************************************
    %     END TIME-STEPPING LOOP
    %***********************************************************************
    
end