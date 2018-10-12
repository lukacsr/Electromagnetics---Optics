%% This code calculates the Bloch Green function in 2D for any given single cell
% Green function. It uses the Bloch boundary conditions.
% The grid is defined in cell 0.
%

clear all
%profile on
format long

% Choose your Green function
% 1 -- 2D free space boundary conditions Green function
% 2 -- 2D Bucket model Green function
% 3 -- 2D mirror boundary conditions on the bottom

Green_function=1;

% Define the number of cells, this must be all the time odd

NrCell=1;

% Define the interval for the wavenumber and the incoming plane wave's angle

% lambda_array = [2280 3000];
lambda_array = 1/(1282*100);
%lambda_array = 9.6432015439*10^-6;
phi= pi;
k=(2*pi)/lambda_array;
kx=k*cos(phi);
ky=k*sin(phi);

% Define the refractive index

n_index_array= [1.9]; % here we set in for ex. c-Si

% Define the resolution

% NX=100;
% NY=NX;
ArrayOfCircles=1;
Zebra=0;

% Define the total energy of the incoming wave

cc=3e8;
hh=6.626070040e-34;

K=sqrt(kx^2+ky^2);

% Define the Bloch momentum

p0=kx;
q0=ky;

%for ll = 1:length(lambda_array);

% !! Husk exit på Stallo!!

for lambda= lambda_array(1):1:lambda_array(end)
    
    switch ArrayOfCircles
        case 1
            Nx=100;             % Resolution in x-direction
            Ny=100;             % Resolution in y-direction
            N=Nx*Ny;
            NCx=1;              % No of disks in x-direction
            NCy=1;              % No of disks in y-direction
            %R=power(10,-6);     % Radius
            R = 10e-6;
            F = 0.0*R;            % Size of frame as fraction of R
            %a=2*NCx*R*(1+F);
            %b=2*NCy*R*(1+F);
            a=2*NCx*R;
            b=2*NCy*R+2*F;
            n_index=1.9;        % Refractiv index of the disk(s)
            AA=[a 0];
    end
    
    switch Zebra
        case 1
            Nx=350;             % Resolution in x-direction
            Ny=40;             % Resolution in y-direction
            N=Nx*Ny;
            T = 500 * 10^-9; % Thickness of film
            a = 12*T;         % Width of film (and test region)
            b = 1.5*T;         % Higth test region
            n_index = 3+0.004*1i;
            AA=[a 0];
    end
    
    % n_index=n_index_array(ll);
    
    dx=a/Nx; % spacing in x-direction
    dy=b/Ny; % spacing in y-direction
    
    x0=0;
    y0=0;
    %x0=F;
    %y0=F;
    Nxadd = (2*x0)./dx;
    Nyadd = (2*y0)./dy;
    Nx = Nx + Nxadd;
    Ny = Ny + Nyadd;
    N = Nx*Ny;
    
    x=zeros(1,N);
    y=zeros(1,N);
    
    %     for i=1:Nx
    %         for l=1:Ny
    %             j=(l-1)*Nx+i;
    %             x(j)=(i-0.5)*dx;
    %             y(j)=(l-0.5)*dy;
    %         end
    %     end
    
    indexi=reshape(((1:1:Nx)'*ones(1,Ny))',1,N);
    indexl=reshape(((1:1:Ny)'*ones(1,Nx)),1,N);
    
    indexj=(indexl-1).*Nx+indexi;
    x(indexj)=((indexi-0.5).*dx); % the Descartes coordinate x
    y(indexj)=((indexl-0.5).*dy); % the Descartes coordinate y
    
    % Establish the potential
    
    v=zeros(1,N,'double');
    Int_min = 0;
    Int_max = 0;
    
    for m=1:N
        
        switch ArrayOfCircles case 1
            for nx=1:NCx
                for ny=1:NCy
                    % Array of circles
                    xcircle=(2*nx-1)*R+x0;
                    ycircle=(2*ny-1)*R+F;
                    rn=sqrt((x(m)-xcircle)*(x(m)-xcircle)+(y(m)-ycircle)*(y(m)-ycircle));
                    if (rn<R)
                        v(m)=1.0-n_index*n_index;
                    end
                end
            end
        end
        
        switch Zebra case 1
            if y(m) > 0.25*T && y(m) < 1.25*T
                v(m)=1.0-n_index*n_index;
                
                if Int_min == 0
                    Int_min = m;
                end
                Int_max = m;
            end
        end
    end
    
    % Define the Bloch momentum p_n and q_n
    
    en1=floor((-a/(2*pi))*(kx+K));
    en2=floor((-a/(2*pi))*(kx-K));
    
    if en1 <= en2
        en=en1:1:en2;
    else
        en=en2:1:en1;
    end
    
    p=p0+en*2*pi/a;
    q=sqrt(K*K-p.*p);
    
    %     en1=floor((-a/(2*pi))*(ky+K));
    %     en2=floor((-a/(2*pi))*(ky-K));
    %
    %     if en1 <= en2
    %         en=en1:1:en2;
    %     else
    %         en=en2:1:en1;
    %     end
    %
    %     q=q0+en*2*pi/a;
    %     p=sqrt(K*K-q.*q);
    
    % calculate the Bloch Green function
    
    Gtilde=zeros(N,N);
    
    jj = reshape(((1:1:N)'*ones(1,N))',1,N*N);
    mm = reshape(((1:1:N)'*ones(1,N)),1,N*N);
    
    for MM=-((NrCell-1)/2):1:((NrCell-1)/2);
        
        Z = k*sqrt((x(jj)-((MM*a)+x(mm))).*(x(jj)-((MM*a)+x(mm)))+(y(jj)-y(mm)).*(y(jj)-y(mm))); % k*abs(r-r')
        
        switch Green_function
            
            case 1 % free space boundary conditions
                
                G = reshape(besselh(0,1,Z),N,N);
                G(isnan(G)) = ((1i/4)-(1/(2*pi))*(log(k*dx/2)+eulergamma)-((1/4)*log(1/2)-3/4+pi/8)/pi)*(dx*dy);
                
            case 2
                
            case 3
                
        end
        
        for ll=1:1:size(p0)
            
            Gp=exp(1i*MM*(p0(ll)*AA(1))).*G;
            Gtilde(1:N,:)=Gtilde(1:N,:)+Gp(1:N,:);
            
        end
        
        clear G Z Gp;
    end
    
    Gtilde(1:N,:)=Gtilde(1:N,:).*reshape(reshape((v'*ones(1,N))',1,N*N),N,N);
    
    %clear v;
    
    % Define the source
    
    eikr=exp(1i*(kx.*x+ky.*y));
    
    %clear x y;
    
    IdentityMat=eye(N,N);
    M=IdentityMat+1i*(k*k*dx*dy*0.25)*Gtilde;
    
    clear Gtilde;
    clear IdentityMat;
    
    %     vMat=reshape(v, [Nx,Ny]);
    %     figure;
    %     pcolor(real(vMat'))
    
    %% Solve Linear Equation by Matlab routine
    
    psi_complex = linsolve(M,eikr');
    %psi_complex_reshaped=reshape(psi_complex,Nx,Ny)';
    %psi_complex_reshaped=fliplr(psi_complex_reshaped);
    %psi_complex_reshaped=flipud(psi_complex_reshaped);
    
    %Xplot=abs(psi_complex_reshaped).*abs(psi_complex_reshaped);
    Xplot=abs(psi_complex).*abs(psi_complex);
    
    IntegralPsi=sum(Xplot(1:N));
    
    %     figure;
    %     pcolor(Xplot');
    %     %imagesc(Xplot);
    %     set(gcf,'Color',[1 1 1]);
    %     string = strcat('Bloch Green function 2D: \lambda =',num2str(lambda),' m;',' Number of cells =',num2str(NrCell));
    %     title(string);
    %     shading interp;
    %     %caxis([0 5]);
    
    phi=angle(psi_complex_reshaped); % gives angles in the range -pi:pi
    
    % Calculate the scattering amplitudes
    
    %     x=x-R;
    %     y=y-R;
    %
    %     for mmm=1:length(en);
    %
    %         Aen(mmm)=1i*(k^2)*exp(1i*(en(mmm)*pi/2))*sum(v'.*psi_complex*dx*dy);
    %         Hen=sqrt(2./(pi.*Z)).*exp(1i.*(Z-(en(mmm)*pi*0.5)-0.25*pi));
    %
    %     end
    %     for ll=-20:1:20
    %
    %     phi_prim=atan2(y,x);
    %     bess=besselj(ll,(kx*x+ky*y));
    %     integral=sum(bess.*exp(-1i*ll*phi_prim).*v.*psi_complex'*dx*dy);
    %     A_l(ll+21)=(-(k*k)/4)*sqrt(2/(pi*k))*exp(-1i*((pi/2)*ll-pi/4)*integral);
    % end
    
    % write the calculated values in files
    
    if length(lambda_array) == 1
        
        filename = strcat('Psi_2D_Bloch_n1.9_NR_cells',num2str(NrCell),'_',num2str(Nx),'x',num2str(Ny),'_',num2str(lambda_array(1)),'_10micron.mat');
        
        save(filename,'Xplot','lambda_array');
        
    else
        
        g=fopen(['2D_Bloch_NR_cells',num2str(NrCell),'_',num2str(lambda_array(1)),'_',num2str(NX),'x',num2str(NY),'_','_n1.9_10micron.txt'],'a');
        fprintf(g,'%f    ',lambda);
        fprintf(g,'%f\n',IntegralPsi);
        
        fclose(g);
        
    end
end

%exit








