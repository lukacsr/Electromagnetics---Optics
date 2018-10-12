%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Scattering_from_a_circular_potential_08122016.m                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The program calculates the integrated abs2(Psi(x,y)) value for different
% constant refractive indexes and different absorption coefficients

%Remember exit if going to Stallo!!


clear all
close all

n_i=[0.0; 0.0003; 0.0013; 0.0053; 0.0105; 0.0158; 0.0211; 0.0263; 0.0316; 0.0368; 0.0421; 0.05];

%% ______ properties that can be changed_________
%____________________disk________________________

for ii=1:length(n_i)
    
    n = 1.4+1i*n_i(ii,1); % Refractive index
    a = 9.26e-06; % radius
    
    
    l = -20:1:20;
    
    %%______ Parameters of the incoming wave_______
    
    nu = 2000:2:2300;
    %nu=2090; %2226 2082 2220
    kArr = nu*100*2*pi;
    
    
    betaArr = kArr*n; % wavenumber inside the disk
    
    
    %% Grid
    
    % if (ShowImages)
    % Plot the wave functions
   
    aa=2; % change the frame
    
    X0=-aa*a;
    Y0=-aa*a;
    Nstep=100; % Choose a even number for section plot!
    step=2*aa*a/Nstep;
    
    for jnu = 1:length(nu);
        
        Psii = zeros(Nstep,Nstep);
        Psi=zeros(Nstep,Nstep);
        absPsii = zeros(Nstep,Nstep);
        inssumabs2Psii=0;
        outsumabs2Psii=0;
        
        k = kArr(jnu);
        beta = betaArr(jnu);
        [A_l,B_l]= constants_planewave(l,a,k,beta,n);
        
        for i=1:Nstep
            x=X0+step*double(i);
            for j=1:Nstep
                y=Y0+step*double(j);
                theta1=atan2(y,x);
                r1=sqrt(x*x+y*y);
                insPsii=0;
                outPsii=0;
                for jk = 1:length(l)
                    if (r1<a)
                        psi=besselj(l(jk),beta*r1)*exp(1i*l(jk)*theta1);
                        Psi(j,i)= B_l(jk)*psi;
                        Psii(j,i) = Psii(j,i)+ Psi(j,i);
                        insPsii=insPsii+Psii(j,i);
                    elseif (r1>a)
                        h1=(1i^l(jk))*besselj(l(jk),k*r1)*exp(1i*l(jk)*theta1);
                        h2=A_l(jk)*besselh(l(jk),1,k*r1)*exp(1i*l(jk)*theta1);
                        Psi(j,i)=h1+h2;
                        Psii(j,i) = Psii(j,i)+ Psi(j,i);
                        outPsii=outPsii+Psii(j,i);
                    end
                    
                end
                insPsii;
                outPsii;
                Psii(j,i);
                inssumabs2Psii=inssumabs2Psii+abs(insPsii).*abs(insPsii);
                outsumabs2Psii=outsumabs2Psii+abs(outPsii).*abs(outPsii);
                absPsii(j,i) = abs(Psii(j,i)).*abs(Psii(j,i));
            end
            
        end
        
        figure;
        set(gcf,'Color',[1 1 1]);
        pcolor(absPsii);
        shading interp;
        colorbar;
        string = strcat('Program III: \nu =',int2str(nu),'cm^{-1},n0=1.4, radius=9.26 micron, n_i=',num2str(imag(n)));
        title(string);
        
        
        B(jnu,:) = sum(absPsii);
        A(jnu) = sum(B(jnu));
        
        % B(jnu,:) = sum(Psii);
        % A(jnu) = sum(B(jnu,:));
        %
        % absplot=abs(A(jnu)).*abs(A(jnu));
        
        if length(nu) > 1
            
            filename = strcat('sum_of_Psii_n_',num2str(real(n)),'_n_i_',num2str(imag(n)),'nu',num2str(nu(1)),'_',num2str(nu(end)),'_9_26micron.txt');
            
            g=fopen(filename,'a');
            
            fprintf(g,'%f    ',nu(jnu));
            
            fprintf(g,'%f    ',inssumabs2Psii);
            
            fprintf(g,'%f    ',outsumabs2Psii);
            
            fprintf(g,'%f\n', A(jnu));
            
            fclose(g);
            
        else
            filename=['Psii_peak_n',num2str(real(n)),'_n_i_',num2str(imag(n)),'nu',num2str(nu),'_9_26micron','.mat'];
            
            save(filename);
        end
    end
end
%___ Plotting the figure____
%     nu = k./(2*pi*100);
%     ZAbsWaveFunction.d=abs(Psii).*abs(Psii);
%
%
%     if(1)
%         figure;
%         set(gcf,'Color',[1 1 1]);
%         pcolor(ZAbsWaveFunction.d);
%         shading interp;
%         colorbar;
%         string = strcat('Program III: \nu =',int2str(nu),'cm^{-1},n0=1.4, radius=9.26 micron, n_i=',num2str(imag(n)));
%         title(string);
%     end
%
%
%
%     if (ShowCrossSection)
%         figure; set(gcf,'Color',[1 1 1]);
%         Section=ZWaveFunction.d(uint16(Nstep/2.0),:); plot(Section);
%         title('Cross section through wave function along x-axis');
%     end

% end


%     nu = k./(2*pi*100); Xplot=abs(psi).*abs(psi);
%
%     figure; pcolor(Xplot');
% %      colormap(hsv)
%     k_string = strcat('S_', int2str(m_ind), '(k) with ', ' nu_'...
%     ,int2str(p), ' = ',  int2str(nu)); title(k_string);
%     set(gcf,'Color',[1 1 1]);
% end

exit
