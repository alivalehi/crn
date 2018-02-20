%clc; 
close all; clear all; warning off;
disp('This program checks some approximations used in equations [Not necessary for main results]');
run_mode = input('Press 1 for large LT, 2 for small LT approx, 3 for plot: ');
if isempty(run_mode), run_mode=0; end
if run_mode == 1 || run_mode==2   %run approximations
    
    
    
    epsilon=1e-4;
    %my params:    
    R = 1; lambda=1; N=50; H=30; b=0.001; a=1-b; 

    %parameter set for present Large and small LT approximation 
    H=16; N=8; lambda=10; Tmin=0.5; Exp_PER=0.2; BER=0.001; R=2000; b=BER; a=1-b;
    Tv = (1/lambda)*[0.01 0.02 0.03 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100]; %good
    LT = lambda*Tv; T=Tv;
    
    eta=H/N; N2=N^2; H2=H^2; NH=N*H; R2=R^2;    P0=exp(-LT); Pn0=1-exp(-LT);


    Etau = T./(1-P0); Etau2 = T.^2 .*((1+P0)./(1-P0).^2); sig2Tau= T.^2 .* ((P0)./(1-P0).^2);
    K2Tau = P0; 
    
    %ES has three terms
    s1 = (exp(-LT)./(R*a^H)) .* (N*LT*a^(-N)).*exp(LT.*a^(-N));
    s2 = (exp(-LT)./(R*a^H)) .* (H).*exp(LT*a^(-N));
    s3 = (exp(-LT)./(R*a^H)) .* (-H); %correction term for no dummy

    Es_dummy = s1+s2;%slotted mode
    Es = (s1+s2+s3)./Pn0; % general mode (no dummy for zero packet)

    %Two terms
    Es_1 = s1./Pn0; %Dominant term
    Es_2 = (s2+s3)./Pn0; %small term
    Es_copy = Es_1 + Es_2;  %check Es_nodummy_c=Es_nodummy

    %approximations for large LT
    if run_mode == 1 % run for LT_large
        Es_approx_large_LT1 = s1./Pn0; %(=ES1)approx for large LT, good for small LT as well
        Es_approx_large_LT2 = s1*(a^H)./Pn0; %almost the same as Es_approx_large_LT1, use a^H=1 for a~=1
        Es_approx_large_LT3 = N*LT.*exp(-LT.* (1-a^(-N))) ./(R*a^(H+N));  %paper ver1
        Es_approx_large_LT4 = N*LT.*exp(LT.* N*b)         ./(R*a^(H+N));  %paper ver2 
        fprintf('Press a Key to see various approximations for E[s] for large Lambda T \n   Lambda T = %s   ...', num2str(LT)); pause;
        Esall = [Es; Es_approx_large_LT1; Es_approx_large_LT2; Es_approx_large_LT3; Es_approx_large_LT4]  

    else
        %approx : small LT , a==>1  (a^H=1-Hb,  exp(LT a^(-kN))=1+ LT(1+akN)
        Es_approx_small_LT1 = ((H+N+b*N*(N+H))/(R*a^H))*exp(-LT);  %good approx
        Es_approx_small_LT2 = ((H+N)/(R*a^H))*exp(-LT);            %approx 2
        fprintf('Press a Key to see various approximations for E[s] for small Lambda T \n   Lambda T = %s   ...', num2str(LT)); pause;
        Esall = [Es; Es_approx_small_LT1; Es_approx_small_LT2]
    end
    %%  Compact notation: sort based on exp(LT*a^(-2*N))
    A1 = (2*N2*LT*a^(-2*N)+2*N2*LT.^2 *a^(-4*N)+ 4*NH*LT*a^(-2*N) + 2*H2).*(a^(-H) * exp(LT*a^(-2*N)));
    A2 = -(N2*LT*a^(-N) + N2*LT.^2*a^(-2*N)+2*NH*LT*a^(-N)+H2).* exp(LT*a^(-1*N));
    AH2_nodummy= H2 *(1-2*a^(-H));   

    Es2 = (A1+ A2 + AH2_nodummy) .* exp(-LT)./(R2 * a^H*(1-exp(-LT))); 
    Es2_dummy = (A1+ A2) .* exp(-LT)./(R2 * a^H); 

    sigma2s = Es2-Es.^2;
    sigma2s_dummy = Es2_dummy-Es_dummy.^2; 

    Ks=sqrt(sigma2s)./Es;  K2s = sigma2s./(Es.^2); %Coef of variation
    Ks_dummy=sqrt(sigma2s_dummy)./Es_dummy;


    %Kingman Equation for GI/GI/1
    rho = Es ./ Etau;
    Ef = (T/2);  % ~=T/2
    Ew = 0.5 * (rho./ (1-rho)) .* Es .* (K2s + K2Tau);  %Kingman Formula
    Ed = Ef + Es + Ew;  %Original quantification

    if run_mode == 1 % run for LT_large
        %Approx for large LT and alpha<<1 ==> a^(-2N)>> a^(-N) >> 1
        Es2_approx_large_LT1 = A1 .* exp(-LT)./(R2 * a^H*(1-exp(-LT))); %not good
        Es2_approx_large_LT2 = (2*N2*LT.^2 *a^(-4*N)).* (a^(-H)*exp(LT*a^(-2*N)))  .* exp(-LT)./(R2 * a^H*(1-exp(-LT)));  %choose one term of A1
        Es2_approx_large_LT3 = (2*N2*LT.^2 *a^(-4*N)).* (a^(-H)*exp(LT*a^(-2*N)))  .* exp(-LT)./(R2 * a^H); %Approx 2: except (1-exp(-LT)), paper version
        Es2_approx_large_LT4 = (2*N2*LT.^2)          .*         exp(2*LT*b*N)                 ./(R2 * a^(2*H+4*N)); %paper ver 2


        fprintf('Press a Key to see various approximations for E[s2] for large Lambda T \n   Lambda T = %s   ...', num2str(LT)); pause;
        Es2all=[Es2; Es2_approx_large_LT1; Es2_approx_large_LT2; Es2_approx_large_LT3; Es2_approx_large_LT4]

        %Kingman Equation

        %approx (1): use approximations
        Etau_1=Etau;  K2Tau_1 = K2Tau; Ef_1=Ef; 
        Es_1 = Es_approx_large_LT1; rho_1 = Es_1 ./ Etau_1;   K2s_1 = Es2_approx_large_LT1./ Es_approx_large_LT1.^2-1; 
        Ew_1 = 0.5 * (rho_1./ (1-rho_1)) .* Es_1 .* (K2s_1 + K2Tau_1);
        Ed_1 = Ef_1 + Es_1 + Ew_1;

        %approx (2): use further approximations
        Etau_2=Etau;  K2Tau_2 = K2Tau; Ef_2=Ef; 
        Es_2 = Es_approx_large_LT2; rho_2 = Es_2 ./ Etau_2;   K2s_2 = Es2_approx_large_LT2./Es_approx_large_LT2.^2 -1; 
        Ew_2 = 0.5 * (rho_2./ (1-rho_2)) .* Es_2 .* (K2s_2 + K2Tau_2);
        Ed_2 = Ef_2 + Es_2 + Ew_2;

        %approx (3): use further approximations
        Etau_3=Etau;  K2Tau_3 = K2Tau; Ef_3=Ef; 
        Es_3 = Es_approx_large_LT3; rho_3 = Es_3 ./ Etau_3;   K2s_3 = Es2_approx_large_LT3./Es_approx_large_LT3.^2 -1; 
        Ew_3 = 0.5 * (rho_3./ (1-rho_3)) .* Es_3 .* (K2s_3 + K2Tau_3);
        Ed_3 = Ef_3 + Es_3 + Ew_3;

        %approx (4): use further approximations
        Etau_4=Etau;  K2Tau_4 = K2Tau; Ef_4=Ef; 
        Es_4 = Es_approx_large_LT4; rho_4 = Es_4 ./ Etau_4;   K2s_4 = Es2_approx_large_LT4./Es_approx_large_LT4.^2 -1; 
        Ew_4 = 0.5 * (rho_4./ (1-rho_4)) .* Es_4 .* (K2s_4 + K2Tau_4);
        Ed_4 = Ef_4 + Es_4 + Ew_4;

        % approx (5): use more simplification and close form format
        K2s_5 = 1+4*N*b;
        Etau_5= T;  K2Tau_5 = 0; Ef_5=T/2; 
        Es_5 = Es_approx_large_LT4; rho_5 = Es_5 ./ Etau_5;   
        Ew_5 = 0.5 * (rho_5./ (1-rho_5)) .* Es_5 .* (K2s_5 + K2Tau_5);
        Ed_5 = Ef_5 + Es_5 + Ew_5;


        fprintf('Press a Key to see various approximations for E[s] for large Lambda T \n   Lambda T = %s   ...', num2str(LT)); pause;
        Edall = [Ed; Ed_1; Ed_2; Ed_3; Ed_4; Ed_5]
        save('large_LT');
    else

%%     %Approx for small LT, small b=>0

        %Original form
        Es2 = (A1+ A2 + AH2_nodummy) .* exp(-LT)./(R2 * a^H*(1-exp(-LT))); 

        %remove LT.^2 terms
        A1c=(2*N2*LT*a^(-2*N)+4*NH*LT*a^(-2*N)+2*H2) *a^(-H) .* exp(LT*a^(-2*N));
        A2c=-(N2*LT*a^(-N)+2*NH*LT*a^(-N)+H2) .* exp(LT*a^(-N));
        A3c=H2*(1-2*a^(-H)); %=AH2_nodummy
        Es2_approx_small_LT0 = (exp(-LT)./((1-exp(-LT))*R2*a^H)).* (A1c + A2c  + A3c);



        %apply approx for some a^N=1-aN, exp(LT a^(-N))=1+LT a^(-N) = 1+LT(1+b*N)=1+LT+LT*b*N
        A1c2=(2*N*(N+2*H)*LT* (1+2*b*N)+2*H2) *a^(-H) .* exp(LT*a^(-2*N));
        A2c2=-(N*(N+2*H) *LT* (1 + b*N)  +H2) .* exp(LT*a^(-N));
        A3c2=H2*(1-2*a^(-H)); %=AH2_nodummy
        Es2_approx_small_LT0 = (exp(-LT)./((1-exp(-LT))*R2*a^H)).* (A1c2 + A2c2  + A3c2);


        %apply approx for all a^N=1-aN, exp(LT a^(-N))=1+LT a^(-N) = 1+LT(1+b*N)=1+LT+LT*b*N
        A1c2=(2*N*(N+2*H)*LT* (1+2*b*N)+2*H2) *(1+b*H) .* (1 + LT*(1 + 2*b*N));
        A2c2=-(N*(N+2*H) *LT* (1 + b*N)  +H2) .*  (1 + LT*(1 + b*N));
        A3c2=H2*(1-2*(1+b*H)); %=AH2_nodummy
        Es2_approx_small_LT0 = (exp(-LT)./((1-exp(-LT))*R2*a^H)).* (A1c2 + A2c2  + A3c2);


        %Not good
        Es2_approx_small_LT1 = (exp(-LT)./((1-exp(-LT))*R2*a^H)).* LT * ((H+N)^2 + b * (2*H^3 + 7*H^2*N + 3*N^3 + 8*N^2*H));
        K2s_small_LT1 = Es2_approx_small_LT1 ./ (Es_approx_small_LT1).^2 -1;

        Es2_approx_small_LT2 = (exp(-LT)./((1-exp(-LT))*R2*a^H)).* LT * (H+N)^2;
        K2s_small_LT2 = Es2_approx_small_LT2 ./ (Es_approx_small_LT2).^2 - 1;
        Es2all = [Es2; Es2_approx_small_LT0; Es2_approx_small_LT1; Es2_approx_small_LT2], disp('Es2 small LT'); pause;


        %Sigma2S approx
        K2sall = [K2s; K2s_small_LT1; K2s_small_LT2], disp('sigma for small LT'); pause;

        %Kingman Equation
        %USE Kingman with Original Parameters
        rho = Es ./ Etau;
        Ef = (T/2);
        Ew = 0.5 * (rho./ (1-rho)) .* Es .* (K2s + K2Tau);  %Kingman Formula
        Ed = Ef + Es + Ew;  %Original quantification

        %approx (2) use approximations
        Etau_1=Etau;  K2Tau_1 = K2Tau; Ef_1=Ef; 
        Es_1 = Es_approx_small_LT1; rho_1 = Es_1 ./ Etau_1;   K2s_1 = Es2_approx_small_LT1 ./ Es_approx_small_LT1.^2 -1; 
        Ew_1 = 0.5 * (rho_1./ (1-rho_1)) .* Es_1 .* (K2s_1 + K2Tau_1);
        Ed_1 = Ef_1 + Es_1 + Ew_1;

        % approx (3) more approximations, (Not good)
        Etau_2=Etau; K2Tau_2 = 0; Ef_2=Ef;
        Es_2 = Es_approx_small_LT2; rho_2 = Es_2 ./ Etau_2;   K2s_2 = Es2_approx_small_LT2 ./ Es_approx_small_LT2.^2 -1; 
        Ew_2 = 0.5 * (rho_2./ (1-rho_2)) .* Es_2 .* (K2s_2 + K2Tau_2);
        Ed_2 = Ef_2 + Es_2 + Ew_2;


        % approx (4): use more simplification and close form format
        Etau_3=1/lambda; K2Tau_3 = 1/lambda^2; Ef_3=T;
        Es_3 = Es_approx_small_LT2; rho_3 = Es_3 ./ Etau_3;   K2s_3 = Es2_approx_small_LT2 ./ Es_approx_small_LT2.^2 - 1; 
        Ew_3 = 0.5 * (rho_3./ (1-rho_3)) .* Es_3 .* (K2s_3 + K2Tau_3);
        Ed_3 = Ef_3 + Es_3 + Ew_3;

        %closed form notation of Ed_3
        Ed_4 = ((N+H)*exp(-LT)/(R*a^H)).*((N+H)*exp(-LT)./(2*(R*a^H-lambda*(N+H)*exp(-LT))) + 1) + T;

        fprintf('Press a Key to see various approximations for E[d] for small Lambda T \n   Lambda T = %s   ...', num2str(LT)); pause;
        Edall = [Ed; Ed_1; Ed_2; Ed_3; Ed_4]
        save('small_LT');

    end
    
elseif run_mode == 3% plot approximations
    
    plot_Ed =1;
    if plot_Ed
        if isempty(dir('small_LT.mat')) || isempty(dir('large_LT.mat')), error('Please first the program for mode=1 and mode=2 to calculate approximations'); end 
        load ('small_LT', 'Ed_4', 'Ed', 'LT');
        Ed_small=Ed_4;    Ed1=Ed;
        load ('large_LT', 'Ed_5');
        Ed_large=Ed_5;

        f=figure; plotdef(0);
        loglog(LT,Ed, 'b.-', 'linewidth', 2); hold on;
        loglog(LT,Ed_small, 'g--o','linewidth', 2);  
        loglog(LT,Ed_large, 'r--S','linewidth', 2); 
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex'); ylabel('Expected Delay: E[d]', 'interpreter', 'latex'); l=legend('Exact E[d]', 'Approx: Small $\lambda T$', 'Approx: Large $\lambda T$'); set(l, 'interpreter', 'latex'); 
        title(''); pause; saveas(gcf, 'Approx_Ed','png');  
    end    
    
    
    plot_ES=1;
    if plot_ES
        
        %%%approximates for large LT
        load ('large_LT');
        figure;  plotdef(0); [~,ind]=find(LT>1);
        semilogy(LT(ind),Es_dummy(ind),'b->');   hold on; 
        semilogy(LT(ind),Es(ind),'g-S');  
        semilogy(LT(ind),Es_approx_large_LT1(ind),'r:*');  
        semilogy(LT(ind),Es_approx_large_LT2(ind),'g:S');  
        semilogy(LT(ind),Es_approx_large_LT3(ind),'m:<');  
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex'); ylabel ('E[s]: Approximates for Large $\lambda T$',  'interpreter', 'latex');
        title(''); pause; saveas(gcf, 'Approx_Es_LargeT','png');  

        %%%approximates for small LT
        load ('small_LT');
        figure;   [~,ind]=find(LT<1); plotdef(0);
        loglog(LT(ind),Es(ind),'g-S');   hold on;   
        loglog(LT(ind),Es_approx_small_LT1(ind),'r:*');  
        loglog(LT(ind),Es_approx_small_LT2(ind),'m:S');  
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex'); ylabel('E[s]: Approximates for Small $\lambda T$', 'interpreter', 'latex'); 
        title(''); pause; saveas(gcf, 'Approx_Es_smallT','png');  
    end




    plot_Es2=1; 
    if plot_Es2 

        load ('large_LT');
        figure; plotdef(0);
        semilogy(LT,Es.^2,'b->'); hold on; 
        semilogy(LT,Es_dummy.^2,'g-S');  
        semilogy(LT,abs(Es2),'c-*');     
        semilogy(LT,abs(Es2_dummy),'r-o'); l = legend('$E[s]^2$', '$E[\tilde{s}]^2$', '$E[s^2]$', '$E[\bar{s}]$'); set(l, 'interpreter', 'latex');
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex');
        ylabel('$(E[s])^2$ and $E[s^2]$', 'interpreter', 'latex'); hold on; %ylabel('$(E[s])^2$ and $E[s^2]$: Comp for {\bf s} and {\bf \tilde{s}}'); hold on; 
        title(''); pause; saveas(gcf, 'Approx_Es2_compare','png');  
 
        figure; plotdef(0);
        semilogy(LT,Es2,'b->');  l=title('E[s^2]: Approx for large \lambda T'); set(l,'Interpreter','Latex'); hold on; 
        semilogy(LT,Es2_approx_large_LT1,'g-S'); 
        semilogy(LT,Es2_approx_large_LT2,'c-o'); 
        semilogy(LT,Es2_approx_large_LT3,'m-x'); 
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex'); ylabel('E[s^2]');
        title(''); pause; saveas(gcf, 'Approx_Es2_largeT','png');   
        
        figure; plotdef(0);
        semilogy(LT,sigma2s,'b->');  hold on; %title('\sigma^2[s]: Comp for {\bf s} and {\bf \tilde{s}}');
        semilogy(LT,sigma2s_dummy,'g-S');  l = legend('$\sigma^2[s]$', '$\sigma^2[{\tilde{s}}]$'); set(l,'Interpreter','Latex');
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex');
        ylabel('$\sigma^2[s]$','Interpreter','Latex');
        pause; title(''); saveas(gcf, 'Approx_Vars_compare','png');  



        figure; plotdef(0);
        semilogy(LT,Ks,'b->'); hold on; %title('K[s] = \sigma[s]/ E[s]: Comp for {\bf s} and {\bf \tilde{s}}');  
        semilogy(LT,Ks_dummy,'g-S');   
        l=legend('$K[s]$', '$K[\tilde{s}]$');  set(l,'Interpreter','Latex'); 
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex');
        ylabel('$K[s] = \sigma[s]/ E[s]$','Interpreter','Latex', 'interpreter', 'latex');
        pause; title(''); saveas(gcf, 'Approx_Ks_compare','png');  

        figure; plotdef(0);
        semilogy(LT,sigma2s,'b->');  hold on; %l=title('$\sigma^2[s]$: Approx for large $\lambda T$'); set(l,'Interpreter','Latex'); hold on; 
        semilogy(LT,Es2_approx_large_LT3-Es_approx_large_LT2.^2,'g-S'); 
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex');
        ylabel('$\sigma^2[s]$: Approx for large $\lambda T$', 'interpreter', 'latex');
        pause; title(''); saveas(gcf, 'Approx_Vars_largeT','png');  

        load ('small_LT'); figure; plotdef(0);
        loglog(LT,Es2,'b->'); hold on; % title('E[s^2]: Approx for small \lambda T'); 
        loglog(LT,Es2_approx_small_LT2,'r-.');  
        xlabel('Packetization Interval: $\lambda T$', 'interpreter', 'latex');
        ylabel('E[s^2]: Approx for small \lambda T', 'interpreter', 'latex');
        ylim([1e-6 1e-2]); pause; title(''); saveas(gcf, 'Approx_Es2_smallT','png');  
    end
end
