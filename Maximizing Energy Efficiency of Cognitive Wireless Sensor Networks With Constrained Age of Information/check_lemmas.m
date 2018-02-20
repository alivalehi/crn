%This program plots Energy Efficiency curves
%Also performs some additional check used in lemmas
%Previous name was: check_lemmas

close all; clc; clear all; 
disp('Some assuring checks for lemmas');


plot_poiss_pdf=0; %Extra Check for Poisson random value generation process [later remove]
chk_lemmas = 0;  %Some Extra Checks in Lemmas [Later remove]
check_Ef=0;    %Some Extra Check for Energy Efficiency [Later remove]
plot_expint=1; %Check exponential integral

clrs = {'b','g','r','m','c','k','y',...
        'b-.','g-.','r-.','m-.','c-.','k-.','y-.', ...
        'b--','g--','r--','m--','c--','k--','y--',...
        'b-*','g-*','r-*','m-*','c-*','k-*','y-*', ...
        'b-x','g-x','r-x','m-x','c-x','k-x','y-x'};
mu = 4; Flen = 1e6; kmax = max(10,5*mu);
N = 8; H=40;  %bits and header
el = 0.5772156649; % Euler–Mascheroni constant
kv = poissrnd(mu, [1, Flen]); kkv=kv(kv~=0); %remove zeros
n=1; beta=0.1; alpha=1-beta;  eta= alpha^n; 

if plot_poiss_pdf
    z=zeros(1,kmax); fk1=z; fk2=z; kval=[0:kmax-1];
    for k= kval  %empirical pdf
        fk1(k+1)=sum(kv==k)/length(kv);
        fk2(k+1) = exp(-mu) * (mu^k)/factorial(k);
    end
    fk3=poisspdf(kval, mu);
    figure; subplot(121); title('pdf (Poisson)');
    stem(kval, fk1,'bo', 'markersize', 8); hold on; plot(kval, fk2,'r*'); plot(kval, fk2,'cS'); legend('pdf: empirical', 'pdf: theoretical');



    %remove zeros
    z=zeros(1,kmax); fkk1=z; fkk2=z; kval=[0:kmax-1];
    for k= kval  %empirical pdf
        fkk1(k+1)=sum(kkv==k)/length(kkv);
        if k ~= 0, fkk2(k+1) = exp(-mu) * (mu^k)/factorial(k) ./ (1-exp(-mu)); end
    end

    subplot(122); title('pdf (Poisson  like)');
    stem(kval, fkk1,'bo', 'markersize', 8); hold on; plot(kval, fkk2,'r*'); legend('pdf: empirical', 'pdf: theoretical');
    pause;
end

%%
if chk_lemmas
    %check lemma 1: E[k^n/eta^k]  for k:poisson;
    fprintf('\n\n*******************************************\n    Check Lemma 3.1 \n*******************************************\n');
    stirling2table = [1  nan(1,10);      %horizontal index k:0,1,.. vertical index n:0,1,...
     0	1   nan(1,9);
     0	1	1   nan(1,8);
     0	1	3	1       nan(1,7);
     0	1	7	6       1       nan(1,6);
     0	1	15	25      10      1       nan(1,5);
     0	1	31	90      65      15      1       nan(1,4);
     0	1	63	301     350     140 	21      1       nan(1,3);
     0	1	127	966     1701	1050	266     28      1       nan(1,2);
     0	1	255	3025	7770	6951	2646	462     36      1   nan;
     0	1	511	9330	34105	42525	22827	5880	750     45	1];

    sum=0;
    for i = 1: n
        sum = sum + stirling2table(n+1, i+1)*(mu/eta)^i;
    end
    xv=(kv.^n) ./ (eta.^ kv); Ex= mean(xv);
    Ex_th=sum*exp(-mu*(1-1/eta));
    fprintf(' E[k^n / eta^k] for poisson k, n:%d, eta:%1.3f is [%d   %d] error:%d\n', n,eta, Ex, Ex_th, abs(Ex- Ex_th)/Ex);
    pause;



    %check lemma 2: E[x(k)=1(k~=0)/eta^k]: ([x(k)=1/eta^k  for k:poisson and nonzer], [x(k)=0 for k=0]);  %h(k)=!(k~=0)
    fprintf('\n\n*******************************************\n    Check Lemma 3.2 \n*******************************************\n');
    xv=1./ (eta.^ kv); xv(kv==0)=0; Ex= mean(xv);
    Ex_th=exp(-mu+mu/eta)-exp(-mu);
    fprintf(' E[1(k~=0) / eta^k] for poisson k, eta:%1.3f is [%d   %d] error:%d\n', eta, Ex, Ex_th, abs(Ex- Ex_th)/Ex);
    pause;



    %check lemma 3: E[x(k)] x(k)=h(k)*k^n/(eta^k)  for k:poisson and 
    fprintf('\n\n*******************************************\n    Check Lemma 3.3 \n*******************************************\n');
    xv=kv.^n./ (eta.^ kv); xv(kv==0)=0; Ex= mean(xv);

    sum=0;
    for i = 1: n
        sum = sum + stirling2table(n+1, i+1)*(mu/eta)^i;
    end
    Ex_th=exp(-mu+mu/eta)*sum;
    fprintf(' E[1 / eta^k] for poisson-like k, eta:%1.3f is [%d   %d] error:%d\n', eta, Ex, Ex_th, abs(Ex- Ex_th)/Ex);
    pause;
end

%%
if check_Ef
    %Check E[1/eta^kk]  for poisson like, kk~=0, used in energy efficiency
    fprintf('\n\n*******************************************\n    Check E[1/eta^kk], First part of Ef \n*******************************************\n');
    xv=1./ (eta.^ kkv); Ex= mean(xv);
    Ex_th=(exp(-mu+mu/eta)-exp(-mu)) /(1-exp(-mu));
    Ex_th=(exp(mu/eta)-1) /(exp(mu)-1); %the same as above
    
    fprintf(' E[1(k~=0) / eta^k] for poisson k, eta:%1.3f is [%d   %d] error:%d\n', eta, Ex, Ex_th, abs(Ex- Ex_th)/Ex);
    pause;
    fprintf('\n\n*******************************************\n    Check E[1/(kk*eta^kk)] : second part of Ef\n*******************************************\n');

    %Assumption 1: K: Poisson, eta=1
    xv=1./ (kkv .* eta.^ kkv); 
    Ex= mean(xv);  %average of 1/k for only positive k~Poiss(mu)

    EXPI = real(-expint(-mu/eta)); %Equal to : EXPI_mu_copy = -1i*pi - expint(-mu);
    Ex_th= (EXPI-log(mu/eta)-el)/(exp(mu)-1);
    
    %Wrong
    %EXPI_mu = real(-expint(-mu)); %Equal to : EXPI_mu_copy = -1i*pi - expint(-mu);
    %Ex_th= (EXPI_mu-log(mu))/(exp(mu/eta)-1);

    fprintf(' E[1(k~=0) / eta^k] for poisson k, eta:%1.3f is [%d   %d] error:[%d]\n', eta, Ex, Ex_th,  abs(2*(Ex- Ex_th)/(Ex+Ex_th)));
    pause;

    %ee    %pause;
    
end
if plot_expint  %plot E1 Ei relation
    %for E1(x) and Ei(x) see help expint
    x = -2:0.01:2; x=x*100
    E1x=expint(x);
    E1x_r=real(E1x); E1x_im=imag(E1x);
    E1x_rc=E1x; E1x_rc(x<0)=E1x(x<0) + 1i*pi; %imaginary part is -ipi
    figure; subplot(221); plot(x,E1x_r,'b'); hold on; plot(x(1:300:end),E1x_rc(1:300:end),'r*'); title('real part of: E1(x)=expint(x)');
    subplot(222); plot(x,E1x_im,'b'); title('imaginary part of: E1(x)=expint(x)');
    
    xp=0:0.001:2;
    E1xn=expint(-xp);Eix = -E1xn -1i*pi;
    Eix_c = real(-E1xn);
    
    subplot(223); plot(xp,Eix,'b'); hold on; plot(xp(1:300:end),Eix_c(1:300:end),'r*'); title('Ei(x)=real(-E1(-x))');
end

%%
save check_lemmas_lastrun