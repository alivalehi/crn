plot_ps=1;  %plot ps
plot_ps3d=1;
plot2d_pmf_s = 1; 
plot_anal_Es = 1;
plot_Ks=1;


if plot_Ks == 1
    clc; close all; clearvars -except plot_Ks plot_ps plot_ps3d plot_anal_Es plot2d_pmf_s
    sim1 = init_sim;
    sim1.N = 16;
    sim1.H = 40;
    sim1.Exp_PER = 0.2;
    sim1.NUM_RUNS = 5;
    sim1.use_framing_to_calc_Kp = 0;
    sim1.maxPack=sim1.maxPack*10;
    sim1.fname='TBP_N16_H40_PER20Test';


    %run the file again
    %plot_data_file('', 'TBP_N16_H40_PER20Test.mat',0)
    legend('Analytical: PER=0.2','Analytical: PER=0','Simulation:PER=0.2','Simulation:PER=0');
    ylabel('K_s = \sigma_s / E[s]');
    
end


if plot_ps == 1
    clc; close all; clearvars -except plot_Ks plot_ps plot_ps3d plot_anal_Es plot2d_pmf_s
    disp('Plot P(T:packet length) and p(S: number of symbols at packet:) ');
    lambda = 1; %symbol rate
    a=1/lambda; %interarrival time

    maxT = 20; %max length of packetization interval 
    maxK = 20; %max Number of bits at the last packetization interval
    maxRet = 20; %max number of retransmission  

    clrs={'b','g','r','y','m','c', 'b:','g:','r:','m:'};
    Tv= 5;% [0.2 0.5 1 2 5 10];

    H=20; N=8; b=2e-2; a=1-b;



    for tid = 1: length(Tv)
        pid=[]; lgnd={};
        T=Tv(tid); LT = lambda * T;   p0 = exp(-LT); q0=1-p0;

        pLen = zeros(1,maxT);
        for j = 1: maxT
            pLen(j) = p0^(j-1)*q0;  %p of zero symbols at j-1 time slot and unzero symbols at the last time slot
        end


        z = zeros(1,maxK); pKinT =z; pK=z; %K arrival in the last interval of length T

        kt=[0:.02:maxK]; %high resolution
        pKt = exp(-LT)*((LT).^kt)./gamma(kt+1) ./(1-p0); %envelope of pK

        kv = 1:maxK; Lv = H + kv*N; 
        PERv = 1-a.^Lv; PSRv=1-PERv;
        pKinT = (exp(-LT)*((LT).^kv)./factorial(kv));
        pK=pKinT/(1-p0); %pK = pKinT * (p(Len=1)+p(Len=2)+...) = pKinT*(1/(1-p0))

        fprintf('1: %s\n', num2str(pK)); pause;
        figure(tid);    stem(pK,[clrs{tid},'.']); hold on;    pid(tid) = plot(kt,pKt, [clrs{tid},':']);  axis([0 12 0 .5]);    lgnd{tid} = ['T\lambda=',num2str(T)]; title('P(service time) '); legend(pid,lgnd);    


        s=[];ps=[];
        for r = 1: maxRet
            s = [s; r*Lv];
            ps = [ps; pK .* ((PERv).^(r-1)) .* PSRv];
        end
        %check = [sum(ps); pK],    pause;


        tv=[0:0.01: maxRet]; st=[]; pst=[]; %envelope, smooth envelope for each k, continous retransmission
        for t = tv
            st=[st;t*Lv];
            pst = [pst; pK .* ((PERv).^(t-1)) .* PSRv];
        end

        %envelope, smooth envelope for each k, continous retransmission
        Lvt=H+kt*N;
        PERvt = 1-a.^Lvt; PSRvt=1-PERvt;
        skt=[]; pskt=[];
        for r = 1: maxRet
            skt = [skt; r*Lvt];
            pskt = [pskt; pKt .* ((PERvt).^(r-1)) .* PSRvt];
        end




        figure; pid=[]; lgnd={};
        for r = 1: 10%maxRet
            pid(t)=stem(s(r,:), ps(r,:),[clrs{r},'.']); hold on; lgnd{r}=['r=',r];
            if r==2, plot(skt(r,:),pskt(r,:), clrs{r}); end %envelope
        end
        pause;

        figure; plotdef(0); pid=[]; lgnd={};
        for k = 1: 10%maxK
            pid(k)= stem(s(:,k), ps(:,k),[clrs{k},'.']); hold on; lgnd{k}=['k=',k]; 
            if k==3, plot(st(:,k),pst(:,k), clrs{k}); end %envelope
        end
        legend(pid,lgnd);
        pause;

        sv=s(:);     psv=ps(:);
        figure; plotdef(0); pid=[]; lgnd={};
        stem(sv,psv,'g.', 'linewidth', 1); hold on;
        r=2;
            stem(s(r,:), ps(r,:),'b*', 'linewidth', 2); hold on; lgnd{1}=['r=',num2str(r)];
            pid(1)=plot(skt(r,:),pskt(r,:), 'b--', 'linewidth', 2); %envelope
        k=3;
            stem(s(:,k), ps(:,k),'ro', 'linewidth', 2); hold on; lgnd{2}=['k=',num2str(k)]; 
            pid(2)= plot(st(:,k),pst(:,k), 'r:', 'linewidth', 2); %envelope
        axis([0 600, 0 0.08]);
        xlabel('service time: s=r(kN+H)/R'); ylabel('propability mass function: f_s(s)');
        legend(pid, lgnd);
        %saveas(gcf,'Fig_ps','fig'); saveas(gcf,'Fig_ps','pdf'); saveas(gcf,'Fig_ps','png'); saveas(gcf,'Fig_ps','bmp'); saveas(gcf,'Fig_ps','eps'); 
    end
end
    
if plot_ps3d
    clc; close all; clearvars -except plot_Ks plot_ps plot_ps3d plot_anal_Es plot2d_pmf_s
end


if plot2d_pmf_s
    clc; close all; clearvars -except plot_Ks plot_ps plot_ps3d plot_anal_Es plot2d_pmf_s
    N= 3; H= 10; lambda = 10;
    clr = {'b.','g.','r.','y.','c.','k.','m.','b*','c*','m*','r*'};

    Kmax =20; Rmax = 100;

    k = [1: Kmax];     %number of symbols in a packet
    l = H+k*N


    r = [1:1:Rmax];   %number of retransmission;
    s = l' * r;
    

    figure;  stem3 (r,l,s,'b.');    pause;
    surf(r,l,s);    pause;

    figure;
    for i = 1:Kmax
        x=s(i,:);
        y=ones(1,Rmax);
        stem(x,y,clr{1+mod(i, 10)},'b.');
        hold on;
        pause;
    end;
end

if plot_anal_Es
    clc; close all; clearvars -except plot_Ks plot_ps plot_ps3d plot_anal_Es plot2d_pmf_s
    N= 3; H= 10; lambda = 10;
    clr = {'b.','g.','r.','y.','c.','k.','m.','b*','c*','m*','r*'};
    %mean service time
    figure;   
    T = 5/lambda;
    Plen_av = lambda * T * N + H;
    
    
    

    for PERav = 0;%[0.2];% 0.1]
        pb = 1- (1 - PERav)^(1/Plen_av);
        qb = 1-pb;
        
        nbits_in_T = Plen_av .* (1./(1-PERav));
        R = (nbits_in_T / T);
        
        
         Tv = [0.1/lambda, 1/lambda,   2/lambda:1/lambda:10/lambda, 20/lambda:10/lambda:100/lambda];
%       Tv = [2/lambda];%[0.1/lambda: 0.1/lambda : 1/lambda,   2/lambda:1/lambda:10/lambda];
        
        fprintf('H:%d   N:%d   lambda:%1.2f  T:%1.2f   TL:%d  pb:%d  Plen_av:%d R:%1.3f \n', ...
            H,N,lambda, T, T*lambda, pb, Plen_av, R);
        pause;
        
        
        LT = lambda * Tv;
        Ek_anal = lambda * Tv;
        Es_nopck0_anal =  N/(R * qb ^ H).* ((H/N+ LT./qb^N) .* exp(-LT*(1-1/qb^N))) - H/(R * qb ^ H)*exp(-LT);
        Es_dummy_anal =  N/(R * qb ^ H).* ((H/N+ LT./qb^N) .* exp(-LT*(1-1/qb^N)));

       
%% One Notation        
%Something might be wrong        
BB = qb;
A1 = (2*(N^2)*LT) / (BB^(2*N+2*H)) + (2*(N^2)*LT.^2) / (BB^(4*N+2*H)) ...
    +(4*H*N*LT) / (BB^(2*N+2*H)) + (2*H^2)/(BB^(2*H));
A2 = (N^2)*LT / BB^(N+H) + (N^2)*LT.^2 / BB^(2*N+H) ...
    +2*H*N*LT / BB^(N+H) + H^2/BB^(H);
A3 = 2*H^2 / BB^(2*H) - H^2 / BB^H;
%Es2_nopck0_anal = (1./R^2) .* (A1 .* exp(-LT*(1-1/BB^(2*N))) - A2 .* exp(-LT*(1-1/BB^N)) + A3 .* exp(-LT));



%% %Different notation
% Accurate Result

AN2_1 = (2*(N^2)/BB^(2*H))*(LT*BB^(-2*N) + LT.^2 *BB^(-4*N)) .* exp(-LT*(1-BB^(-2*N))) ;
AN2_2 = -((N^2)/BB^(H))*(LT*BB^(-N) + LT.^2 *BB^(-2*N)) .* exp(-LT*(1-BB^(-1*N))) ;
AN2 = (AN2_1+AN2_2) / R^2;

AKHN_1 =  (4*H*N/BB^(2*H))*(LT*BB^(-2*N)) .* exp(-LT*(1-BB^(-2*N))) ;
AKHN_2 = -(2*H*N/BB^(H))*(LT*BB^(-N)) .* exp(-LT*(1-BB^(-1*N))) ;
AKHN = (AKHN_1 + AKHN_2) / R^2;

AH2_1 = (2*H^2/BB^(2*H)) .* (exp(-LT*(1-BB^(-2*N)))-exp(-LT)) ;
AH2_2 = -(1*H^2/BB^(1*H)) .* (exp(-LT*(1-BB^(-1*N)))-exp(-LT));
AH2 = (AH2_1 + AH2_2) / R^2;
Es2_nopck0_anal = AN2 + AKHN + AH2;
%%

Es2_dummy_anal =  (A1 .* exp(-LT*(1-1/BB^(2*N))) - A2 .* exp(-LT*(1-1/BB^N)))./(R^2) ;
        
Es_ifkp_nopck0_anal = Es_nopck0_anal ./ (1-exp(-LT));        
Es2_ifkp_nopck0_anal = Es2_nopck0_anal ./ (1-exp(-LT));        


testEs2_TL_approach0 = (1/R^2)*(H+N).^2* (2/(BB^(2*(H+N))) - 1/(BB^(H+N)));

        %simulation
        for i = 1: length(Tv)
            fprintf('.');
            T = Tv(i);
            
            maxPn = 20000;  %20000 
            maxSn = 100000; %100000
            minSn = 10000;
            
            xlen = max(minSn, min(maxSn,  floor(maxPn*T*lambda)));
            
            x = random('exp', (1/lambda),[1, xlen]);  %interarrivals
            [tp, kp, fiv] = f_perform_framing(1, x, T, 0, (1/lambda), 0, 0);
            
             Ek_sim(i) = mean(kp);

            [wnvd, svd, rvd] = f_q_evolution(tp, kp, N, H, pb, R,1); %1: send_dummy_for_k0
            
           
nzind = find(kp ~= 0);
tp = tp(nzind);
kp = kp(nzind);

[wnvnd, svnd, rvnd] = f_q_evolution(tp, kp, N, H, pb, R,0); %0: no_dummy_for_k0
% mean(wnvnd)
% pause;
% 
% [mean(tp(2:end)-tp(1:end-1)) mean(svnd)]
% pause;
            
            
            
            
            Es_dummy_sim(i) = mean(svd);
            Es2_dummy_sim(i) = mean(svd.^2);% - ((mean(sv)).^2);
            
            Es_nopck0_sim(i) = mean(svnd);
            Es2_nopck0_sim(i) = mean(svnd.^2);% - ((mean(sv)).^2);
            
            svp = svnd(svnd ~= 0);
            Es_ifkp_nopck0_sim(i) = mean(svp);
            Es2_ifkp_nopck0_sim(i) = mean(svp.^2);% - ((mean(sv)).^2);

            
        end

        
        figure;
        plot(Tv*lambda, Es_nopck0_anal, 'b', 'linewidth', 2);
        hold on;
        plot(Tv*lambda, Es_dummy_anal, 'g', 'linewidth', 2);
        plot(Tv*lambda, Es_ifkp_nopck0_anal, 'r', 'linewidth', 2);
        
        plot(Tv*lambda, Es_nopck0_sim, 'b--', 'linewidth', 2);
        plot(Tv*lambda, Es_dummy_sim, 'g--', 'linewidth', 2);
        plot(Tv*lambda, Es_ifkp_nopck0_sim, 'r--', 'linewidth', 2);
        legend('Anal: No pack for 0', 'Anal: Dummy pack for 0', 'Anal: Es2 for k>0', 'Simul: No pack for 0', 'Simul: Dummy pack for 0', 'Sim: Es for k>0');
        ylabel('E[s]');
        xlim([0,6]);
        pause;
        

        figure;
        plot(Tv*lambda, Es2_nopck0_anal, 'b', 'linewidth', 2);
        hold on;
        plot(Tv*lambda, Es2_dummy_anal, 'g', 'linewidth', 2);
        plot(Tv*lambda, Es2_ifkp_nopck0_anal, 'r', 'linewidth', 2);
       
        plot(Tv*lambda, Es2_nopck0_sim, 'b--', 'linewidth', 2);
        plot(Tv*lambda, Es2_dummy_sim, 'g--', 'linewidth', 2);
        plot(Tv*lambda, Es2_ifkp_nopck0_sim, 'r--', 'linewidth', 2);
        legend('Anal: No pack for 0', 'Anal: Dummy pack for 0', 'Anal: Es2 for k>0', 'Simul: No pack for 0', 'Simul: Dummy pack for 0', 'Sim: Es2 for k>0');
        ylabel('E[s2]');
        xlim([0,6]);
        
    end
end
