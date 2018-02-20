clc; clear all; close all
fprintf('This program generates plots for paper. \n Note 1-some of the plots maybe in plot_results. \n Note 2-set flags in the code for desired plots\n'); 

clrs = {'b','g','r','y','c','k','m','b-.','g-.','r-.','y-.','c-.','k-.','m-.','b--','g--','r--','y--','c--','k--','m--',...
    'b-*','g-*','r-*','y-*','c-*','k-*','m-*','b-x','g-x','r-x','y-x','c-x','k-x','m-x'};


run_comp_per = 0;    %plot E[d] for different PER choices
plot_comp_per = 0;


run_snr = 0;     %plot Topt versus SNR for different PER choices [1,3:uncoded system] and [2,3:coded system] and a given H/N and Rch [
plot_snr = 1;

run_comp_H_coded = 0;     %plot Topt versus H/N for different PER choices [1,3:uncoded system] and [2,3:coded system]
plot_comp_H_coded = 0;

run_comp_H_uncoded = 0; %for uncoded different PER
plot_comp_H_uncoded = 0;


run_compare_Ks_analytical_simulation=0;   %1 runs for uncoded system, 2 runs for coded system  3:runs for both systems
plot_compare_Ks_analytical_simulation=0;   %plots main  results [compares analysis to MC simulation for K[s] and E[w] and E[d]]

plot_Ef=0;     %Plot Energy Efficiency Curves



%%

if run_snr
    
    N=8; H=40;

    load('BER_r60', 'nRuns', 'EsNodBv', 'uncodedBERVec', 'convBERVec')  
    mkdir dataSNR

    figure; semilogy(EsNodBv, uncodedBERVec, 'b'); hold on; semilogy(EsNodBv, convBERVec{1}, 'g');  semilogy(EsNodBv, convBERVec{2}, 'r');  legend('Uncoded', 'Coded[1]', 'Coded[2]'); xlabel('SNR'); ylabel('BER'); pause;
    
    
    BERimpCA = uncodedBERVec ./convBERVec{1}(1,:); %conv rate 1/2
    BERimpCB = uncodedBERVec ./convBERVec{2}(1,:); %conv rate 3/4
    
    snr_ind_v1 = find(uncodedBERVec>0); 
    snr_ind_v2 = find(BERimpCA>1 & isfinite(BERimpCA)  & (convBERVec{1}(1,:)>0));%find(isfinite(BERimpCA));
    snr_ind_v3 = find(BERimpCB>1 & isfinite(BERimpCB)  & (convBERVec{2}(1,:)>0));%find(isfinite(BERimpCA));

    clear sims
    for i = sort(unique([snr_ind_v1,snr_ind_v2, snr_ind_v3])) 
        
        sim = init_sim; 
        sim.Rch=  400;
        sim.overwrite_active=1; sim.plot_active=0;
        sim.datapath = [cd , '\dataSNR\'];
        sim.Tv  = (1/sim.lambda)*[0.001:0.001:0.01, 0.02:0.01:10, 10.1:0.1: 30, 31:1:200]; sim.Tv2=sim.Tv;
            

        sim.lambda = 10;   
        sim.T_TB = 5 * (1/sim.lambda);    
        sim.N = N;
        sim.H = H;
            
        sim.BER = uncodedBERVec(i); 
        sim.simulation = 0;  sim.analysis = 1; 
        sim = calc_sim_params(sim); 
        
        if ismember(i, snr_ind_v1)
            simU = sim; simU=fqueu_TB(simU); simsU{i} = simU;
        end
        
        if ismember(i, snr_ind_v2)
            simC1=sim;   simC1.CodedSystem = 1; 
            simC1.impBERD = BERimpCA(i); simC1.impBERH = BERimpCA(i);  %reduction in BER  effectiveBER = BER/BERimp
            simC1.CTD = 1; simC1.CTH=1; 
            simC1.CRD = 1/2; simC1.CRH=1/2;      %coding rate
            simC1.nTbitsD =3; simC1.nTbitsH =3;  %number of Tail Bits
            simC1=fqueu_TB(simC1); simsC1{i} = simC1;
        end
        
        if ismember(i, snr_ind_v3) 
            simC2=sim; 
            simC2.CodedSystem = 1; 
            simC2.impBERD = BERimpCB(i); simC2.impBERH = BERimpCB(i);  %reduction in BER  effectiveBER = BER/BERimp
            simC2.CTD = 1; simC2.CTH=1;          %0:uncoded 1:covEnc  2:Turbo  3:LDPC
            simC2.CRD = 3/4; simC2.CRH=3/4;      %coding rate
            simC2.nTbitsD =2; simC2.nTbitsH =3;  %number of Tail Bits
            simC2=fqueu_TB(simC2); simsC2{i} = simC2;
        end
        
        if ismember(i,snr_ind_v2) && ismember(i, snr_ind_v3)
            simC3=sim; 
            simC3.CodedSystem = 1; 
            simC3.impBERD = BERimpCB(i); simC3.impBERH = BERimpCA(i);  %reduction in BER  effectiveBER = BER/BERimp
            simC3.CTD = 1; simC3.CTH=1;          %0:uncoded 1:covEnc  2:Turbo  3:LDPC
            simC3.CRD = 3/4; simC3.CRH=1/2;      %coding rate
            simC3.nTbitsD =2; simC3.nTbitsH =3;  %number of Tail Bits
            simC3=fqueu_TB(simC3); simsC3{i} = simC3;
        end
    end
    save ('data_SNR.mat');
end

if plot_snr
    load ('data_SNR.mat','simU', 'simsU', 'simsC1', 'simsC2','simsC3', 'H', 'N',  'EsNodBv', 'snr_ind_v1', 'snr_ind_v2', 'snr_ind_v3');
    clrs = {'b','g','r','c','m';    'b:','g:','r:','c-o','m-o';    'b--','g--','r--','c--','m--'}; %clrs=clrs';
    
    pid=0; p=[]; lg={};
    nSNR=max(sort(union(snr_ind_v1,snr_ind_v2)) );
    ToptU=nan(1,nSNR); ToptC1=nan(1,nSNR); ToptC2=nan(1,nSNR); ToptC3=nan(1,nSNR);
    EdU=nan(nSNR, length(simU.Tv2)); EdC1=EdU; EdC2=EdU; EdC3=EdU;
    k = 0;
    for j=sort(unique([snr_ind_v1,snr_ind_v2, snr_ind_v3]))
        k = k+1;
        if exist('simsU', 'var') && ismember(j, snr_ind_v1)
            simP = simsU{j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda; EdU(j,:)=Ed; Ed(~isfinite(Ed))=+1e10; 
            [minval, ind] = min(Ed); if ind > 1, ToptU(j) = Tv(ind); else ToptU(j)=nan; end
        end
        if exist('simsC1', 'var') && ismember(j, snr_ind_v2)
            simP = simsC1{j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda;  EdC1(j,:)=Ed; Ed(~isfinite(Ed))=+1e10;
            [minval, ind] = min(Ed);  if ind > 1, ToptC1(j) = Tv(ind); else ToptC1(j)=nan; end
        end
        if exist('simsC2', 'var') && ismember(j, snr_ind_v3)
            simP = simsC2{j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda; EdC2(j,:)=Ed; Ed(~isfinite(Ed))=+1e10; 
            [minval, ind] = min(Ed); if ind > 1, ToptC2(j) = Tv(ind); else ToptC2(j)=nan; end
        end
        if exist('simsC3', 'var') && (ismember(j, snr_ind_v2) && ismember(j, snr_ind_v3))
            simP = simsC3{j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda; EdC3(j,:)=Ed; Ed(~isfinite(Ed))=+1e10; 
            [minval, ind] = min(Ed); if ind > 1, ToptC3(j) = Tv(ind); else ToptC3(j)=nan; end
        end

    end
    if 0
        figure; plotdef(0); 
        p(1) = plot(EsNodBv(snr_ind_v1), ToptU(snr_ind_v1), 'b'); hold on;
        p(2) = plot(EsNodBv(snr_ind_v2), ToptC1(snr_ind_v2), 'g'); hold on;
        p(3) = plot(EsNodBv(snr_ind_v2), ToptC2(snr_ind_v2), 'c'); 
        p(4) = plot(EsNodBv(snr_ind_v2), ToptC3(snr_ind_v2), 'r');
        legend(p, {'UnCoded', 'CS1', 'CS2', 'CS3'});
        xlabel('Channel SNR (dB)');
        ylabel('Optimum Packetization Time: D \lambda');
        title('');  %saveas(gcf, 'SNR_TOPT', 'fig');   
        pause;
    end
    
    for col = [1059]
        close; 
        figure; plotdef(0); snrvplot=EsNodBv([1:nSNR]); 
        EdUp = EdU(:,col); EdC1p = EdC1(:,col)'; EdC2p= EdC2(:,col); EdC3p = EdC3(:,col)'; 
        %smth='rlowess'; EdC1p= smooth(EdC1(:,col)', smth)'; EdC2p= smooth(EdC2(:,col)', smth)'; EdC3p = smooth(EdC3(:,col)', smth)';
        p(1) = plot(snrvplot,EdUp, 'b'); hold on;   m1=min(EdUp); plot(snrvplot, m1*ones(1,length(snrvplot)), 'b--', 'linewidth',1.5);
        p(2) = plot(snrvplot, EdC1p, 'g'); hold on; m1=min(EdC1p); plot(snrvplot, m1*ones(1,length(snrvplot)), 'g--','linewidth',1.5);
        p(3) = plot(snrvplot, EdC2p, 'r'); hold on; m1=min(EdC2p); plot(snrvplot, m1*ones(1,length(snrvplot)), 'r--', 'linewidth',1.5);
        p(4) = plot(snrvplot, EdC3p, 'c'); hold on; m1=min(EdC3p); 
        
        
        legend(p, {'UnCoded', 'CS1', 'CS2', 'CS3'});
        xlabel('Channel SNR (dB)');  axis([-6 12 0 120])
        ylabel('Expected End-to-End Delay: E[d]');
        %title(['col:',num2str(col)]);  saveas(gcf, ['FigCodedEd_col', num2str(col)], 'png');   
        saveas(gcf, 'FigCodedEd', 'png');   
    end
end


%%
if run_comp_H_coded
    mkdir dataH
    NUM_RUNS = 10; N = 8;    
    Hv = [1:1:50,60:5:150];  BERv = [0.001 0.005 .01];  %For Uncoded and coded system

    clear sims
    for i = 1:length(BERv)   
        for j=1:length(Hv)  
            sim = init_sim; 
            sim.Rch=400;     
            sim.overwrite_active=1; sim.plot_active=0;
            sim.datapath = [cd , '\dataH\'];
            sim.Tv  = (1/sim.lambda)*[0.001:0.001:0.01, 0.02:0.01:10, 10.1:0.1: 30, 31:1:200]; sim.Tv2=sim.Tv;
            

            sim.lambda = 10;   
            sim.T_TB = 5 * (1/sim.lambda);    
            sim.N = N;
            sim.H = Hv(j);
            
            sim.BER = BERv(i); 
            sim.simulation = 0;  sim.analysis = 1; 
            sim = calc_sim_params(sim); 
            
            %run for uncoded system
            simU = sim; simU=fqueu_TB(simU); simsU{i,j} = simU;
                
            %run for coded system
            simC1=sim; 
            simC1.CodedSystem = 1; 
            simC1.impBERD = 10; simC1.impBERH = 10;  %reduction in BER  effectiveBER = BER/BERimp
            simC1.CTD = 1; simC1.CTH=1; 
            simC1.CRD = 1/2; simC1.CRH=1/2; %coding rate
            simC1.nTbitsD =2; simC1.nTbitsH =2;  %number of Tail Bits
            simC1=fqueu_TB(simC1); simsC1{i,j} = simC1;

            simC2=sim; 
            simC2.CodedSystem = 1; 
            simC2.impBERD = 10; simC2.impBERH = 100;  %reduction in BER  effectiveBER = BER/BERimp
            simC2.CTD = 1; simC2.CTH=1; 
            simC2.CRD = 1/2; simC2.CRH=1/3; %coding rate
            simC2.nTbitsD =2; simC2.nTbitsH =3;  %number of Tail Bits
            simC2=fqueu_TB(simC2); simsC2{i,j} = simC2;
        end
    end
    save ('data_comp_H_coded.mat');
end

if plot_comp_H_coded
    %load('data_comp_H_bak.mat')
    load ('data_comp_H_coded.mat','simsU', 'simsC1', 'simsC2', 'BERv',  'Hv');
    clrs = {'b','g','r','c','m';    'b:','g:','r:','c-o','m-o';    'b--','g--','r--','c--','m--'}; %clrs=clrs';
    
    k = 0;
    figure; plotdef(0); pid=0; p=[]; lg={};
    for i = [1, 3]%length(BERv)%
        ToptU=nan(1,length(Hv)); ToptC1=nan(1,length(Hv)); ToptC2=nan(1,length(Hv));
        for j=1:length(Hv)  
            k = k+1;
            simP = simsU{i,j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda; Ed(~isfinite(Ed))=+1e10; 
            [minval, ind] = min(Ed); if ind > 1, ToptU(j) = Tv(ind); else ToptU(j)=nan; end
            simP = simsC1{i,j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda; Ed(~isfinite(Ed))=+1e10; 
            [minval, ind] = min(Ed);  if ind > 1, ToptC1(j) = Tv(ind); else ToptC1(j)=nan; end

            simP = simsC2{i,j}; Tv = simP.Tv2 * simP.lambda;
            Ed = simP.Ed* simP.lambda; Ed(~isfinite(Ed))=+1e10; 
            [minval, ind] = min(Ed); if ind > 1, ToptC2(j) = Tv(ind); else ToptC2(j)=nan; end
            
        end
        pid=pid+1; p(pid) = semilogy(Hv, ToptU, clrs{1,i}); lg{pid} = sprintf('UnCoded: BER:%1.1d', BERv(i)); hold on;   
        pid=pid+1; p(pid) = semilogy(Hv, ToptC1, clrs{2,i}); lg{pid} = sprintf('CS[1]'); hold on;
        pid=pid+1; p(pid) = semilogy(Hv, ToptC2, clrs{3,i}); lg{pid} = sprintf('CS[2]');
    end
    legend(p, lg); xlim([0 60]);   xlabel('Number of Header Bits');    ylabel('Optimum Packetization Time: D \lambda'); xlim([0 100]);
    title('');    %saveas(gcf, 'PER', 'fig');    
end

%%
if run_comp_H_uncoded
    mkdir dataH2
    NUM_RUNS = 10; N = 8;    
    
    Hv = [1:1:50,60:10:150]; Exp_PERv = [0.2 0.4 0.5 0.6]; %Use for uncoded 
    
    clear sims
    for i = 1:length(Exp_PERv)    %for i = 1:length(Exp_PERv)
        
        for j=1:length(Hv)  
            sim = init_sim; 
            
            sim.overwrite_active=1; sim.plot_active=0;
            sim.datapath = [cd , '\dataH2\'];
            sim.Tv  = (1/sim.lambda)*[0.001:0.001:0.01, 0.02:0.01:10, 10.1:0.1: 30, 31:1:200]; sim.Tv2=sim.Tv;

            sim.lambda = 10;   
            sim.T_TB = 5 * (1/sim.lambda);    
            sim.N = N;
            sim.H = Hv(j);
            
            sim.Exp_PER = Exp_PERv(i);
            sim.simulation = 0;  sim.analysis = 1; 
            sim = calc_sim_params(sim); 
            
            %run for uncoded system
            simU = sim; simU=fqueu_TB(simU); simsU{i,j} = simU;
        end
    end
    save ('data_comp_H_uncoded.mat');
end

if plot_comp_H_uncoded
    load ('data_comp_H_uncoded.mat','simsU', 'Exp_PERv',  'Hv');
    clrs = {'b','g','r','c','m';    'b:','g:','r:','c-o','m-o';    'b--','g--','r--','c--','m--'}; %clrs=clrs';
    
    k = 0;
    figure; plotdef(0); pid=0; p=[]; lg={};
    for i = 1:length(Exp_PERv)
        ToptU=nan(1,length(Hv)); ToptC1=nan(1,length(Hv)); ToptC2=nan(1,length(Hv));
        for j=1:length(Hv)  
            k = k+1;
            if  exist('simsU', 'var')
                simP = simsU{i,j}; Tv = simP.Tv2 * simP.lambda;
                Ed = simP.Ed* simP.lambda; Ed(~isfinite(Ed))=+1e10; 
                [minval, ind] = min(Ed); if ind > 1, ToptU(j) = Tv(ind); else ToptU(j)=nan; end
            end
           
        end
        if  exist('simsU', 'var')
            pid=pid+1; p(pid) = semilogy(Hv, ToptU, clrs{1,i}); lg{pid} = sprintf('UnCoded: PER:%1.1d',Exp_PERv(i)); hold on;   
        end
 
    end
    legend(p, lg); xlim([0 30]);
    xlabel('Number of Header Bits');    ylabel('Optimum Packetization Time: D \lambda');
    title('');    %saveas(gcf, 'PER', 'fig');    
end


%%
if run_comp_per
    NUM_RUNS = 10;
    N = 16;
    H = 30;
    Exp_PERv = [0, 0.12, 0.14, 0.15];    
    BERv = [0, 5e-4, 8e-4, 1e-3];%   , 2e-3 5e-3];
    
    clear simsper
    
    pi = 0;
    for i=1:length(BERv)%(Exp_PERv)  
        simr = init_sim;
        
        simr.lambda = 10;   
        simr.T_TB = 5 * (1/simr.lambda);    
simr.Tv = (1/simr.lambda)*[[1:9]*1e-4, [1:9]*1e-3, 0.01:0.01:100];%[1 2 5]*1e-2, 0.1, 0.2, 0.5, 0.8, 1:0.5:30, 31:1:40, 45:5:70, 80:10:100, 150, 200]; %for BER use
simr.Tv2 = simr.Tv;
        simr.N = N;
        simr.H = H;
        simr.Exp_PER = [];%Exp_PERv(i);
        simr.BER = BERv(i);
        simr.simulation = 0;
        simr.analysis = 1;
        simr.save_active=0;
simr.Rch = 250;
        simr=fqueu_TB(simr);
        simsper(i) = simr;
     end
     save ('data_comp_BER','simr','simsper', 'Exp_PERv', 'BERv');
 end
 if  plot_comp_per
    load('data_comp_BER','simr','simsper','BERv');
    clrs = {'k','b','r','g','c','m', 'y'};
    %clrs2 = {'k-.','b-o','r-x','g-*','k-S','m', 'k','k:','m:','m--','c','c:','c--','b:','b--','g:','g--','r:','r--'};
    figure; plotdef(0);
    for i=1:length(BERv)%(Exp_PERv)  
        simr = simsper(i);
        fprintf('I:%d, BER:%1.4f  Rch:%d   Rbv:%s\n',i, simr.BER, simr.Rch, num2str(simr.Rbv));
        Tv = simr.Tv2 * simr.lambda;
        Ed = simr.Ed * simr.lambda;
        
        
        [a] = find(isnan(Ed)==0);
        if isempty(a), Ed(1:end)=1e5; else %all NAN 
            ind = a(1) - 1;
            if (ind > 0) && (ind < length(Ed))
                Ed(ind) = 1e5;
            end
        end        
        
        p(i) = plot(Tv, Ed, clrs{i}); hold on;
        lgnd{i} = sprintf('BER:%0.0g', simr.BER);  
       
    end
    legend(p, lgnd);    xlim([0 20]);     ylim([0 60])
    xlabel('Packetization Time: T \lambda');
    ylabel('Average Delay: E[d] \lambda');
    titl=sprintf('Average End to End delay for N=%d H=%d \\lambda:%d', simr.N, simr.H, simr.lambda);
    title(titl);    title('');    %saveas(gcf, 'PER', 'fig');    saveas(gcf, 'PER', 'pdf');
end



if run_compare_Ks_analytical_simulation
    mkdir data
    disp(' This is a test call with given params to generate comparsion between simulation and analysis');
    disp('The following are the parameters for uncoded system');
    
    sim1 = init_sim;

    sim1.minSym = 10 * 50000;
    sim1.maxSym = 10* 500000;
    sim1.maxPack = 10* 50000;


    sim1.plot_active = 1;
    sim1.N = 16;
    sim1.H = 40;
    sim1.Exp_PER = 0.2;
    sim1.NUM_RUNS = 100;
    sim1.use_framing_to_calc_Kp = 0;
    sim1.Tv = [0.1, 0.5, 1:1: 3, 3.1:0.1: 20, 21:1:30, 35:5:50, 60:10:100] * sim1.a;
    sim1.fname='Es_Ew_comp.mat';
    sim1.save_active=1;
    sim1.analysis = 1;
    sim1.simulation = 1;
    fqueu_TB(sim1);
    
end 
if plot_compare_Ks_analytical_simulation
    
    plot_data_file('data\', 'Es_Ew_comp.mat', 0); xlim([0 2]);
    
    load data\Es_Ew_comp
    plot_analysis = 1;
    plot_simulation = 1;
    if plot_analysis
        figure; plotdef(0);
        plot(sim.Tv2, sqrt((sim.Es2)-((sim.Es).^2)) ./ sim.Es, 'm', 'linewidth', 1); hold on;
        plot(sim.Tv2, sqrt((sim.Es2_0)-((sim.Es_0).^2)) ./ sim.Es, 'g', 'linewidth', 1); hold on;
        pause;
    end
    if plot_simulation
        simsigEs = sqrt(sim.s2av(lpi,:)-sim.sav(lpi,:).^2);
        simcoeffEs = simsigEs./ sim.sav(lpi,:);
        simcoeffEs_0 = sqrt((sim.Es2_0)-((sim.Es_0).^2)) ./ sim.Es;


        %downsample
        ind = [1,3, 5, 25, find(mod(sim.Tv * 2, 2) == 0)];
        Tvds = sim.Tv(ind);
        simcoeffEs = simcoeffEs(ind);
        simcoeffEs_0 = simcoeffEs_0(ind);

        plot(Tvds, simcoeffEs, 'm--o', 'linewidth', 2, 'markersize', 8); hold on
        plot(Tvds, simcoeffEs_0, 'g--*','linewidth', 2 ,'markersize', 8); hold on; %we are sure that there is a perfect match between analythical and theory
    end

    xlim([0,5]);    ylim([0,1.5]);    ylim([0,1]);    
    xlabel('Packetization Time: T');
    legend('Analytical: PER=0.2','Analytical: PER=0','Simulation:PER=0.2','Simulation:PER=0');
    ylabel('K_s = \sigma_s / E[s]');
    pause;  %         saveas(gcf, 'Ks2', 'fig');         saveas(gcf, 'Ks2', 'pdf');
end
    


if plot_Ef
mu = 4; Flen = 1e6; kmax = max(10,5*mu);
N = 8; H=40;  %bits and header
el = 0.5772156649; % Euler–Mascheroni constant
   
run_simulation=1;
use_Er=0;   use_average_r=0; 
    lambda = 1;
    Tv = (1/lambda) * [0:0.001:0.009, 0.01:0.01:20];
    Tvp = (1/lambda)* [0.001:0.001:0.009, 0.01:0.01:0.09, 0.1:0.1:0.9, 1: 0.5: 20];  %for simulations
    Tvp = (1/lambda)* [0: 0.5: 20];  %for simulations
    muv = lambda * Tv; 
    fe = figure; plotdef(0);%figure(2); plotdef(0);
    pid=[]; lgnd={}; i=0;
    for beta = [0.02 0.01  0.002 0]%[0.02 0.01 0.005 0.002 0.001 0.0001]
        alpha=1-beta;
        eta= alpha^N;    i=i+1; 


        Ef1 = (1/alpha^H) * (exp(muv/eta)-1) ./(exp(muv)-1);
        Ef2 = (1/alpha^H) * (H/N) *(real(-expint(-muv/eta))-log(muv/eta)-el)./(exp(muv)-1);
        Ef=Ef1+Ef2; scale_factor=min(Ef);
        Efs = Ef ./ scale_factor;  %normalize
        
        %%
        if run_simulation %&& (alpha == 0.99) %calculate experimental
            Flen = 1e3; 
            
            j=0;   fprintf('Simulation %d/%d: ', i, length(alpha));
            for Tp = Tvp
                fprintf('.');
                mup = lambda * Tp; j=j+1;
                k = poissrnd (mup,1,Flen);
                kk = k(k>0); %remove zeros
                PSR = alpha .^ (N * kk + H); %packet sucess rate
                Er_th = 1./PSR; %E[r|k] 

                if use_Er  %approximation
                    r = Er_th;
                elseif use_average_r
                    r2=zeros(100,length(PSR));
                    for nr = 1: 100
                       r2(nr,:) = geornd(PSR) + 1;
                    end
                    Er_p = sum(r2./100);
                    fprintf('\nE[r]  Analytical: %s \n      Practical: %s \n', num2str(round(Er_p)), num2str(round(Er_th)) ); %[Er_p; Er_th]
                    fprintf('Error : %s \n', num2str(abs(2*(Er_p-Er_th)./(Er_p+Er_th))));

                    r=Er_p;
                else
                    r = geornd(PSR) + 1;  %in matlab starts with 0 P(1-P)^x
                end
                
                Efp(j) = mean(r .* (kk * N + H) ./ (kk*N));
            end
            %remove nan, inf
            Tvpp = Tvp(isfinite(Efp));
            Efpp = Efp(isfinite(Efp));
            
            fprintf('\n');
            %normalize with the same factor
            
            %Efpps=Efpp/scale_factor;
            %num2str([Ef; Efp]), pause;
            figure(fe)
            semilogy(Tvpp, Efpp,[clrs{i},'-.'], 'linewidth',3); hold on;
            
        end
        
        figure(fe);
        pid(i)=semilogy(Tv, Ef,clrs{i}); hold on; lgnd{i}=['BER: ',num2str(1-alpha)];
        
        %find optimal T
        [minEf, minind] = min(Ef); optTv=Tv(minind); optmu=lambda*optTv;
        if minind < length(Tv), semilogy([optTv optTv], [1e-2 minEf], [clrs{i},'--^'], 'markersize',10); end
        
        
    end
    
    
    plot([5 5],[0.01 100], 'k--'); plot([10 10],[0.01 100], 'k--');
    %ylim([0.1 100]); 
    xlabel('Packetization Time: \lambda T'); ylabel('Normalized Energy Per Bit');
    legend(pid,lgnd); ylim([0.9 20]);
    saveas(gcf,'EF','fig');    saveas(gcf,'EF','png');    saveas(gcf,'EF','eps');    saveas(gcf,'EF','pdf');
end

%%


