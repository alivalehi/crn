%This program runs simulations for different parameter setups
%see also call_framing_delay for similar test setups

%clear all; close all; clc;
quick_run=1;
if ispc && quick_run, disp('quick run activated !!!!');
    %     pause;
end
disp('b');
%#######################################################################
%     SIMULATION OPTIONS
%#######################################################################
%multicase result
if 0
    sim = init_sim(2);
    % sim.SU =[1 5 10 50];
    sim.SU =1;
    sim.CodedSystem=0;  sim.cognitive=1;
    sim.control.simulation = 0;
    %  sim.ch.Rch = 10000;
    sim.ch.BER= 8e-3;
    sim.trafficmode = 1;
    
    if(sim.trafficmode == 1) %heavy
        sim.ch.Rch = 1000;
        sim.lambda = 5;
    else
        sim.ch.Rch = 1000;
        sim.lambda = 1;
    end
    BER_vect = [0   2e-3:1e-3:1e-2  ];
    %BER_vect = [ 2e-2:1e-2:1e-1 2e-1:1e-1:5e-1];
    % BER_vect = [0 1e-8:1e-8:1e-7  1.1e-7:1e-8:1e-6  1.1e-6:1e-7:1e-5  1.1e-5:1e-6:1e-4  1.1e-4:1e-5:1e-3  1.1e-3:1e-4:1e-2   1.1e-2:1e-3:1e-1 1.1e-1:1e-2:5e-1];
    H_vect = [70];
    N_vect = [5:2:20];
    Rch_vect = [1000 1500 1000 1500 500 500];
    lambda_vect = [10 15 30 60 10 20];
    chNotAvRatio_vect = [1 2 1 2 10];
    chAvPlenRatio_vect = [1 2 2 1 10];
 %   for i = 1: length(BER_vect)
        for j = 1: length(N_vect)
            %for k = 1: length(Rch_vect)
             %   for l = 1: length(chNotAvRatio_vect)
                    try
                      %  BER=BER_vect(i);
                        sim.H = H_vect(1);
                        sim.N = N_vect(j);
                      %  sim.ch.Rch = Rch_vect(k);
                      %  sim.lambda = lambda_vect(k);
                      %  sim.CRN.chNotAvRatio = chNotAvRatio_vect(l);
                      %  sim.CRN.chAvPlenRatio = chAvPlenRatio_vect(l);
                        if quick_run
                            sim.control.plot_active = 1;
                            sim.control.debug_active = 0;
                            sim.control.runtime_plot_active=1;
                            
                            sim.run.NUM_RUNS = 4;
                            %run.quick_run_coef = 0.1;
                            sim.Kv = [1:1:30];
                        else
                            sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
                            sim.run.NUM_RUNS = 4;
                            sim.run.quick_run_coef = 1;  %only for very quick test
                            sim.Kv = [1:1:10,    12:2:20, 25:5:50];
                        end
                        sim = fqueu_NB(sim);
                    catch
                        % Nothing to do
                    end
               % end
            %end
        end
%    end
    
end
%Uncoded basic system
if 1
    sim = init_sim(2);
    %sim.SU =[15 30 45];
    sim.SU =1;
    sim.CodedSystem=0;  sim.cognitive=1;
    %  sim.ch.Rch = 10000;
    sim.ch.BER= 9e-3;
    sim.trafficmode = 1;
    sim.H = 40;
    sim.N = 4 ;
    if(sim.trafficmode == 1) %heavy
        sim.ch.Rch = 700;
        sim.lambda = 5;
    else
        sim.ch.Rch = 1000;
        sim.lambda = .1;
    end
    if quick_run
        sim.control.plot_active = 1;
        sim.control.debug_active = 0;
        sim.control.runtime_plot_active=1;
        
        sim.run.NUM_RUNS = 4;
        %run.quick_run_coef = 0.1;
        sim.Kv = [1:1:30];
    else
        sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
        sim.run.NUM_RUNS = 4;
        sim.run.quick_run_coef = 1;  %only for very quick test
        sim.Kv = [1:1:10,    12:2:20, 25:5:50];
    end
    sim = fqueu_NB(sim);
end
if 0% compare with fixed
    clear all
    %BER_vect = [0 1e-8:1e-8:1e-7  2e-7:1e-7:1e-6  2e-6:1e-6:1e-5  2e-5:1e-5:1e-4  1e-4:1e-4:1e-3  2e-3:1e-3:1e-2   2e-2:1e-2:1e-1 2e-1:1e-1:5e-1];
    
    if 1
         BER_vect = [0   2e-3:1e-3:1e-2  ];
      %  BER_vect = [ 2e-2:1e-2:1e-1 2e-1:1e-1:5e-1];
      %   BER_vect = [0 1e-8:1e-8:1e-7  1.1e-7:1e-8:1e-6  1.1e-6:1e-7:1e-5  1.1e-5:1e-6:1e-4  1.1e-4:1e-5:1e-3  1.1e-3:1e-4:1e-2   1.1e-2:1e-3:1e-1 1.1e-1:1e-2:5e-1];
        Ed_vect = nan(1, length(BER_vect)); K_vect = nan(1, length(BER_vect));
        for i = 1: length(BER_vect)
            BER=BER_vect(i);
            sim = init_sim(2);
            sim.trafficmode = 2;
            sim.SU =1;
            sim.CRN.Approx = 1; 
            sim.H = 128;
            sim.N = 8 ;
            sim.lambda = .1;
            sim.control.simulation=0;  sim.control.save_active=0;
            sim.CodedSystem=0;  sim.cognitive=1;
            sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
            sim.run.NUM_RUNS = 4;
            sim.Kv = [1:1:30];
            sim.AKv_step = 1; %soft graph
            sim.ch.Rch = 100;
            sim.ch.BER=BER;
            sim = fqueu_NB(sim);
            
            Ed = sim.Res.A_CRN1.Ed; KvA=sim.KvA;  %find optimal K
            if any(isfinite(Ed))
                [~, ind] = min(Ed);
                if ~isempty(ind),
                    Ed_vect(i)= Ed(ind); K_vect(i) = KvA(ind);
                end
            end
        end
        
        [BER_vect; K_vect]
        sim.control.plot_active = 1;
        sim.control.debug_active = 0;
        sim.control.runtime_plot_active=1;
        sim.control.simulation = 1;
        sim.run.NUM_RUNS = 4;
        %run.quick_run_coef = 0.1;
        K_vect(isnan(K_vect)) = [];
        for j = 1: length(BER_vect)
            sim.Kv = K_vect(j);
            BER=BER_vect(j);
            sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
            sim.ch.BER=BER;
            sim = fqueu_NB(sim);
            optimumdelay(j) = mean(sim.Res.S_CRN.Davsym);
        end
        
        K =  [K_vect,1,5,10,15];
        K = unique(K);
        legend_text = cell(1,length(K));
        for i = 1: length(K)
            sim.Kv = K(i);
            legend_text{i} = sprintf('k=%1.0d',K(i));
            for j = 1: length(BER_vect)BER=BER_vect(j);
                sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
                sim.ch.BER=BER;
                sim = fqueu_NB(sim);
                delay(i,j) = mean(sim.Res.S_CRN.Davsym);
            end
        end
        for j = 1: length(K_vect)
            avg_delay(j) =mean(delay(:,j));
        end
        address = ['\PCodeN\NData\BER_K\data_ber_k.mat'];
        save ([cd,address]);
        
        for i = 1: length(K)
            semilogy(BER_vect,delay(((length(K)+1)-i),:));
            hold on;
            
        end

        %BER_vect(10)=[];
        semilogy(BER_vect,avg_delay);
        semilogy(BER_vect,optimumdelay);
        
        legend_text = [legend_text,'Average','optimized'];
        legend(legend_text);
        grid on;
    end
    
    %     load ([cd,'\NData\BER_K\data_ber_k.mat']);
    %     %stem(BER_vect,K_vect)
    %     L = min(find(isnan(K_vect))); BER_vect=BER_vect(1:L); K_vect= K_vect(1:L);
    %     [K,ind]= unique(K_vect); B=BER_vect(ind);
    %     disp('limits for BER to change K: '); [B; K], pause;
    %     %Bmid = 0.5 * ([B, 0] + [0, B]); Bmid=Bmid(1:end-1);
    %
    %
    %
    %     xlabel('Channel Bit Error Probability \beta'); ylabel('Optimal Framing Parameter: k');
end
if 0% delay BER
    for BER = [0   2e-3:1e-3:1e-2  ];
        sim = init_sim(2);
        sim.CodedSystem=0;  sim.cognitive=1;
        sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
        sim.run.NUM_RUNS = 4;
        sim.Kv = [1:1:30];
        %%%   sim.AKv_step = 0.01; %soft graph
        sim.ch.Rch = 100; sim.ch.BER=1e-5; sim.control.simulation=1;
        sim.ch.BER=BER;
        sim = fqueu_NB(sim);
    end
end


if 0% BER-optK
    BER_vect =  [0   2e-3:1e-3:1.5e-2  ];
    %  BER_vect = [0 1e-8:1e-8:1e-7  2e-7:1e-7:1e-6  2e-6:1e-6:1e-5  2e-5:1e-5:1e-4  1e-4:1e-4:1e-3  2e-3:1e-3:1e-2   2e-2:1e-2:1e-1 2e-1:1e-1:5e-1];
    % BER_vect = [0   1e-4:1e-4:1e-3  2e-3:1e-3:1e-2   2e-2:1e-2:1e-1 ];
    if 1
        %     BER_vect = [0 1e-8:1e-8:1e-7  1.1e-7:1e-8:1e-6  1.1e-6:1e-7:1e-5  1.1e-5:1e-6:1e-4  1.1e-4:1e-5:1e-3  1.1e-3:1e-4:1e-2   1.1e-2:1e-3:1e-1 1.1e-1:1e-2:5e-1];
        Ed_vect = nan(1, length(BER_vect)); K_vect = nan(1, length(BER_vect));
        for i = 1: length(BER_vect)
            BER=BER_vect(i);
            sim = init_sim(2);
            sim.control.simulation=0;  sim.control.save_active=0;
            sim.CodedSystem=0;  sim.cognitive=1;
            sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
            sim.run.NUM_RUNS = 4;
            sim.Kv = [1:1:30];
            %sim.AKv_step = 1; %soft graph
            sim.ch.Rch = 100; sim.ch.BER=1e-5;
            sim.ch.BER=BER;
            sim = fqueu_NB(sim);
            
            Ed = sim.Res.A_CRN1.Ed; KvA=sim.KvA;  %find optimal K
            if any(isfinite(Ed))
                [~, ind] = min(Ed);
                if ~isempty(ind),
                    Ed_vect(i)= Ed(ind); K_vect(i) = KvA(ind);
                end
            end
        end
        [BER_vect; K_vect]
        save ([cd,'\NData\BER_K\data_ber_k.mat']);
    end
    
    %load ([cd,'\NData\BER_K\data_ber_k.mat']);
    %stem(BER_vect,K_vect)
    L = min(find(isnan(K_vect))); BER_vect=BER_vect(1:L); K_vect= K_vect(1:L);
    [K,ind]= unique(K_vect); B=BER_vect(ind);
    disp('limits for BER to change K: '); [B; K], pause;
    %Bmid = 0.5 * ([B, 0] + [0, B]); Bmid=Bmid(1:end-1);
    
    
    %bar(Bmid,K)
    figure; plotdef(0)
    for i=1:length(K)-1
        l=B(i); u=B(i+1); k=K(i);
        area([l, u], [k, k]); hold on;
    end
    xlabel('Channel Bit Error Probability \beta'); ylabel('Optimal Framing Parameter: k');
end




%Coded basic system
if 0
    for i=1:2
        sim = init_sim(2);
        sim.cognitive=1;
        
        if i==1
            sim.CodedSystem=1;  sim.coding.impBERD = 10; sim.coding.impBERH = 10;  sim.coding.CRD = 1/2; sim.coding.CRH=1/2; sim.coding.nTbitsD =3; sim.coding.nTbitsH =3;
            sim.ch.BER=1e-5;  sim.ch.Rch = 100;   sim.ch.Rch0 = 5;
        else
            sim.CodedSystem=1;  sim.coding.impBERD = 10; sim.coding.impBERH = 100;  sim.coding.CRD = 1/2; sim.coding.CRH=1/3; sim.coding.nTbitsD =3; sim.coding.nTbitsH =3;
            sim.ch.BER=1e-5;  sim.ch.Rch = 75;   sim.ch.Rch0 = 5;
        end
        
        
        if quick_run
            sim.control.plot_active = 1;
            sim.control.debug_active = 0;
            sim.control.runtime_plot_active=1;
            
            sim.run.NUM_RUNS = 4;
            sim.run.quick_run_coef = 0.1;
            sim.Kv = [1:1:10,    12:5:50];
        else
            sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
            sim.run.NUM_RUNS = 4;
            sim.run.quick_run_coef = 1;  %only for very quick test
            sim.Kv = [1:1:10,    12:2:20, 25:5:50];
        end
        
        sim = fqueu_NB(sim);
        if ispc && quick_run, plot_data_fileNB([cd, '\NData\'],sim.fname, 0); end
    end
end




if 0
    clear all
    chNotAvRatiov = [0.1 1 5 10 100];    %run for Ed
  % chNotAvRatiov = [0.01:0.01:0.1 0.2:0.1:1 2:1:100];   %run for ES vs u/v
   chNotAvRatiov = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20 50, 100];   %sparse run for ES vs u/v
     
    Rchv = 1000000 * ones(1,length(chNotAvRatiov));
 
    if 1
        for i = 1: length(chNotAvRatiov)
            sim = init_sim(2);
            sim.CodedSystem=0;  sim.coding.impBERD = 10; sim.coding.impBERH = 10; sim.coding.CTD = 1; sim.coding.CTH=1; sim.coding.CRD = 1/2; sim.coding.CRH=1/2; sim.coding.nTbitsD =3; sim.coding.nTbitsH =3;
            sim.control.simulation=0;
                 sim.trafficmode = 1;
            sim.datapath = regexprep([cd , '/NData/'], '\', '/'); addpath(sim.datapath);
            sim.run.NUM_RUNS = 8;
            sim.run.quick_run_coef = 1;
            sim.Kv = [1:1:50];
            sim.AKv_step = 0.01; %soft graph
%             if(sim.trafficmode == 1) %heavy
%                 sim.ch.Rch = 1000;
%                 sim.lambda = 5;
%             else
%                 sim.ch.Rch = 1000;
%                 sim.lambda = 1;
%             end
            sim.control.plot_active = 0;
            sim.control.save_active = 1;
             
            sim.ch.BER=1e-5;  sim.ch.Rch = Rchv(i);
             
            if length(chNotAvRatiov)>10,  %run for ES vs u/v
                sim.control.save_active=0;
                sim.AKv_step = 1;
            end
             
            sim.cognitive=1;
            sim.CRN.chNotAvRatio = chNotAvRatiov(i);
            sim.CRN.chAvPlenRatio = 0; sim.CRN.Ch_AvailMeanTime = 0.01;  %not variable Ch_AvailMeanTime
             
            %chNotAvRatiov(i), a=sim.CRN.chNotAvRatio, disp('param in call'); pause;
            sim = fqueu_NB(sim);
             
            ES(i,:) = sim.Res.A_CRN1.Es;
            ED(i,:) = sim.Res.A_CRN1.Ed;
        end
        save NData\CRN_ES\CRD_data
    end
     
end