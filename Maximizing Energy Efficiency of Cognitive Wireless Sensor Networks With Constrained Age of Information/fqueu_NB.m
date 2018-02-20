function [sim] = fqueu_NB(sim)
%This program gets simulation parameters and simulates the queue
% input:  struct sim
% output: struct sim
% Procedure steps  
% step 1: find all the simulation parameters from the given simulation parameters [calc_sim_params]
% step 2: check if the simulation is already performed 
% step 3: run analysis to find expected average delay using kingman formula [fanalyzesim]
% step 4: rum simulations [framing and queue evolution] for different packetization interval Kv [f_perform_framing, f_q_evolution]
    
%sim.CRN.Ch_AvailMeanTime , disp('calc_param'); pause;
sim = calc_sim_params(sim);

H = sim.H;  N = sim.N; lambda = sim.lambda;  a=1/lambda; Kv=sim.Kv;
%BER = sim.ch.BER; Rch = sim.ch.Rch; Plenv = sim.Plenv; PERv = sim.PERv; 
NUM_RUNS=sim.run.NUM_RUNS;

t=clock;
fprintf('\n\n***********************************************************************************\n');
fprintf(' Date:%d-%d-%d   Time: %d:%d \n', floor(t(1:5)));
fprintf('***********************************************************************************\n');
    
[sim.fname] = calc_fname(sim);
dfname = [sim.datapath, sim.fname];

fprintf('Simulation (Number Based) File:%s \n Params [sample rate:%f Info rate:%1.3f]  [N:%d] [H:%d]  [Coding Rate D:%1.3f  %1.3f] [Tail Bits D:%1.3f H:%1.3f] [BER:%d impN:%d  impH:%d]  [Channel Bit rate:%1.3f]\n', ...
    sim.fname, lambda, lambda*N, N, H, sim.coding.CRD, sim.coding.CRH, sim.coding.nTbitsD, sim.coding.nTbitsH, sim.ch.BER, sim.coding.impBERD, sim.coding.impBERH, sim.ch.Rch);
if sim.cognitive
    fprintf(' Cognitive Params chNotAvRatio:%f  chAvPlenRatio:%f   Ch_AvailMeanTime:%f   Ch_NotAvailMeanTime:%f\n', sim.CRN.chNotAvRatio, sim.CRN.chAvPlenRatio, sim.CRN.Ch_AvailMeanTime, sim.CRN.Ch_NotAvailMeanTime); end

SimulationAlreadyPerformed = 0;
if ~strcmp(sim.fname, 'test.mat') && ~isempty(dir(dfname)) && ~sim.control.overwrite_active
    load(dfname, 'sim');
    if isfield(sim, 'Res')
        if isfield(sim.Res, 'S_SN')
            if isfield(sim.Res.S_SN, 'Davsym')
                SimulationAlreadyPerformed = 1; fprintf('Simulation is already performed File:%s. \n', dfname);
            end
        end
    end
end

if ~ SimulationAlreadyPerformed    
    fprintf('Simulation started, Date: %d-%d-%d  Time: %d:%d:%d File:%s\n', floor(t(1:6)), sim.fname);
    %********************************************************
    %  Analytical Calculation    
    %********************************************************
    if sim.control.analysis || sim.control.simulation
        sim = fanalyzesimNB(sim);
        if 0 && sim.control.debug_active, save ('temp_analysis', 'sim'); disp('temp_analysis  saved'); pause; end
        if sim.control.save_active, save (dfname, 'sim'); end
    end

    
    %********************************************************
    %  Simulations to calc Expected Waiting Time
    %********************************************************
    if sim.control.simulation
        fprintf('Simulation of waiting time for time based using packet length and retransmission distribution.\n');
    
        nLP = length(sim.run.Last_Percentage);  %to be more accurate, the average delay is calculated for the last packets [e.g 25%, 50%, 75%, ...] to avoid tranition effect in intitial queue development
        nKv = length(Kv);
    nSU = length(sim.SU)
        %intitialization
        z1=zeros(1,nKv); z2= zeros(nLP,nKv); z3=zeros(nLP,nKv,NUM_RUNS); z4= zeros(nLP,nKv,nSU,NUM_RUNS);
        
        CRN=[]; SN=[];  %result struct for single network and cognitive radio
        
        %sim.nIntervals=z1; sim.nSymbols=z1; sim.nPackets=z1;  %SN.Tav=z1; 
        
        SN.fav = z1; SN.wav = z2; SN.sav = z2; SN.rav = z2; %sim.w2av = z2; sim.s2av = z2; sim.r2av = z2;
        wav_per_run = z4; sav_per_run = z4; rav_per_run = z4; w2av_per_run = z4; s2av_per_run = z4; r2av_per_run = z4;
        CRNwav_per_su = z4; CRNsav_per_su = z4; CRNrav_per_su = z4; CRNw2av_per_su = z4; CRNs2av_per_su = z4; CRNr2av_per_su = z4;
        
        CRN.fav = z1; CRN.wav = z2; CRN.sav = z2; CRN.rav = z2; %sim.w2av = z2; sim.s2av = z2; sim.r2av = z2;
        CRNwav_per_run = z3; CRNsav_per_run = z3; CRNw2av_per_run = z3; CRNs2av_per_run = z3; 
counter = 0;
%for Kind = [18,20,21,22,30]
for SU = sim.SU 
counter = counter+1;
        for Kind = 1: nKv 
            K = Kv(Kind);   % run for packetization interval
            ts=clock; 
            fprintf('\n Running for K:%1.3f  (%d/%d) Time:%d-%d-%d  %d:%d\n', K, Kind, nKv, floor(ts(1:5)));

            
            
            nsyms001= sim.run.Nsym; %nsyms defines the number of symbols for each run
            nsyms002= min(nsyms001, min(sim.run.maxSym, sim.run.maxPack*K)); 
           % nsyms= max(nsyms002, max(sim.run.minSym, sim.run.minPack*K)); 
            nsyms =  nsyms001;
            %repeat the simulation for NUM_RUN times for more accurate results 
            z=zeros(1, NUM_RUNS); Tav=z; T2av=z;  nSymbols=z; nPackets=z;  Fi_sum= 0; %Kav=z;  
            
            for nruns = 1:NUM_RUNS   %parfor
                fprintf('.');
                
                %calculate time of packets [tp]  and number of symbols in each packet [kp], note kp can be zero [handled later]
                if (sim.UFCPP)  
                    if (sim.input_process == 1)              %Poisson process with rate 1/a, 
                        int_arrivals = random('exp', a,[SU,nsyms]);   %generate interarrival times
                    elseif (sim.input_process == 0)          %Deterministic inputs
                        int_arrivals = a * ones (SU,nsyms);           %interarrival times are constant
                    end
                    [tp, kp, fiv, tindex] = f_perform_framing(sim.Framing_mode, int_arrivals, K, sim.control.debug_active); %lambda
                    Fi_sum= Fi_sum + mean(mean(fiv));
                else  % use theoretical approach                
                    if (sim.input_process == 1)       %Poisson process with rate 1/a,  
                        pack_int_arrivals = gamrnd(K, 1/lambda, [1,nsyms]); %shape =K, scale=1/lambda=a rate=lambda
                        tp(1)=pack_int_arrivals(1); for ii=2:length(pack_int_arrivals), tp(ii)=tp(ii-1)+pack_int_arrivals(ii); end
                        kp = K * ones(1, nsyms);
                    else
                        error('Err: Not Developed for non-Poisson input');
                     end
                end

                
% tp, pause;
% tp(2:end)-tp(1:end-1), pause;
% mean(tp(2:end)-tp(1:end-1)), pause;

                nSymbols(nruns) = sum(kp); nPackets(nruns) =length(kp);
                Tav(nruns) = mean(tp(2:end)-tp(1:end-1));       %average inter-packet times
                T2av(nruns) = mean((tp(2:end)-tp(1:end-1)).^2); %average inter-packet times
                
                queue.tp=tp; queue.kp=kp;
                [queue] = f_q_evolution(queue, sim, Kind,SU);% N, H, BER, Rch, send_dummy_for_k0);
%save temp2; disp('temp2 save'); pause;

                if 0 && sim.control.runtime_plot_active
                    plot_queue(queue, Kind)                     
                end

                %calculate statistics at lp last percentage of vectors
                for SUn = 1: SU
                    tmp.CRNwnv = queue.CRNwnv(find((tindex>(SUn-1)*length(kp))&(tindex<=SUn*length(kp))));
                    tmp.CRNsv = queue.CRNsv(find((tindex>(SUn-1)*length(kp))&(tindex<=SUn*length(kp))));
                    
                for lpi = 1: nLP
                    lp = sim.run.Last_Percentage(lpi);
                    stp = 1 + floor(((100-lp)/100) * length(queue.wnv)/SU); %start point

                    wav_per_run(lpi, Kind,SUn, nruns) = mean(queue.wnv(stp:end));   
                    sav_per_run(lpi, Kind,SUn, nruns) = mean(queue.sv(stp:end));
                    rav_per_run(lpi, Kind,SUn, nruns) = mean(queue.rv(stp:end));

                    w2av_per_run(lpi, Kind,SUn, nruns) = mean(queue.wnv(stp:1:end).^2);   
                    s2av_per_run(lpi, Kind,SUn, nruns) = mean(queue.sv(stp:1:end).^2);
                    r2av_per_run(lpi, Kind,SUn, nruns) = mean(queue.rv(stp:1:end).^2);
                    
                    if sim.cognitive
                        CRNwav_per_su(lpi, Kind,SUn, nruns) = mean(tmp.CRNwnv(stp:end));   
                        CRNsav_per_su(lpi, Kind,SUn, nruns) = mean(tmp.CRNsv(stp:end));

                        CRNw2av_per_su(lpi, Kind,SUn, nruns) = mean(tmp.CRNwnv(stp:1:end).^2);   
                        CRNs2av_per_su(lpi, Kind,SUn, nruns) = mean(tmp.CRNsv(stp:1:end).^2);
                    end
                end   
               
            end %for each su
            for lpi = 1: nLP
             if sim.cognitive
                    CRNwav_per_run(lpi, Kind,nruns) = mean(CRNwav_per_su(lpi, Kind,:,nruns));
                    CRNsav_per_run(lpi, Kind,nruns) = mean(CRNsav_per_su(lpi, Kind,:,nruns)); 
                    CRNw2av_per_run(lpi, Kind,nruns) = mean(CRNw2av_per_su(lpi, Kind,:,nruns));
                    CRNs2av_per_run(lpi, Kind,nruns) = mean(CRNs2av_per_su(lpi, Kind,:,nruns)); 
             end
            end
            end  %for nruns = 1:NUM_RUNS  
                        for lpi = 1: nLP
                SN.wav(lpi, Kind) = mean(wav_per_run(lpi, Kind,:));
                SN.sav(lpi, Kind) = mean(sav_per_run(lpi, Kind,:)); 
                SN.rav(lpi, Kind) = mean(rav_per_run(lpi, Kind,:)); 

                SN.w2av(lpi, Kind) = mean(w2av_per_run(lpi, Kind,:));
                SN.s2av(lpi, Kind) = mean(s2av_per_run(lpi, Kind,:)); 
                SN.r2av(lpi, Kind) = mean(r2av_per_run(lpi, Kind,:)); 
                if sim.cognitive
                    CRN.wav(lpi, Kind) = mean(CRNwav_per_run(lpi, Kind,:));
                    CRN.sav(lpi, Kind) = mean(CRNsav_per_run(lpi, Kind,:)); 
                    CRN.w2av(lpi, Kind) = mean(CRNw2av_per_run(lpi, Kind,:));
                    CRN.s2av(lpi, Kind) = mean(CRNs2av_per_run(lpi, Kind,:)); 
                end

            end

            sim.run.nsyms(Kind)=nsyms;
            %average over runs
            SN.nSymbols(Kind)=mean(nSymbols);
            SN.nPackets(Kind)=mean(nPackets);
            SN.Tav(Kind)=mean(Tav);  SN.T2av(Kind)=mean(T2av); %average over runs

          
            if (sim.input_process == 1) && (sim.UFCPP)
                SN.fav(Kind) = Fi_sum/ NUM_RUNS;
            else
                SN.fav(Kind) = (K-1)/(2*sim.lambda);
            end
            

            
            SN.Dav = SN.wav + SN.sav;  % average packet delay
            SN.Davsym = repmat(SN.fav,[nLP,1]) + SN.wav + SN.sav;  % average sample delay
            sim.Res.S_SN = SN; A_SN = sim.Res.A_SN; 
        
            if sim.cognitive
                %Copy unchanged fileds from SN to CRN 
                CRN.fav = SN.fav; CRN.rav=SN.rav; CRN.r2av=SN.r2av;   
                CRN.nSymbols = SN.nSymbols;  CRN.nPackets=SN.nPackets; CRN.Tav=SN.Tav;  CRN.T2av=SN.T2av; 
                CRN.Dav = CRN.wav + CRN.sav;  % average packet delay
                CRN.Davsym = repmat(CRN.fav,[nLP,1]) + CRN.wav + CRN.sav;  % average sample delay
                sim.Res.S_CRN = CRN;
                A_CRN = sim.Res.A_CRN1; 
            end
            
            
            
            %SN.nSymbols(Kind)=mean(nSymbols);
            %SN.nPackets(Kind)=mean(nPackets);
            %sim.run.nsyms(Kind)
            
            
            te=clock; fprintf('Results MeanInterval[A:%f S:%f] Packet_formation[A:%f S:%f]  ServiceTime[A:%f S:%f]   WaitingTime[A:%f S:%f]  Delay[A:%f S:%f]    [Nsym:%d=>%d=>%d=>%d Npacket:%d]   Time:[%d-%d %d:%d:%d to %d-%d-%d]   \n', ...
                A_SN.ET(Kind),SN.Tav(Kind),        A_SN.Ef(Kind), SN.fav(Kind),    A_SN.Es(Kind),SN.sav(end,Kind),       A_SN.Ew(Kind),SN.wav(end,Kind),       A_SN.Ed(Kind),SN.Davsym(end,Kind),  nsyms001,nsyms002,nsyms, SN.nSymbols(Kind),SN.nPackets(Kind),   sim.run.nsyms(Kind), floor(ts(1:5)), floor(te(4:5)));    
            if sim.cognitive, fprintf('\t\t\t\t COGNITIVE Results ServiceTime[A:%f S:%f]   WaitingTime[A:%f S:%f]  Delay[A:%f S:%f] \n', ...
                A_CRN.Es(Kind), CRN.sav(end,Kind),   A_CRN.Ew(Kind), CRN.wav(end,Kind),    A_CRN.Ed(Kind), CRN.Davsym(end,Kind));    end
        
            if sim.control.save_active
                save (strcat(dfname,'su=',int2str(SU))); 
            end  %save per each K

        end
        if(exist('sim.Res.S_CRN'))
 MultiSU.CRN{counter}=  sim.Res.S_CRN;   
        end
end
        if sim.control.save_active
            save (strcat(dfname,'su=',int2str(sim.SU),'.mat'));
            save (dfname, 'sim'); 
            if(length(sim.SU)>1)
               sim.fname = {strcat(dfname,'su=',int2str(sim.SU),'.mat')};
            end
        end
        if ispc && sim.control.runtime_plot_active,     plot_data_fileNB(sim.datapath,sim.fname, 0,''); end

    end
end