function [sim] = fqueu_TB(sim)

%This program gets simulation parameters and simulates the queue
% input:  struct sim
% output: struct sim
% Procedure steps  
% step 1: find all the simulation parameters from the given simulation parameters [calc_sim_params]
% step 2: check if the simulation is already performed 
% step 3: run analysis to find expected average delay using kingman formula [fanalyzesim]
% step 4: rum simulations [framing and queue evolution] for different packetization interval Tv [f_perform_framing, f_q_evolution]

sim = calc_sim_params(sim);
read_sim

a = 1/lambda;            %lambda: Symbol Rate, a: inter-arrival time

t=clock;
fprintf('\n\n***********************************************************************************\n');
fprintf(' Date:%d-%d-%d   Time: %d:%d \n', floor(t(1:5)));
fprintf('***********************************************************************************\n');
fprintf('Parameters: N:%d  H:%d   Rate:%1.3f  Mean Inbterval Time:%1.3f  Packetization (Time):%1.3f \n', sim.N,sim.H,sim.lambda,1/sim.lambda, sim.T_av);
%fprintf('Framing Mode:%s     Packetization (Time):%1.3f   Packetization (K):%d\n', FramingModeStr{sim.Framing_mode}, sim.T_TB, sim.NinP_PB);
fprintf('Average Packet Length:%d  Info rate:%1.3f  BER:%d Packet error rate:%d  \n', ...
sim.Plen_av, lambda*N, sim.ch.BER, sim.ch.Exp_PER); 
fprintf('Service Rate:%1.3f  Mean Service time:%1.3f  Channel Bit rate:%1.3f\n', ...
sim.ch.Rch/sim.Plen_av, sim.Plen_av/Rch,  Rch);
t = clock;


if isempty(sim.fname) && sim.CodedSystem
    fname=sprintf('TB_L%d_N%d_H%d_BERU%1.3f_CS%d%d_CR%d%d_nT%d%d_impBER%d%d.mat', lambda, N,H, sim.ch.BER*10000, sim.coding.CTD, sim.coding.CTH, sim.coding.CRD, sim.coding.CRH, sim.coding.nTbitsD, sim.coding.nTbitsH, sim.coding.impBERD, sim.coding.impBERH);  
elseif isempty(sim.fname) && ~sim.CodedSystem
    fname=sprintf('TB_L%d_N%d_H%d_BER%1.3f.mat', lambda, N,H,sim.ch.BER*10000);  %fname=sprintf('TB_L%d_N%d_H%d_PER%d.mat', lambda, N,H,sim.ch.Exp_PER*1000);
else
    fname = sim.fname;
end
dfname = [sim.datapath, fname];

if ~strcmp(fname, 'test.mat') && ~isempty(dir(dfname)) && ~sim.control.overwrite_active
    fprintf('Simulation is already done. \n');
else
    fprintf('Simulation started, Date: %d-%d-%d  Time: %d:%d:%d\n', floor(t(1:6)));
    
    %********************************************************
    %  Analytical Calculation    
    %********************************************************
    if sim.control.analysis || sim.control.simulation
        sim = fanalyzesim(sim);
    end
    
    
    %********************************************************
    %  Simulations to calc Expected Waiting Time
    %********************************************************
    if sim.control.simulation
        fprintf('Calculate evolution of waiting time for time based using packet length and retransmission distribution.\n');
        fprintf('Time Based Method.\n');
        
        NLP = length(sim.Last_Percentage);  %to be more accurate, the average delay is calculated for the last packets [e.g 25%, 50%, 75%, ...] to avoid tranition effect in intitial queue development
        LTV = length(Tv);

        %figure(1);

        %intitialization
        z1=zeros(1,LTV); z2= zeros(NLP,LTV); z3=zeros(NLP,LTV,sim.run.NUM_RUNS);
        
        
        sim.nIntervals=z1; sim.nSymbols=z1; sim.nPackets=z1; sim.Tav=z1; sim.Fisum = z1;
        sim.wav = z2; sim.sav = z2; sim.rav = z2; sim.w2av = z2; sim.s2av = z2; sim.r2av = z2;
        wav_per_run = z3; sav_per_run = z3; rav_per_run = z3; w2av_per_run = z3; s2av_per_run = z3; r2av_per_run = z3;


        for Tind = 1: LTV
            T = Tv(Tind);   % run for packetization interval
            fprintf('\n Running for T:%1.3f  (%d/%d)\n', T, Tind, LTV);

            sim.Fisum(Tind) = 0;
            z=zeros(1, sim.run.NUM_RUNS); nIntervals=z; nSymbols=z; nPackets=z;  Tav=z;
            
            for nruns = 1:sim.run.NUM_RUNS   %repeat the simulation for NUM_RUN times for more accurate results 
                fprintf('.');
                
                %xlen defines the number of symbols for each run
                % it is lower bounded by minSym to avoid few symbols and low accuracy [effective for small T]
                % and is upper bounded by maxPacks to avoid extremely large           [effective for large T]

                %xlen = max(sim.run.minSym, min(sim.run.maxSym,  ceil(sim.run.maxPack*max(1,T*lambda))));
                xlen = max(sim.run.minSym, min(sim.run.maxSym, ceil(sim.run.maxPack*((T*lambda)/(1-exp(-T*lambda)))))); %is correct , later check
                
                %switch to this one after carefully check
                sim.run.nsyms(Tind)=xlen; %number of symbols
                
                
                
                %calculate time of packets [tp]  and number of symbols in each packet [kp], note kp can be zero [handled later]
                if (UFCPP)  
                    if (input_process == 1)              %Poisson process with rate 1/a, 
                        x = random('exp', a,[1,xlen]);   %generate interarrival times
                    elseif (input_process == 0)          %Deterministic inputs
                        x = a * ones (1,xlen);           %interarrival times are constant
                    end
                    [tp, kp, fiv] = f_perform_framing(Framing_mode, x, T, sim.control.debug_active); %lambda, 
                    sim.Fisum(Tind) = sim.Fisum(Tind) + mean(fiv);
                else  % use theoretical approach                
                    if (input_process == 1)       %Poisson process with rate 1/a,  
                        ninterval = ceil(maxPack * max(1, 1/ (1-exp(-T*lambda))));
                        tp = T * (1:1: ninterval);
                        kp = poissrnd(lambda * T, 1, ninterval);
                    else
                        error('Err: Not Developed for non-Poisson input');
                    end
                end

                nIntervals(nruns) = length(tp);
                nSymbols(nruns) = sum(kp);
                
                %delete zero length packets
                nonzero_ind = find(kp ~= 0);
                tp = tp(nonzero_ind);
                kp = kp(nonzero_ind);

                nPackets(nruns) =length(tp);
                Tav(nruns) = mean(tp(2:end)-tp(1:end-1));
                
                
                [wnv, sv, rv] = f_q_evolution(tp, kp, sim, Tind); %was [wnv, sv, rv] = f_q_evolution(tp, kp, N, H, BER, Rch, send_dummy_for_k0);
                

                %calculate statistics at lp last percentage of vectors
                for lpi = 1: NLP
                    lp = sim.Last_Percentage(lpi);
                    length(wnv);
                    stp = 1 + floor(((100-lp)/100) * length(wnv)); %start point

                    wav_per_run(lpi, Tind, nruns) = mean(wnv(stp:1:end));   
                    sav_per_run(lpi, Tind, nruns) = mean(sv(stp:1:end));
                    rav_per_run(lpi, Tind, nruns) = mean(rv(stp:1:end));

                    w2av_per_run(lpi, Tind, nruns) = mean(wnv(stp:1:end).^2);   
                    s2av_per_run(lpi, Tind, nruns) = mean(sv(stp:1:end).^2);
                    r2av_per_run(lpi, Tind, nruns) = mean(rv(stp:1:end).^2);


                end        
            end
            
            %average over runs
            sim.nIntervals(Tind) = mean(nIntervals);
            sim.nSymbols(Tind)=mean(nSymbols);
            sim.nPackets(Tind)=mean(nPackets);
            sim.Tav(Tind)=mean(Tav);

            for lpi = 1: NLP
                sim.wav(lpi, Tind) = mean(wav_per_run(lpi, Tind,:));
                sim.sav(lpi, Tind) = mean(sav_per_run(lpi, Tind,:)); 
                sim.rav(lpi, Tind) = mean(rav_per_run(lpi, Tind,:)); 

                sim.w2av(lpi, Tind) = mean(w2av_per_run(lpi, Tind,:));
                sim.s2av(lpi, Tind) = mean(s2av_per_run(lpi, Tind,:)); 
                sim.r2av(lpi, Tind) = mean(r2av_per_run(lpi, Tind,:)); 

            end

        end

        sim.Dav = sim.wav + sim.sav;  % average delay

        if (input_process == 1) && (UFCPP)
            sim.Davsym = sim.Dav + ones(NLP,1) * sim.Fisum / sim.run.NUM_RUNS;
        else
            sim.Davsym = sim.Dav + ones(NLP,1) * Tv/2;
        end
    end

    if sim.control.save_active
        t = clock;
        fprintf('data saved to: %s\n', dfname);
        save (dfname)
    end
    if sim.control.plot_wnv_active
        plot_data_fileTB(fname, 0);
    end
end

