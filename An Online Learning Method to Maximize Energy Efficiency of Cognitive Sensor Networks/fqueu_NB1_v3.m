function optimum_k = fqueu_NB1_v2(sim)
%This program gets simulation parameters and simulates the queue
% input:  struct sim
% output: struct sim
% Procedure steps
% step 1: find all the simulation parameters from the given simulation parameters [calc_sim_params]
% step 2: check if the simulation is already performed
% step 3: run analysis to find expected average delay using kingman formula [fanalyzesim]
% step 4: rum simulations [framing and queue evolution] for different packetization interval Kv [f_perform_framing, f_q_evolution]

%sim.CRN.Ch_AvailMeanTime , disp('calc_param'); pause;
%sim = calc_sim_params(sim);

%H = sim.H;  N = sim.N; lambda = sim.lambda;  a=1/lambda; Kv=sim.Kv;
%BER = sim.ch.BER; Rch = sim.ch.Rch; Plenv = sim.Plenv; PERv = sim.PERv;
NUM_RUNS=sim.run.NUM_RUNS;
lambda = sim.lambda;
t=clock;
fprintf('\n\n***********************************************************************************\n');
fprintf(' Start: Date:%d-%d-%d   Time: %d:%d \n', floor(t(1:5)));
fprintf('***********************************************************************************\n');

sim.Av_length =[];
sim.B_length =[];
alislot = [];

%********************************************************
%  Simulations to calc Energy
%********************************************************
if sim.control.simulation
    fprintf('Simulation of Energy\n');
    
    nSU = length(sim.SU);
    %intitialization
    
    %for Kind = [18,20,21,22,30]
    for SU = sim.SU
        
        
        nsyms =  sim.run.Nsym;
        %repeat the simulation for NUM_RUN times for more accurate results
        z=zeros(1, NUM_RUNS); Tav=z; T2av=z;  nSymbols=z; nPackets=z;  Fi_sum= 0; %Kav=z;
        t = [];
        for nruns = 1:NUM_RUNS   %parfor
            fprintf('.');
            t(1)=0;
            t_t(1)=0;
            for j=1:length(lambda)
                a=1/lambda(j);
                int_arrivals = random('exp', a,[SU,nsyms]);   %generate interarrival times
                int_arrivals = [int_arrivals,int_arrivals];
                for i = 1:nsyms-1,
                    t(i+1) = t(i) + int_arrivals(i+1); %calc arrival times from interarrival times
                end
                % t_t = [t_t,t];
                t_t =t;
                fprintf('symbol generation set %d done!\n',j);
                %          [tp, kp, fiv, tindex] = f_perform_framing(sim.Framing_mode, int_arrivals, K, sim.control.debug_active); %lambda
            end
            fprintf('symbol generation done!\n');
            % u_rnd =randi([sim.uminsim.umax],1,sim.Nexperts); %rand(1,sim.Nexperts)*10;
            %v_rnd = randi([sim.vmin,sim.vmax],1,sim.Nexperts);%rand(1,sim.Nexperts)*10;sim.v .*ones(1,sim.Nexperts);
            u_rnd =datasample([sim.umin:0.1:sim.umax],sim.Nexperts); %rand(1,sim.Nexperts)*10;
            v_rnd = datasample([sim.vmin:0.1:sim.vmax],sim.Nexperts);%rand(1,sim.Nexperts)*10;sim.v .*ones(1,sim.Nexperts);
            
            Berr = [1e-4:1e-5:1e-3];
            %  BER_rnd = datasample(Berr,sim.Nexperts);
           % u_rnd(end) =1;
           % v_rnd(end) =1;
            %BER_rnd(end) =  sim.ch.BER;
            BER_rnd = ones(1,sim.Nexperts).*sim.ch.BER;
            %    BER_rnd = rand(1,sim.Nexperts);
            % k = [sim.k_min:sim.k_max];
            k = 11;    %Was ---------------------------------------------------> k = sim.k;
            %            k=[];
            %             for o=sim.k_min :sim.k_max
            %                 k = [k,ones(1,1e2)*o];
            %             end
            LastSymbol_flag=0;
            sim = channel_perform(sim);
            fprintf('Channel perform done!\n');
            ll=1;
            n=sim.n;
            t_t_copy = t_t;
            Nrho = 1;
            Kpredict = NaN;
            EE_best = Inf;
            KANAL = [];
            
            vtrack  = sim.vtrack;
            utrack  = sim.utrack;
            uvtrack = sim.slotuv;
            allv = zeros(1,length(uvtrack));
            allv(1:2:end) = vtrack;
            allv(2:2:end) = vtrack;
            KPRED = zeros(1,length(uvtrack));
            KPRED11 = [];
            k_pso = [];
            v_pso = [];
            u_pso =[];
            xxRTtrack = [];
            t_t_copy = t_t;
            sim_copy = sim;
            k_static = 19;
            k_static = ones(1,n).*k_static;
            segment = 1;
            
            while(~LastSymbol_flag)
                %simulate Energy  with fixed k,system real v ,BER
                k = ones(1,n).*k;
                [aa,NN1,NN2,NN3,EE_real,t_t,sim] = EE_simulator(k,sim,t_t,n);
                fprintf('simulation for iteration %d with dynamic k done!\n',ll);
                NumberofCollision(ll) = mean(NN1);
                %commented 09/13/2017    [dummy1,dummy2,dummy3,dummy4,EE_real_static,t_t_copy,sim_copy] = EE_simulator(k_static,sim_copy,t_t_copy,n);
                fprintf('simulation for iteration %d with static k done!\n',ll);
                fprintf('small:%f busy:%f Av:%f\n',mean(NN1),mean(NN2),mean(NN3));
                
                
                if(segment == 1)
                    EE_real_window = EE_real ;
                    n= sim.windowsshift;%round(sim.n/3);
                    
                    EE_real = mean(EE_real_window);
                    EE_real_track = EE_real;
                else
                    EE_real_window = [EE_real_window, EE_real];
                    
                    
                    for yu =1:sim.windowsshift
                        
                        EE_real_window(yu)=[];
                    end
                    EE_real = mean(EE_real_window);
                    EE_real_track = [EE_real_track EE_real];
                end
                
                %commented 09/13/2017 EE_real_static = mean(EE_real_static);
                %commented 09/13/2017 EE_real_static_track(ll) = EE_real_static;
                tic
                for ww=1:200
                    %  ww
                    %Analyticaly calculate Energy for random generated points
                    [EE_prediction_vector,P_s,perr] = EE_Analisor(k,1,sim,u_rnd,v_rnd,BER_rnd,0);
                    EE_real_vector = ones(1,sim.Nexperts).*EE_real;
                    P_s_track(ll,:) = P_s;
                    % find best match for Energy
                    diff  = abs(EE_real_vector-EE_prediction_vector);
                    diff_P_s= abs(NumberofCollision(end)-P_s);
                    diff_perr = abs(mean(aa)-perr);
                    [aps,bps]=min(diff_P_s);
                    [ape,bpe]=min(diff_perr);
                    P_s_diff_track(ww)= aps;
                    P_err_diff_track(ww)= ape;
                    P_s_copy = P_s;
                    perr_copy = perr;
                    % [a,b] = sort(diff);
                    %[a,b] = sort((diff_P_s./sum(diff_P_s)));%+(diff_perr./sum(diff_perr)));
                    [a,b] = sort((diff_P_s));
                    if(length(unique(a))==1)
                        ali=2;
                    end
                    c(ww,:)=b(1:5);
                    %                 v_rnd(b)
                    %                 u_rnd(b)
                    %                 BER_rnd(b)
                    [EE_diff_best,Best_prediction] = min(diff);
                    vpredict(ww) =  v_rnd(Best_prediction);
                    upredict(ww) =  u_rnd(Best_prediction);
                    BERpredict(ww) = BER_rnd(Best_prediction);
                    
                    %update random parameters
                    
                    [u_rnd,v_rnd,BER_rnd,sim] = expert_update(u_rnd,v_rnd,BER_rnd,Best_prediction,sim,EE_diff_best,diff);
                    
                    
                    
                    
                    
                    %find optimum k for best parameters which are chosen last
                    %part
                    k_t = [sim.k_min:sim.k_max];
                    [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,upredict(end),vpredict(end),BERpredict(end),1);
                    [xx,k_suggest]= min(EE_prediction_vector1);
                    kkkkk(ww) =  k_suggest;
                    k=k(1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%temporary
                    if(isnan(Kpredict))
                        Kpredict = k;
                        kfloat =k_suggest;
                    end
                    if(EE_real_track(end)<EE_best)
                        k_mbest = k(1);
                        EE_best = EE_real_track(end);
                    end
                    %   toc
                end
                k_pso = [k_pso ; kkkkk];
                v_pso = [v_pso;vpredict];
                u_pso = [u_pso;upredict];
                perr_copy(bpe)
               P_s_copy(bps)
                P_s_finaldiff_track(segment)= aps;
                P_err_finaldiff_track(segment)= ape; 
                P_s_track1(segment)= P_s_copy(bps);
                P_err_track(segment)= perr_copy(bpe);
                P_err_real_track(segment)= NumberofCollision(end);
                P_s_real_track1(segment)= mean(aa);
                fprintf('PSO done!\n');
                EE_diff_best_track(ll) = EE_diff_best;
                % u_rnd =randi([sim.umin,sim.umax],1,sim.Nexperts);
                % v_rnd = randi([sim.vmin,sim.vmax],1,sim.Nexperts);
                u_rnd =datasample([sim.umin:0.1:sim.umax],sim.Nexperts); %rand(1,sim.Nexperts)*10;
                v_rnd = datasample([sim.vmin:0.1:sim.vmax],sim.Nexperts);%rand(1,sim.Nexperts)*10;sim.v .*ones(1,sim.Nexperts);
                
             %   BER_rnd = datasample(Berr,sim.Nexperts);%Berr(randperm(length(Berr),sim.Nexperts));
                sim = value_reset( sim );
              %  u_rnd(end) = upredict(end);
               % v_rnd(end) =vpredict(end);
                % BER_rnd(end) = BERpredict(end);
                BER_rnd = ones(1,sim.Nexperts).*sim.ch.BER;
                                     u_rnd(end) =1;
                                      v_rnd(end) =1;
                %                     BER_rnd(end-1) =  sim.ch.BER;
                
                k_last = k;
                k = k_suggest; % was k  = round((sim.omega)*Kpredict(end)+(1-sim.omega)*k_suggest);
                %   kfloat = (1-sim.alpha).*kfloat+sim.beta.* (Kpredict(end)-k_mbest)+sim.gamma.*(Kpredict(end)-k_suggest);
                %                if(ll>1)
                %                 k= round((k+Kpredict(end))/2);
                %                end
                %                 if((abs(k-k_last)>3)&&length(Kpredict)>2)
                %                     k = k_last;
                %                 end
                if(k<1)
                    k=1;
                end
                if(k>sim.k_max)
                    k = sim.k_max;
                end
                Kpredict(ll) = k;
                %   xx(ll) = xx;
                length(t_t);
                
                fprintf('segment :%d u:%f v:%f BER:%f k:%f\n',segment,upredict(end),vpredict(end),BERpredict(end),k);
                segment = segment+1;
                %################################################################################################
                %                                  Real time simulation
                %################################################################################################
                if(exist('C:\Temp\matlab\alldata.mat', 'file') == 2)
                    delete 'C:\Temp\matlab\alldata.mat';
                end
                % save test -regexp ^(?!(t|int_arrivals|t_t|t_t_copy)$).
                save('C:\Temp\matlab\alldata.mat','-regexp', '^(?!(t|int_arrivals|t_t|t_t_copy)$).');
         %   if(segment>25)
%                 if(sim.variable == 0)%was0
%                     close all;
%                     plot(EE_real_track);
%                     hold on;
%                     plot(EE_real_track+EE_diff_best_track);
%                     title('best particle energy track');
%                     
%                     figure;
%                     uuu = sim.u;
%                     vvv = sim.v;
%                     BBB = sim.ch.BER;
%                     [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uuu,vvv,BBB,1);
%                     [xxRT,kAnalRT]= min(EE_prediction_vector1RT);
%                     k_last
%                     plot(Kpredict);
%                     hold on;
%                     plot(.05+(kAnalRT.*ones(1,length(Kpredict))));
%                     [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k,Best_prediction,sim,uuu,vvv,BBB,1);
%                     title('k track');
%                     figure;
%                     plot(EE_real_track);
%                     hold on;
%                     plot(.05+(EE_prediction_vector1RT*ones(1,length(EE_real_track))));
%                     title('real and simulation');
%                     figure;
%                     plot(EE_diff_best_track);
%                     title('Difference');
%                     pause(.001);
%                 else
%                     close all;
% %                     Kanallll(1) =0;
% %                     plot(100.*[1:length(EE_real_track)],EE_real_track);
% %                     hold on;
% %                     % plot(EE_real_static_track);
% %                     %  hold off
% %                     %legend('dynamic','static');
% %                     %title('static vs dynamic - dynamic v');
% %                     %  plot(EE_real_track+EE_diff_best_track);
% %                     % title('best particle energy track');
% %                     
% %                     
% %                     BBB = sim.ch.BER;
% %                     alislot22 = unique(uvtrack(sim.startslot:sim.endslot));
% %                     alislot = [alislot, nan, alislot22];
% %             
% %                     [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uvtrack(2),uvtrack(1),BBB,1); %temporary
% %                     for op=[sim.startslot:2:sim.endslot];
% %                         %try
% % %                             if(mod(sim.startslot,2)~=0)
% % %                               [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uvtrack(op+1),uvtrack(op),BBB,1);
% % %                             else
% % %                                 [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uvtrack(op),uvtrack(op+1),BBB,1);
% % %                             end
% %                             [xxRT,kAnalRT]= min(EE_prediction_vector1RT);
% %                             Kanallll(end+1) = kAnalRT;
% %                         %catch
% %                         %    ali=1;
% %                         end
% %                     end
% %                     xxRTtrack =   [xxRTtrack, xxRT];
% %                     plot(100.*[1:length(xxRTtrack)],xxRTtrack);
% %                     
% %                     %legend('dynamic','static');
% %                     %title('static vs dynamic - dynamic v');
% %                     plot(100.*[1:length(EE_real_track)],EE_real_track+EE_diff_best_track);
% %                     hold off
% %                     title('best particle energy track');
% %                     
% %                     Kanallll(1)=[];
% %                     figure;
% %                    % KANAL = [KANAL ,datasample(Kanallll,100,'Replace',false) ];
% %                      
%                 if(segment == 2)
% %                    KANAL = [KANAL ,datasample(Kanallll,sim.n,'Replace',false) ];
%                     KPRED11(1:sim.n) = .05+(k_last.*ones(1,sim.n));
%                 else
% %                    KANAL = [KANAL ,datasample(Kanallll,sim.windowsshift,'Replace',false) ];
%                     lastindex = length(KPRED11);
%                     for wq=1:sim.windowsshift
%                     KPRED11(lastindex+wq) = .05+(k_last);%.05+(k_last.*ones(1,sim.windowsshift));
%                     end
%                     
%                 end  
%                 
%                         kanal = sim.kanal;
%                     Kanalplot = kanal(sim.slottrack);
%                     kanal = [];
%                     %KPRED(sim.startslot:sim.endslot) = .05+(k_last.*ones(1,length(sim.startslot:sim.endslot)));
%                     %    KPRED = [KPRED , .05+(Kpredict(end).*ones(1,length(Kanallll)))];
%                     
%                     %Kanalplot = [Kanallll,zeros(1,length(KPRED)-length(KANAL)) ];
%                     plot(KPRED11);
%                     hold on;
%                     plot(Kanalplot);
%                     hold off
%                     title('k track dynamic v');
%                     
%                     %                     for o=2:length(KPRED1)
%                     %                     if(KPRED1(o)==0)
%                     %                         KPRED1(o) = KPRED1(o-1);
%                     %                     end
%                     %                     end
%                     %                     [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k,Best_prediction,sim,uuu,vvv,BBB,1);
%                     %                     figure;
%                     %                     plot(EE_real_track);
%                     %                     hold on;
%                     %                     plot(.05+(EE_prediction_vector1RT*ones(1,length(EE_real_track))));
%                     %                     title('real and simulation');
%                     %
% %                     figure;
% %                     plot(P_err_finaldiff_track);
% %                     title('Difference Perr');
% %                      figure;
% %                     plot(P_s_finaldiff_track);
% %                      title('Difference P_s');
% %                     figure;
% %                     plot(P_err_diff_track);
% %                     title('Difference PSO Perr');
% %                      figure;
% %                     plot(P_s_diff_track);
% %                      title('Difference PSO P_s');
%                      
%                        figure;
%                     plot(P_s_track1);
%                     hold on
%                     plot(P_s_real_track1);
%                     title('Difference Perr');
%                         figure;
%                     plot(P_err_track);
%                     hold on
%                     plot(P_err_real_track);
%                     title('Difference Perr');
%                  
% %                      figure;
% %                     plot(EE_diff_best_track);
% %                     title('Difference dynamic v');
%                     pause(.001);
%                     
%                 end
        %     end
                %################################################################################################
                %                                END------------Real time simulation
                %################################################################################################
                if(sim.channel_change_flag==1)
                    sim.channel_change_flag =0;
                    
                    %   LastSymbol_flag = 1;
                    uuu = sim.u;
                    vvv = sim.v;
                    BBB = sim.ch.BER;
                    [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,uuu(sim.Nparamset-1),vvv(sim.Nparamset-1),BBB(sim.Nparamset-1),1);
                    [xx,kAnal]= min(EE_prediction_vector1);
                    vars = {'Kpredict(end)','BERpredict(end)','upredict(end)','vpredict(end)','EE_real_track(end)','EE_diff_best_track(end)','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    
                    save(['result_for' num2str(Nrho) '.mat']);
                    Nrho = Nrho+1;
                    optimum_k = k;
                    fprintf('u:%f v:%f BER:%f kcalculated:%f kanal:%f\n',upredict(end),vpredict(end),BERpredict(end),k,kAnal);
                    % plot(EE_real_track)
                    % hold on;
                    %plot(EE_real_track+EE_diff_best_track)
                    %hold off, semilogy(EE_real_track), semilogy(EE_real_track+EE_diff_best_track);
                    % semilogy(EE_real_track),hold all ,semilogy(EE_real_track+EE_diff_best_track);
                    %figure
                    %plot(Kpredict),hold on, plot(ones(1,length(Kpredict))*kAnal)
                    % plot(EE_prediction_vector1)
                    ll = 1;
                end
                if((length(t_t)<k*n)||ll>sim.Itlimit)
                    
                    LastSymbol_flag = 1;
                    sim.channel_change_flag =0;
                    
                    %   LastSymbol_flag = 1;
                    uuu = sim.u;
                    vvv = sim.v;
                    BBB = sim.ch.BER;
                    [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,uuu,vvv,BBB,1);
                    [xx,kAnal]= min(EE_prediction_vector1);
                    if(length(Kpredict)>1)
                        perc = length(Kpredict);
                        number = ceil(.50*perc);
                        Kpredictend = mean(Kpredict(number:end));
                    else
                        Kpredictend = Kpredict(end);
                    end
                    BERpredictend = BERpredict(end);
                    upredictend = upredict(end);
                    vpredictend = vpredict(end);
                    EE_real_trackend = EE_real_track(end);
                    EE_diff_best_trackend = EE_diff_best_track(end);
                    %  k = mean();
                    
                    %  vars = {'Kpredict(end)','BERpredict(end)','upredict(end)','vpredict(end)','EE_real_track(end)','EE_diff_best_track(end)','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    vars = {'Kpredict','BERpredict','upredict','vpredictend','EE_real_trackend','EE_diff_best_trackend','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    save(['result_for' num2str(sim.i) '.mat'],vars{:});
                    Nrho = Nrho+1;
                    optimum_k = k;
                    fprintf('u:%f v:%f BER:%f kcalculated:%f kanal:%f\n',upredict(end),vpredict(end),BERpredict(end),round(mean(Kpredict)),kAnal);
                    
                    ll = 1;
                    %                     fprintf('u:%f v:%f BER:%f k:%f\n',upredict(end),vpredict(end),BERpredict(end),k);
                    %                     k = [sim.k_min:sim.k_max];
                    %                     EE_prediction_vector1 = EE_Analisor(k,Best_prediction,sim,sim.u,sim.v,sim.ch.BER,0);
                    %                     [xx,k]= min(EE_prediction_vector1)
                    %                     plot(EE_prediction_vector1)
                end
            end
            
            % pause;
            % tp, pause;
            % tp(2:end)-tp(1:end-1), pause;
            % mean(tp(2:end)-tp(1:end-1)), pause;
            % figure
            %  plot(Kpredict)
        end
    end
end