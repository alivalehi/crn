function optimum_k = fqueu_NB1(sim)
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
            u_rnd =randi([sim.umin,sim.umax],1,sim.Nexperts); %rand(1,sim.Nexperts)*10;
            v_rnd = randi([sim.vmin,sim.vmax],1,sim.Nexperts);%rand(1,sim.Nexperts)*10;sim.v .*ones(1,sim.Nexperts);
            Berr = [5e-5:1e-5:5e-3];
            BER_rnd = datasample(Berr,sim.Nexperts);
            %             u_rnd(end) =sim.u;
            %             v_rnd(end) =sim.v;
            %             BER_rnd(end) =  sim.ch.BER;
            %  BER_rnd = ones(1,sim.Nexperts).*sim.ch.BER;
            %    BER_rnd = rand(1,sim.Nexperts);
            % k = [sim.k_min:sim.k_max];
            k = 10;    %Was ---------------------------------------------------> k = sim.k;
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
            while(~LastSymbol_flag)
                %simulate Energy  with fixed k,system real v ,BER
                k = ones(1,n).*k;
                [aa,NN1,NN2,NN3,EE_real,t_t,sim] = EE_simulator(k,sim,t_t,n);
                fprintf('simulation for iteration %d done!\n',ll);
                NNN1(ll) = mean(NN1);
                fprintf('small:%f busy:%f Av:%f\n',mean(NN1),mean(NN2),mean(NN3));
                EE_real = mean(EE_real);
                EE_real_track(ll) = EE_real;
                %Analyticaly calculate Energy for random generated points
                [EE_prediction_vector,P_s,perr] = EE_Analisor(k,1,sim,u_rnd,v_rnd,BER_rnd,0);
                EE_real_vector = ones(1,sim.Nexperts).*EE_real;
                P_s_track(ll,:) = P_s;
                % find best match for Energy
                diff  = abs(EE_real_vector-EE_prediction_vector);
                [a,b] = sort(diff);
                if(length(unique(a))==1)
                    ali=2;
                end
                c(ll,:)=b;
                %                 v_rnd(b)
                %                 u_rnd(b)
                %                 BER_rnd(b)
                [EE_diff_best,Best_prediction] = min(diff);
                vpredict(ll) =  v_rnd(Best_prediction);
                upredict(ll) =  u_rnd(Best_prediction);
                BERpredict(ll) = BER_rnd(Best_prediction);
                EE_diff_best_track(ll) = EE_diff_best;
                %update random parameters
                %######################################### Added atJuly 13
                %instead of following line
                %  [u_rnd,v_rnd,BER_rnd,sim] = expert_update(u_rnd,v_rnd,BER_rnd,Best_prediction,sim,EE_diff_best,diff);
                %#########################################################
                u_rnd =randi([sim.umin,sim.umax],1,sim.Nexperts);
                v_rnd = randi([sim.vmin,sim.vmax],1,sim.Nexperts);
                BER_rnd = datasample(Berr,sim.Nexperts);%Berr(randperm(length(Berr),sim.Nexperts));
                u_rnd(end) = upredict(end);
                v_rnd(end) =vpredict(end);
                BER_rnd(end) = BERpredict(end);
                
                
                
                %find optimum k for best parameters which are chosed last
                %part
                k_t = [sim.k_min:sim.k_max];
                [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,upredict(end),vpredict(end),BERpredict(end),1);
                [xx,k_suggest]= min(EE_prediction_vector1);
                k=k(1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%temporary
                if(isnan(Kpredict))
                    Kpredict = k;
                    kfloat =k_suggest;
                end
                if(EE_real_track(end)<EE_best)
                    k_mbest = k(1);
                    EE_best = EE_real_track(end);
                end
                k_last = k;
                k = k_suggest; % was k  = round((sim.omega)*Kpredict(end)+(1-sim.omega)*k_suggest);
                %   kfloat = (1-sim.alpha).*kfloat+sim.beta.* (Kpredict(end)-k_mbest)+sim.gamma.*(Kpredict(end)-k_suggest);
                % k= round((k+Kpredict(end))/2);
                if((abs(k-k_last)>3)&&length(Kpredict)>2)
                    k = k_last;
                end
                if(k<1)
                    k=1;
                end
                if(k>sim.k_max)
                    k = sim.k_max;
                end
                Kpredict(ll) = k;
                %   xx(ll) = xx;
                length(t_t);
                
                fprintf('iteration :%d u:%f v:%f BER:%f k:%f\n',ll,upredict(end),vpredict(end),BERpredict(end),k);
                ll = ll+1;
                %################################################################################################
                %                                  Real time simulation
                %################################################################################################
                if(sim.variable == 0)
                    close all;
                    plot(EE_real_track);
                    hold on;
                    plot(EE_real_track+EE_diff_best_track);
                    title('best particle energy track');
                    
                    figure;
                    uuu = sim.u;
                    vvv = sim.v;
                    BBB = sim.ch.BER;
                    [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uuu,vvv,BBB,1);
                    [xxRT,kAnalRT]= min(EE_prediction_vector1RT);
                    k_last
                    plot(Kpredict);
                    hold on;
                    plot(.05+(kAnalRT*ones(1,length(Kpredict))));
                    [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k,Best_prediction,sim,uuu,vvv,BBB,1);
                    title('k track');
                    figure;
                    plot(EE_real_track);
                    hold on;
                    plot(.05+(EE_prediction_vector1RT*ones(1,length(EE_real_track))));
                    title('real and simulation');
                    figure;
                    plot(EE_diff_best_track);
                    title('Difference');
                    pause(.001);
                else
                    close all;
                    Kanallll(1) =0;
                    plot(EE_real_track);
                    hold on;
                    plot(EE_real_track+EE_diff_best_track);
                    title('best particle energy track');
                    
                    vtrack  = sim.vtrack;
                    utrack  = sim.utrack;
                    uvtrack = sim.slotuv;
                    BBB = sim.ch.BER;
                    alislot22 = unique(uvtrack(sim.startslot:sim.endslot))
                    alislot = [alislot, nan, alislot22];
                    for op=[sim.startslot:2:sim.endslot];
                        try
                            if(mod(sim.startslot,2)~=0)
                                [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uvtrack(op+1),uvtrack(op),BBB,1);
                            else
                                [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k_t,Best_prediction,sim,uvtrack(op),uvtrack(op+1),BBB,1);
                            end    
                            [xxRT,kAnalRT]= min(EE_prediction_vector1RT);
                            Kanallll(end+1) = kAnalRT;
                        catch
                            ali=1;
                        end
                    end
                    Kanallll(1)=[];
                    figure;
                    plot(Kpredict);
                    hold on;
                    plot(.05+(kAnalRT*ones(1,length(Kpredict))));
                    title('k track');
                    
                    
                    
                    %                     [EE_prediction_vector1RT,P_sRT,perrRT] = EE_Analisor(k,Best_prediction,sim,uuu,vvv,BBB,1);
                    %                     figure;
                    %                     plot(EE_real_track);
                    %                     hold on;
                    %                     plot(.05+(EE_prediction_vector1RT*ones(1,length(EE_real_track))));
                    %                     title('real and simulation');
                    %
                    figure;
                    plot(EE_diff_best_track);
                    title('Difference');
                    pause(.001);
                    
                end
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
            figure
            plot(Kpredict)
        end
    end
end