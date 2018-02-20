function optimum_k = fqueu_NB(sim)
%This program gets simulation parameters and simulates the queue
% input:  struct sim
% output: struct sim
% Procedure steps
% step 1:
% step 2:
% step 3:
% step 4:

%sim.CRN.Ch_AvailMeanTime , disp('calc_param'); pause;
%sim = calc_sim_params(sim);

%H = sim.H;  N = sim.N; lambda = sim.lambda;  a=1/lambda; Kv=sim.Kv;
%BER = sim.ch.BER; Rch = sim.ch.Rch; Plenv = sim.Plenv; PERv = sim.PERv;
NUM_RUNS=sim.run.NUM_RUNS;
lambda = sim.lambda;
t=clock;
fprintf('\n\n***********************************************************************************\n');
fprintf(' Date:%d-%d-%d   Time: %d:%d \n', floor(t(1:5)));
fprintf('***********************************************************************************\n');




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
        
        for nruns = 1:NUM_RUNS   %parfor
            %##############################################################################################################
            % Step 1: generating symbols
            %##############################################################################################################
            
            fprintf('.');
            t=0;
            t_t=0;
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
            %##############################################################################################################
            % Step 2: generating channel
            %##############################################################################################################
            u_rnd =randi([sim.umin,sim.umax],1,sim.Nexperts); %rand(1,sim.Nexperts)*10;
            v_rnd = randi([sim.vmin,sim.vmax],1,sim.Nexperts);%rand(1,sim.Nexperts)*10;sim.v .*ones(1,sim.Nexperts);
            Berr = [1e-5:1e-6:5e-2];%changed from Berr = [5e-5:1e-6:5e-4]; 7/5/2017
            BER_rnd = Berr(randperm(length(Berr),sim.Nexperts));
            %             u_rnd(end) =sim.u;
            %             v_rnd(end) =sim.v;
            %             BER_rnd(end) =  sim.ch.BER;
            
            BER_rnd = ones(1,sim.Nexperts).*sim.ch.BER;
            %    BER_rnd = rand(1,sim.Nexperts);
            % k = [sim.k_min:sim.k_max];
            k = sim.k;
            first_k = k;
            %            k=[];
            %             for o=sim.k_min :sim.k_max
            %                 k = [k,ones(1,1e2)*o];
            %             end
            LastSymbol_flag=0;
            sim = channel_perform(sim);
            fprintf('Channel perform done!\n');
            %##############################################################################################################
            % Step 3: simulation and main learning algorithm
            %##############################################################################################################
            lll=1;
            ll=1;
            %iterationCounter = 1;
            n=sim.n;
            t_t_copy = t_t;
            Nrho = 1;
            Kpredict = NaN;
            EE_best = Inf;
            %             if(sim.method == 2)
            %                 k = ones(1,n).*k;
            %                 [aa,NN1,NN2,NN3,EE_real,t_t,sim] = EE_simulator(k,sim,t_t,n);
            %             end
            
            while(~LastSymbol_flag)
                %k = randi([sim.k_min,sim.k_max],1); %just for test should be deleted
                %  if(sim.method == 1)
                k = ones(1,n).*k;
                [aa,NN1,NN2,NN3,EE_real,t_t,sim] = EE_simulator(k,sim,t_t,n);
                fprintf('simulation for iteration %d done!\n',lll);
               pause(0.0001);
                %   end
                NNN1(ll) = mean(NN1);
                fprintf('small:%f busy:%f Av:%f\n',mean(NN1),mean(NN2),mean(NN3));
                EE_real = mean(EE_real);
                EE_real_track(ll) = EE_real;

           
         
          for iterationCounter =1:sim.Itlimit
                    
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
                    [u_rnd,v_rnd,BER_rnd,sim] = expert_update(u_rnd,v_rnd,BER_rnd,Best_prediction,sim,EE_diff_best,diff);
                    %find optimum k for best parameters which are chosed last
                    %part
                    k_t = [sim.k_min:sim.k_max];
                    [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,upredict(end),vpredict(end),BERpredict(end),1);
                    [xx,k_suggest]= min(EE_prediction_vector1);
                    %Kpredict(ll) = k_suggest;
                    k=k(1);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%temporary
                    if(isnan(Kpredict))
                        Kpredict = k;
                        kfloat =k_suggest;
                    end
                    if(EE_real_track(end)<EE_best)
                        k_mbest = k(1);
                        EE_best = EE_real_track(end);
                    end
                    %k = k_suggest;
                    %k  = round((1-sim.omega)*Kpredict(end)+(sim.omega)*k_suggest);
                    ll = ll+1;
                    Kpredictl(ll) =  k_suggest;
                end
                unique(Kpredictl)
                Kpredictl
                %                 if(1)
                %                     if(ll<round(.6*sim.Itlimit))
                %                         k  = round((1-sim.omega)*mode(Kpredict)+(sim.omega)*k_suggest);
                %                         [EE_prediction_vector5,P_s,perr] = EE_Analisor(k,Best_prediction,sim,upredict(end),vpredict(end),BERpredict(end),1);
                %                         [EE_prediction_vector6,P_s,perr] = EE_Analisor(k+1,Best_prediction,sim,upredict(end),vpredict(end),BERpredict(end),1);
                %                         [EE_prediction_vector4,P_s,perr] = EE_Analisor(k-1,Best_prediction,sim,upredict(end),vpredict(end),BERpredict(end),1);
                %                         kres = [k-1 k k+1];
                %                         [eeres,keeres] = min([EE_prediction_vector4,EE_prediction_vector5,EE_prediction_vector6]);
                %                         k = kres(keeres);
                %                         %                     [m1,ff1] = mode(Kpredict);
                %                         %                     Kpredictcopy = Kpredict ;
                %                         %                     Kpredictcopy(Kpredict==m1) = [];
                %                         %                     [m2,ff2] = mode(Kpredictcopy);
                %                         %                     % [kk,ff] = mode(Kpredict);
                %                         %                     if((ff1-ff2)<=.2*ff1)
                %                         %                         kkk = .5*m1+.5*m2;
                %                         %                     else
                %                         %                         kkk= mode(Kpredict);
                %                         %                     end
                %                         %                     k  = round((1-sim.omega)*kkk+(sim.omega)*k_suggest);
                %                         % k1  = floor((1-sim.omega)*mode(Kpredict(ll-1:end))+(sim.omega)*k_suggest);
                %                         % k2  = ceil((1-sim.omega)*mode(Kpredict(ll-1:end))+(sim.omega)*k_suggest);
                %                         % k3 = [k1 k2];
                %                         % k = k3(randperm(2,1));
                %                         sim.omega = sim.omega*.99;
                %                     else
                %                         k  = round((1-sim.omega)*mean(Kpredict)+(sim.omega)*k_suggest);
                %                         sim.omega = sim.omega*.99;
                %                     end
                %                     %   kfloat = (1-sim.alpha).*kfloat+sim.beta.* (Kpredict(end)-k_mbest)+sim.gamma.*(Kpredict(end)-k_suggest);
                %                     %  k= round(kfloat+Kpredict(end));
                 u_rnd =randi([sim.umin,sim.umax],1,sim.Nexperts); 
                 v_rnd = randi([sim.vmin,sim.vmax],1,sim.Nexperts);
%                  u_rnd(end) = upredict(end);
%                  v_rnd(end) =vpredict(end);
               %   BER_rnd = Berr(randperm(length(Berr),sim.Nexperts));
                k = k_suggest;
                if(k<1)
                    k=1;
                end
                if(k>sim.k_max)
                    k = sim.k_max;
                end
                %          end
                % Kpredict(ll) = k;
                %   xx(ll) = xx;
                length(t_t);
                fprintf('iteration :%d u:%f v:%f BER:%f k:%f \n',lll,upredict(end),vpredict(end),BERpredict(end),k);
                %          fprintf('iteration :%d u:%f v:%f BER:%f k:%f E_k:%f E_k+1:%fE_k-1:%f mink:%f\n',ll,upredict(end),vpredict(end),BERpredict(end),k,EE_prediction_vector5,EE_prediction_vector6,EE_prediction_vector4,kres(keeres));
                lll = lll+1;
                Kpredict(lll) = k;
                if(sim.channel_change_flag==1)
                    datetime();
                    sim.channel_change_flag =0;
                    Kpredict = [first_k Kpredict];
                    %   LastSymbol_flag = 1;
                    uuu = sim.u;
                    vvv = sim.v;
                    BBB = sim.ch.BER;
                     [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,uuu,vvv,BBB,1);
                    %[EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,uuu(sim.Nparamset-1),vvv(sim.Nparamset-1),BBB(sim.Nparamset-1),1);
                    [xx,kAnal]= min(EE_prediction_vector1);
                    mode(Kpredict)
                    %vars = {'Kpredict(end)','BERpredict(end)','upredict(end)','vpredict(end)','EE_real_track(end)','EE_diff_best_track(end)','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    vars = {'Kpredict','BERpredict','upredict','vpredict','EE_real_track','EE_diff_best_track','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    
                    save(['result_for' num2str(Nrho) '.mat']);
                    Nrho = Nrho+1;
                    optimum_k = k;
                    fprintf('u:%f v:%f BER:%f kcalculated:%f kanal:%f\n',upredict(end),vpredict(end),BERpredict(end),k,kAnal);
                    plot(EE_real_track)
                    hold on;
                    plot(EE_real_track+EE_diff_best_track)
                    hold off, semilogy(EE_real_track), semilogy(EE_real_track+EE_diff_best_track);
                    semilogy(EE_real_track),hold all ,semilogy(EE_real_track+EE_diff_best_track);
                    figure
                    plot(Kpredict),hold on, plot(ones(1,length(Kpredict))*kAnal)
                    plot(EE_prediction_vector1)
                    ll = 1;
                end
                if((length(t_t)<k*n))%||ll>sim.Itlimit)
                    datetime()
                    LastSymbol_flag = 1;
                    sim.channel_change_flag =0;
                    mode(Kpredict)
                    %   LastSymbol_flag = 1;
                    uuu = sim.u;
                    vvv = sim.v;
                    BBB = sim.ch.BER;
                    [EE_prediction_vector1,P_s,perr] = EE_Analisor(k_t,Best_prediction,sim,uuu,vvv,BBB,1);
                    [xx,kAnal]= min(EE_prediction_vector1,[],1);% min(EE_prediction_vector1);
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
                    Kpredict = [first_k Kpredict];
                    %  vars = {'Kpredict(end)','BERpredict(end)','upredict(end)','vpredict(end)','EE_real_track(end)','EE_diff_best_track(end)','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    % vars = {'Kpredict','BERpredict','upredict','vpredictend','EE_real_trackend','EE_diff_best_trackend','kAnal','k','EE_prediction_vector1','uuu','vvv'};
                    vars = {'Kpredict','BERpredict','upredict','vpredict','EE_real_track','EE_diff_best_track','kAnal','k','EE_prediction_vector1','uuu','vvv'};
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
            hold on
            plot(kAnal*ones(1,length(Kpredict)));
        end
    end
end