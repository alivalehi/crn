%This program runs simulations for different parameter setups
%see also call_framing_delay for similar test setups

clear all
%clc
%close all
%#######################################################################
%     SIMULATION OPTIONS
%#######################################################################

%sim.Nsym = 10000;
%sim.NUM_RUNS = 10;
%sim.H = 20;
%sim.Exp_PER = 0.1;


%some example setups:
%N=16, H= 20 ==>                0.1 <P< 0.5
%N=16, H= 40 ==> P = 0.2,       0.1 <P< 0.5
%N=16, H= 100 ==> P = 0.4,0,5   0.1 <P< 0.5
%N=16, H= 200 ==> P = 0.5,      0.2 <P< 0.6


runall = 0;
run_comp_PER = 0;
run_test = 1;  %single simulation to cpmpare simulation and analytical



if run_test 
    
    sim1 = init_sim;
    
    sim1.plot_active = 0;
    sim1.debug_active = 0;
    sim1.use_framing_to_calc_Kp = 1;  %simulate input packet process

    
    sim1.N = 16;   
    sim1.H = 40;   %or 30
    %sim1.Exp_PER = 0.2; %impose appropriate channel rate if BER is not defined explicitely
    
    %parameters for long run time
    sim1.Tv = [1:0.2:40, 41:1:50] * sim1.a; %may take a long run
    sim1.Tv2 = sim1.Tv;                     %used for analytical calculations and plots
    sim1.fname='TempTest.mat'; %was sim1.fname='Es_Ew_comp.mat';
    

    sim1.minSym =  50;        
    sim1.maxSym =  50000;
    sim1.maxPack = 5000;
    sim1.NUM_RUNS = 100;

    %Set channel parameters
    sim1.Exp_PER=-1;
    sim1.Rch=300;
    sim1.BER=0.001;
    
    
    fqueu_TB(sim1);
    sim1.plot_active = 1;
    plot_data_file(sim1.datapath, sim1.fname , 0);
    legend('E[f]: (Analytical)','E[w]: (Analytical)','E[d]: (Analytical)', 'E[f]: (Simulation)','E[w]: (Simulation)','E[d]: (Simulation)');
    %ylim([0 12]);  xlim([0 4]); saveas(gcf,'Fd_sym2','fig'); 

end    






run_specific_test = 0;
if run_specific_test 

    
    sim1 = init_sim;
    sim1.N = 16;
    sim1.H = 40;
    sim1.Exp_PER = 0.2;
    sim1.NUM_RUNS = 5;
    sim1.use_framing_to_calc_Kp = 0;
    sim1.maxPack=sim1.maxPack*10;
    sim1.fname='TBP_N16_H40_PER20Test';
    fqueu_TB(sim1);
    
    
pause;
pause;
pause;
    
    sim1 = init_sim;
    sim1.N = 8;
    sim1.H = 10;
    sim1.Exp_PER = 0.2;
    sim1.NUM_RUNS = 2;
    sim1.use_framing_to_calc_Kp = 1;
    %sim1.Tv = 0.01/sim1.lambda;
    sim1.maxPack=sim1.maxPack/10;
    sim1.fname='TBP_N8_H10_PER20Test'
    fqueu_TB(sim1);
pause;

end    


if run_comp_PER && ~run_test
    NUM_RUNS = 10;
    sim.simulation = 0;
    sim.analysis = 1;
    N = 16;
    H = 30;
    %Exp_PERv = [0, 0.1, 0.12, 0.14, 0.15, 0.16, 0.18, 0.2];
    Exp_PERv = [0, 0.12, 0.14, 0.15];
    for i=1:length(Exp_PERv)  
        sim = init_sim;
        sim.N = N;
        sim.H = H;
        sim.Exp_PER = Exp_PERv(i);
        sim.simulation = 0;
        sim.analysis = 1;
        fqueu_TB(sim);
    end
end




%#######################################################################
%     SIMULATE FOR DIFFERENT H
%#######################################################################

% sim = init_sim;
% sim.Nsym = 1000;
% fqueu_TB(sim, 1, 0);
if runall && ~run_test
    disp('pc runs for N = 16');

    NUM_RUNS = 10;
    Nsym = 50000;
    
    Nv = [8 16];
    Hv = [10 20 30 40 50 100 200];
    Exp_PERv = [0, 0.8  0.5  .2  .6  .1 .4 0.01, 0.001];
    
    N = 16;
    Exp_PER = 0.2;
    for i=1:length(Hv)  
        sim = init_sim;
        sim.N = N;
        sim.H = Hv(i);
        %sim.Nsym = Nsym;
        sim.Exp_PER = Exp_PER;
        fqueu_TB(sim);
    end
    
    N = 16;
    H = 20;
    for i=1:length(Exp_PERv) 
        sim = init_sim;
        sim.N = N;
        sim.H = H;
        sim.Exp_PER = Exp_PERv(i);
        fqueu_TB(sim);
    end

    N = 16;
    H = 40;
    for i=1:length(Exp_PERv)    
        sim = init_sim;
        sim.N = N;
        sim.H = H;
        sim.Exp_PER = Exp_PERv(i);
        fqueu_TB(sim);
    end
    
    N = 16;
    H = 100;
    for i=1:length(Exp_PERv)    
        sim = init_sim;
        sim.N = N;
        sim.H = H;
        sim.Exp_PER = Exp_PERv(i);
        fqueu_TB(sim);
    end
    
    N = 16;
    H = 200;
    for i=1:length(Exp_PERv)    
        sim = init_sim;
        sim.N = N;
        sim.H = H;
        sim.Exp_PER = Exp_PERv(i);
        fqueu_TB(sim);
    end
    
end


%%
if runall && ~run_test
    disp('mac runs');

    Hv = [10 20 30 40 50 100 200];
    Exp_PERv = [0, 0.8  0.5  .2  .6  .1 .4 0.01, 0.001];
    
    N = 8;
    Exp_PER = 0.2;
    for i=1:length(Hv)  %investigate H effect
        sim = init_sim;
        sim.N = N;
        sim.H = Hv(i);
        sim.Exp_PER = Exp_PER;
        fqueu_TB(sim);
    end
    
    N = 8;
    for H = [10, 20, 40, 50, 200]
        for i=1:length(Exp_PERv) 
            sim = init_sim;
            sim.N = N;
            sim.H = H;
            sim.Exp_PER = Exp_PERv(i);
            fqueu_TB(sim);
        end
    end
    
end
    