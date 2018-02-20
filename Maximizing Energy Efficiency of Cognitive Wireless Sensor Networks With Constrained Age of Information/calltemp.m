%This program runs simulations for different parameter setups
%see also call_framing_delay for similar test setups

clear all
close all
%clc
%#######################################################################
%     SIMULATION OPTIONS
%#######################################################################

%sim.Nsym = 10000;
%sim.NUM_RUNS = 10;
%sim.H = 20;
%sim.Exp_PER = 0.1;

run_test = 1;
if run_test 
    testsim = init_sim;
    sim.plot_active = 1;
    testsim.H = 30;
    testsim.N = 16;
    testsim.Nsym = 1000;
    testsim.Exp_PER = 0.2;
    testsim.NUM_RUNS = 10;
    testsim.Tv = [1, 1.5, 2, 3:3:10, 15:5:20, 30:10:100] * testsim.a;

    %fqueu_TB(testsim, 1, 0);
    fqueu_TB(testsim);

    plot_data_file(testsim.fname, 0);

end    

%#######################################################################
%     SIMULATE FOR DIFFERENT H
%#######################################################################

% sim = init_sim;
% sim.Nsym = 1000;
% fqueu_TB(sim, 1, 0);
if ispc && ~run_test
    disp('pc runs');
    pause;
    H = [10 20 30 40 50];
    N = 16;
    Nsym = 50000;
    Exp_PER = 0.2;
    NUM_RUNS = 20;
    for i=2:5
        sim(i) = init_sim;
        sim(i).N = N;
        sim(i).H = H(i);
        sim(i).Nsym = Nsym;
        sim(i).Exp_PER = Exp_PER;
        sim(i).NUM_RUNS = NUM_RUNS;
        %fqueu_TB(sim(i), 1, 0);
        fqueu_TB(sim(i));
        save dall_pc
    end


    H = [10 20 30 40 50];
    N = 8;
    Nsym = 10000;
    Exp_PER = 0.2;
    NUM_RUNS = 50;

    for i=6:10
        sim(i) = init_sim;
        sim(i).N = N;
        sim(i).H = H(i-5);
        sim(i).Nsym = Nsym;
        sim(i).Exp_PER = Exp_PER;
        sim(i).NUM_RUNS = NUM_RUNS;
        %fqueu_TB(sim(i), 1, 0);
        fqueu_TB(sim(i));
        save dall_pc
    end
end


if ~ispc && ~run_test
    disp('mac runs');
    pause;

    H = 30;
    N = 16;
    Nsym = 50000;
    Exp_PER = [0, 0.05, 0.1, 0.2: 0.1: 0.5];
    NUM_RUNS = 20;
    for i=11:17
        sim(i) = init_sim;
        sim(i).N = N;
        sim(i).H = H;
        sim(i).Nsym = Nsym;
        sim(i).Exp_PER = Exp_PER(i-10);
        sim(i).NUM_RUNS = NUM_RUNS;
        %fqueu_TB(sim(i), 1, 0);
        fqueu_TB(sim(i));
        save dall_mac
    end


    H = 30;
    N = 8;
    Nsym = 10000;
    Exp_PER = [0, 0.05, 0.1, 0.2: 0.1: 0.5];
    NUM_RUNS = 50;
    for i=21:25
        sim(i) = init_sim;
        sim(i).N = N;
        sim(i).H = H;
        sim(i).Nsym = Nsym;
        sim(i).Exp_PER = Exp_PER(i-20);
        sim(i).NUM_RUNS = NUM_RUNS;
        %fqueu_TB(sim(i), 1, 0);
        fqueu_TB(sim(i));
        save dall_mac
    end
end
