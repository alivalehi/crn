function [sim] = init_sim(mode)
%mode:1 Time Based   [All symbols arrive at a time-slot are bundled into a packet]
%mode:2 Number Based [K Consecutive arrival symbols are bundled into a packet]
sim.SU = 1;
sim.datapath = regexprep([cd , '/PCodeN/NData/'], '\', '/');  %was sim.datapath = regexprep([cd , '/data/'], '\', '/')
sim.fname = '';
addpath(sim.datapath);

%sim.CRN.Ch_AvailMeanTime = int16.empty;
sim.control.save_active = 1;
sim.control.plot_active = 0;
sim.control.debug_active = 0;
sim.control.runtime_plot_active=0;
sim.control.overwrite_active = 0;
sim.control.analysis = 1; %set always to 1
sim.control.simulation = 1;


sim.Framing_mode = mode;  %1:Time Based, 2:Number Based;
sim.input_process = 1; %1:Poisson, 0:Deterministic
sim.UFCPP=1;%use_framing_to_calc_packet_params = 1; %0 means directly use poisson
sim.H = 128; 
sim.N = 8;
sim.lambda = .01;        %input symbol rate [symbol/sec]


sim.run.Nsym = 1e3;      %number of symbols to run simulations
sim.run.minSym = 1e3;    %minimum number of symbols to run simulations
sim.run.maxSym = 1e7;  %maximum number of symbols to run simulations
sim.run.minPack = 1e2;   %minimum number of packets to run simulations
sim.run.maxPack = 1e7;   %minimum number of packets to run simulations
sim.run.NUM_RUNS = 4;  
%sim.run.quick_run_coef = 0.01;

sim.run.Last_Percentage = [10 25 50 75 90];%[10 25 50 75 100];  %Let queue stabilize and then calc average waiting and service times 

if mode == 2 %number based method
    sim.Kv = [1:1:10,    12:2:20, 25:5:100,   200:100: 500];
    sim.ch.BER = 1e-4;
    sim.Kmin = 5;
    sim.AKv_step = 1;
else
    sim.send_dummy_for_k0 = 0;
    sim.Tmin = 5 * (1/sim.lambda); %5 times average interarrival duration
    % Channel rate is chosen to be capable of handling traffic in this T
    % (T_TB = 5*a): T<5a ==> Unstable Queue due to large overhead;
    sim.Tv = (1/sim.lambda)*[[1 2 5]*1e-4, [1 2 5]*1e-3, [1 2 5]*1e-2, [1 2 5 8] *0.1,    1:0.5:30,        31:1:40,       45:5:70,         80:10:100,150,200]; %for BER use
    sim.Tv2=sim.Tv;
end
sim.ch.Rch =[]; %not defined, calc based on Tmin or Kmin
sim.ch.Rch0 =[]; %not defined, calc based on Tmin or Kmin

sim.ch.Exp_PER = -1; %not defined
%sim.ch.BER = -1;     %not defined

sim.CodedSystem = 0; 
sim.coding.impBERD = 10; sim.coding.impBERH = 10;  %reduction in BER  effectiveBER = BER/BERimp [Use Refrences or Use Monte Carlo Simulations to calc this parameter]
sim.coding.CTD = 1; sim.coding.CTH=1; %0:uncoded 1:covEnc  2:Turbo  3:LDPC
sim.coding.CRD = 1/2; sim.coding.CRH=1/2; %coding rate
sim.coding.nTbitsD =3; sim.coding.nTbitsH =3;  %number of Tail Bits



sim.cognitive = 1;
sim.CRN.chNotAvRatio = 5;  %not available mean time to available mean time ratio [0 means single network]
sim.CRN.chAvPlenRatio =2 ; 
sim.CRN.Approx = 0;   %use approximation of chAvailTime,chNotAvailTime >> Plen/Rch

%FramingModeStr = {'Time Based', 'Number Based'};
%sim = fill_dummy_sim_params(sim);
%fprintf('Parameters: N:%d  H:%d   Rate:%1.3f  Mean Inbterval Time:%1.3f \n', sim.N,sim.H,sim.lambda,1/sim.lambda);
%fprintf('Framing Mode:%s     Packetization (Time):%1.3f   Packetization (K):%d\n', FramingModeStr{sim.Framing_mode}, sim.T_TB, sim.NinP_PB);

%sim.method=2;
end