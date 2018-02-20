%read content of sim structure

send_dummy_for_k0 = sim.send_dummy_for_k0;
maxSym = sim.maxSym;
minSym = sim.minSym;
maxPack = sim.maxPack;

H = sim.H;
N = sim.N;
%Tmin_TB = sim.T_TB;
%Nsym = sim.Nsym;

Tv = sim.Tv;

lambda = sim.lambda;     %Average Interarrival Time  ?
NUM_RUNS = sim.NUM_RUNS; 
input_process = sim.input_process;

%Channel Parameters
BER = sim.BER;
Rch = sim.Rch;

Plenv = sim.Plenv;    %Packet Length
PERv = sim.PERv;
%Rb_av = sim.Rb_av;
%sim.info_bit_rate = sim.lambda * sim.N;

Framing_mode = sim.Framing_mode;
use_framing_to_calc_Kp = sim.use_framing_to_calc_Kp;

plot_wnv_active = sim.plot_queu_active;
plot_active = sim.plot_active;
debug_active = sim.debug_active;
save_active = sim.save_active;

Last_Percentage=sim.Last_Percentage; %razi added 6-20-2015