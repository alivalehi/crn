function [sim] = calc_sim_params(sim)

%This program calculates some of the simulation parameters
% Inputs:
%    symbol process rate
%    Minimum packetization interval: Tmin
%    PER
%    Packetization time vector: Tv
%
% Outputs:
%    symbol process average interarrival time
%    Average packetization interval = 2 * Tmin
%    BER
%    Average packet length
%    Channel Rate
%
% Considerations: Fixed BER, Tmin ==> Higher input throughput ==> defines the minimum acceptable channel rate (for error free channel: Rb_min, for error channel: Rch_min)


if ~isfield(sim, 'CodedSystem')
    sim.CodedSystem = 0;
    %sim.coding.BERD = -1; sim.coding.BERH = -1;
    sim.coding.CTD = 0; sim.coding.CTH=0; %0:uncoded 1:covEnc  2:Turbo  3:LDPC
    sim.coding.CRD = 1; sim.coding.CRH=1; %coding rate
    sim.coding.nTbitsD =0; sim.coding.nTbitsH =0;  %number of Tail Bits
end

if ~sim.CodedSystem   %using coding rate of 1 and no tail bits
    sim.coding.impBERD =1; sim.coding.impBERH; sim.coding.CRD = 1; sim.coding.CRH = 1; sim.coding.CRD = 1; sim.coding.nTbitsH=0; sim.coding.nTbitsD=0;
end

consider_data_tail_bits_in_header=1;
if consider_data_tail_bits_in_header
    sim.coding.nTbitsH = sim.coding.nTbitsH+sim.coding.nTbitsD; sim.coding.nTbitsD=0; %Tails bits in header
end

sim.a = 1/sim.lambda;
if (sim.Framing_mode==2)
    Kv = sim.Kv; KvA = [min(Kv): sim.AKv_step: max(Kv)]; sim.KvA=KvA; %for soft graph
end
%***************************************************************************************************************%
%  To present results in terms of PER, First we consider an expected PER, then we calculate corresponding BER
%  to achieve this PER for T_av=2Tmin
%  Note: there is no approximation for this assumption, but note that the actual PER depends on (Rb, N, and H)
%***************************************************************************************************************%
if sim.Framing_mode == 1
    if isempty(sim.ch.Rch) || (sim.ch.Rch==-1)
        sim.ch.Rch = 1*sim.ch.Rch_min;                         %Channel Rate, not defined
    end
    
    
    %actual values for TV
    if ~isempty(sim.ch.Rch) && (sim.ch.BER >=0) && (sim.ch.BER <= 0.5) && (sim.coding.impBERD >=1) && (sim.coding.impBERH >= 1) && (sim.coding.CRD <= 1) && (sim.coding.CRD >= 1/8) && (sim.coding.CRH <= 1) && (sim.coding.CRH >= 1/8)
        %Rch should be given for coded system
        
        sim.Plen_Tmin = sim.Tmin * sim.lambda * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Expected average packet length for each T
        sim.Rb_Tmin = sim.Plen_Tmin/sim.Tmin;                   %Expected value of input bit rate for T=Tmin
        
        sim.PER_Tmin = 1 - (1 - sim.ch.BER) .^ (sim.Plen_Tmin); %Exact expected value of PER for Tmin
        sim.ch.Rch_min = (1/(1-sim.PER_Tmin))* sim.Rb_Tmin;     %The rate required for stability in the worst case corresponding to Tmin (maximum input rate)
        
        if isempty(sim.ch.Rch) || (sim.ch.Rch==-1)    %if Channel Rate not defined
            sim.ch.Rch = 2*sim.ch.Rch_min;    %was sim.ch.Rch/sim.ch.Rch_min=1.1 > 1 to ensure stability
        end
        
        sim.T_av = 2 * sim.Tmin;  %Assumption: the optimal time is twice Tmin
        sim.Plen_av = sim.T_av * sim.lambda * (sim.N /sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Expected average packet length for each T
        
        if (sim.ch.BER == -1)&& (sim.ch.Exp_PER ~= -1) %only Ex_PER is provided
            sim.ch.BER = 1 - (1 - sim.ch.Exp_PER) ^ (1/sim.Plen_av);  %is BER is not explicitly defined use EXP_PER to find the corresponding BER
        elseif (sim.ch.BER == -1)&& (sim.ch.Exp_PER == -1)
            error ('Neither BER Not PER Defined !!!');
        end
        
        
        
        sim.Plenv = sim.Tv * sim.lambda * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Expected average packet length for each T
        sim.Rbv = sim.Plenv ./ sim.Tv;                   %Input bit rate (including header bits) for each T
        sim.PERv = 1 - ((1 - sim.ch.BER/sim.coding.impBERD) .^ (sim.Tv * sim.lambda * (sim.N / sim.coding.CRD))+ sim.coding.nTbitsD)*((1 - sim.ch.BER/sim.coding.impBERH) .^ ((sim.H / sim.coding.CRH) + sim.coding.nTbitsH));
    else
        error('Error in Encoder Parameter!!!');
    end
    
elseif sim.Framing_mode == 2  %number based
    %save temp, pause;
    
    if (sim.ch.BER >=0) && (sim.ch.BER <= 0.5) && (sim.coding.impBERD >=1) && (sim.coding.impBERH >= 1) && (sim.coding.CRD <= 1) && (sim.coding.CRD >= 1/8) && (sim.coding.CRH <= 1) && (sim.coding.CRH >= 1/8) %~
        sim.Plen_min = sim.Kmin * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Packet length for each T
        sim.PER_min = 1 - ((1 - sim.ch.BER/sim.coding.impBERD) .^ (sim.Kmin * (sim.N / sim.coding.CRD))+ sim.coding.nTbitsD)*((1 - sim.ch.BER/sim.coding.impBERH) .^ ((sim.H / sim.coding.CRH) + sim.coding.nTbitsH));
        
        %sim.Plenv    =   sim.Kv * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Packet length for each T
        %sim.PERv    = 1 - ((1 - sim.ch.BER/sim.coding.impBERD) .^   (sim.Kv * (sim.N / sim.coding.CRD))+ sim.coding.nTbitsD)*((1 - sim.ch.BER/sim.coding.impBERH) .^ ((sim.H / sim.coding.CRH) + sim.coding.nTbitsH));
        sim.Plenv    =   KvA * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Packet length for each T
        sim.PERv    = 1 - ((1 - sim.ch.BER/sim.coding.impBERD) .^   (KvA * (sim.N / sim.coding.CRD))+ sim.coding.nTbitsD)*((1 - sim.ch.BER/sim.coding.impBERH) .^ ((sim.H / sim.coding.CRH) + sim.coding.nTbitsH));
        
        %lets assume channel is capable of handling Kmin sambols
        if isempty(sim.ch.Rch)
            input_traffic = sim.Plen_min ./ (sim.Kmin ./ sim.lambda);  %on average at eeach k ./ sim.lambda, there is one packet to send
            output_traffic = input_traffic / (1-sim.PER_min);
            sim.ch.Rch = 2 * output_traffic; %Channel rate is adjusted for K > KinP_PB: Lower K might result in unstable queue
        end
        

    else
        save errdata_calc_sim_params
        error('Error in Encoder Parameter!!!');
    end
end
        if sim.cognitive
            if ~isfield(sim.CRN, 'chNotAvRatio'), sim.CRN.chNotAvRatio = 1; end
            if ~isfield(sim.CRN, 'chAvPlenRatio'), sim.CRN.chAvPlenRatio = 40; end
        else
            sim.CRN.chNotAvRatio = 0;
            sim.CRN.chAvPlenRatio = 1;
        end
        if(sim.ch.channel_AV_B_proportional)
            if ~isfield(sim.CRN,'Ch_AvailMeanTime')&& sim.CRN.chAvPlenRatio~= 0
                sim.CRN.Ch_AvailMeanTime = sim.CRN.chAvPlenRatio * sim.Plenv/sim.ch.Rch;
            end
            
            if isempty(sim.CRN.Ch_AvailMeanTime) && sim.CRN.chAvPlenRatio~= 0
                sim.CRN.Ch_AvailMeanTime = sim.CRN.chAvPlenRatio * sim.Plenv/sim.ch.Rch;
            end
            sim.CRN.Ch_NotAvailMeanTime = sim.CRN.chNotAvRatio * sim.CRN.Ch_AvailMeanTime;
        else
            sim.CRN.Ch_AvailMeanTime = sim.CRN.chAvPlenRatio * ones(size(sim.Plenv/sim.ch.Rch));
            sim.CRN.Ch_NotAvailMeanTime = sim.CRN.chNotAvRatio * ones(size(sim.CRN.Ch_AvailMeanTime));
        end
