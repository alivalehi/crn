function [queue] = f_q_evolution(queue, sim, TKind)  %was function [wnv, sv, rv] = f_q_evolution(tp, kp, N, H, BER, Rch, send_dummy_for_k0)
%This function receives input packets, simulation parameter, policy of channel
% and simulate queue of packets in a network
% Inputs:
%    queue: contains tp(timestamp of end of packets),kp(number of symbols in each packet)
%    sim: simulation parameter
%    TKind: Index to find right parameter of the channel it is useful when
%    the expected availability of each slot of the channel is proportional
%    to the length of packets
% Output:
%        queue contains waiting time and service time
try
tp = queue.tp; kp = queue.kp; queue.fail=0;

if sim.Framing_mode==1
    %remove packets with zero symbols
    if ~sim.send_dummy_for_k0
        nzind = find(kp ~= 0);
        tp = tp(nzind);
        kp = kp(nzind);
    end
    npackets = length(kp);
    
    if sim.CodedSystem
        Plen = kp * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH;
        PER = 1 - ((1 - sim.ch.BER/sim.coding.impBERD) .^ (kp * (sim.N / sim.coding.CRD))+ sim.coding.nTbitsD)*((1 - sim.ch.BER/sim.coding.impBERH) .^ ((sim.H / sim.coding.CRH) + sim.coding.nTbitsH));
    else
        Plen = kp * sim.N + sim.H;
        PER = 1 - (1 - sim.ch.BER) .^ (Plen);
    end
else
    %packet lengths are constant in this mode
    npackets = length(kp);
    Plen = sim.Plenv(TKind);
    PER = sim.PERv(TKind);
end

rv = 1+geornd(1-PER,1,npackets);     %retransmission in ARQ
sv = (Plen .* rv)/sim.ch.Rch;           %service time
sv_1tranmition = (Plen /sim.ch.Rch);           %service time for one transmission
if sim.Framing_mode==2, sv1=sv_1tranmition; end  %service time for one transmission

wnv = zeros(1,npackets);
wnv(1)=0;
for p = 1: npackets-1
    inter_pack_time = tp(p+1)-tp(p);
    wnv(p+1) = max(0, wnv(p) + sv(p) - inter_pack_time);  %inter-arrival time waiting time according to Lindley equation
end
queue.rv = rv; queue.sv = sv; queue.sv1 = sv_1tranmition;  queue.wnv = wnv;

CRNwnv = []; CRNsv =[];
if sim.cognitive
    if isfield(sim, 'CRN')
        if isfield(sim.CRN, 'Ch_AvailMeanTime') && isfield(sim.CRN, 'Ch_NotAvailMeanTime')

            a = sim.CRN.Ch_AvailMeanTime(TKind); b=sim.CRN.Ch_NotAvailMeanTime(TKind);
         %   a = sim.v; b=sim.u;
            nch = 100* ceil(max(tp) / (a+b));
            nch = min(1e8, max(1e5, nch));
            sim = channel_perform(sim,nch);
            chf = sim.chf;
            cht = sim.cht;
            chs = sim.chs;
            CRNwnv = zeros(1,npackets); CRNsv = zeros(1,npackets);
            CRNwnv(1) = 0;
            
            %save temp1, disp('temp1 saved'); pause(2);
            npackets
            tp(npackets+1)=tp(npackets)+1e-10; %just to avoid out of bound error
            %begining of the actual transmission process
            for p = 1: npackets
                if sim.Framing_mode==1, sv1=sv_1tranmition(p); end
                p
                %calculate service time
                sst0 = tp(p) + CRNwnv(p);   %service start time for all retransmissions
                sst = sst0;               %service start time
                
                for r = 1 : rv(p)   %try for all retransmission
                    
                    send_flag = 0; ret_cont=0;
                    while (~send_flag && ret_cont<100)
                        ret_cont=ret_cont+1;
                        cht_s1 = find(cht <= sst,1,'last');
                        cht_s2 = find(cht > sst,1,'first');
                        cht_e1 = find(cht <= sst+sv1,1,'last');
                        cht_e2 = find(cht > sst+sv1,1,'first');
                        
                        if sim.control.debug_active, fprintf('p:%d r:%d retry:%d  arrival:%1.4f wt:%1.4f   startp:%1.4f  pt[%1.4f-%1.4f], start(%d-%d) end(%d-%d) ', p,r,ret_cont, tp(p), CRNwnv(p), sst0, sst, sst+sv1, cht_s1, cht_s2, cht_e1, cht_e2); end
                        
                        if isempty(cht_e2)  %extend channel availabity process and retry
                            chf = [chf, chf(2:end)];
                            cht = [cht, max(cht)+ cht(2:end)];
                            if sim.control.debug_active, fprintf('  ch extended \n'); end
                            
                        elseif  ((cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 == cht_e1) && (chf(cht_s1) == 0))  %packet transmission fits into one availale time slot
                            sst = sst+ sv1;  %set start time for next retransmission
                            send_flag = 1;   %mark the packet as transmitted
                            if sim.control.debug_active, fprintf('  OK \n'); end
                        elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (chf(cht_s1) == 1)  %packet transmission arrive at busy time slot
                            next_av_ch = find((cht > sst) & (chf == 0)& (chs>sv1),1,'first');
                            sst = cht(next_av_ch);
                        elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 ~= cht_e1)  && (chf(cht_s1) == 0)  %packet transmission does not fits into one availale time slot
                            next_av_ch = find((cht > sst) & (chf == 0)& (chs>sv1),1,'first');
                            sst = cht(next_av_ch);
                        else
                            display('Unknown condition');
                            pause;
                        end
                    end  %end while
                end %end for r = 1 : rv(p)
                if ret_cont > 100
                    queue.fail=1;
                    break; %FAIL
                else
                    CRNsv(p) = sst - sst0;
                    inter_pack_time = tp(p+1)-tp(p);
                    CRNwnv(p+1) = max(0, CRNwnv(p) + CRNsv(p) - inter_pack_time);  %inter-arrival time
                    if sim.control.debug_active, fprintf('  st(p):%1.4f   T(p+1):%1.4f wt(p):%1.4f    wt(p+1):%1.4f \n',CRNsv(p), inter_pack_time,  CRNwnv(p), CRNwnv(p+1)); end
                end
            end  %for p = 1: npackets
            CRNwnv=CRNwnv(1:end-1);
            queue.CRNwnv = CRNwnv;  queue.CRNsv = CRNsv;  queue.cht = cht;  queue.chf = chf;
        end
    end 
end
catch
    display('Please change the simulation parameters');
end

