load temp1
sim.debug_active=1;

tp(npackets+1)=tp(npackets)+1e-10; %just to avoid out of bound error

tp(1:15)
kp(1:15)
rv(1:15)
sv1
[cht(1:15); chf(1:15)]

rv = rv * 2; queue.rv=rv;
sv1=sv1/1;

        for p = 1: npackets
            if sim.Framing_mode==1, sv1=sv_1tranmition(p); end
            
            %calculate service time
            sst0 = tp(p) + wnv2(p);   %service start time for all retransmissions
            sst = sst0;               %service start time
            
            for r = 1 : rv(p)   %try for all retransmission
                send_flag = 0; ret_cont=0;
                while (~send_flag && ret_cont<100)
                    ret_cont=ret_cont+1;
                    cht_s1 = max(find(cht <= sst)); cht_s2 = min(find(cht > sst));
                    cht_e1 = max(find(cht <= sst+sv1)); cht_e2 = min(find(cht > sst+sv1));
if sim.debug_active, fprintf('p:%d r:%d retry:%d  arrival:%1.4f wt:%1.4f   startp:%1.4f  pt[%1.4f-%1.4f], start(%d-%d) end(%d-%d) ', p,r,ret_cont, tp(p), wnv2(p), sst0, sst, sst+sv1, cht_s1, cht_s2, cht_e1, cht_e2); end


                    if isempty(cht_e2)  %extend channel availabity process and retry
                        chf = [chf, chf(2:end)];
                        cht = [cht, max(cht)+ cht(2:end)];
if sim.debug_active, fprintf('  ch extended \n'); end
                    elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 == cht_e1) && (chf(cht_e2) == 1)  %packet transmission fits into one availale time slot  
                        sst = sst+ sv1;  %set start time for next retransmission
                        send_flag = 1;   %mark the packet as transmitted
if sim.debug_active, fprintf('  OK \n'); end
                    
                    else %for all other conditions check the first upcoming available channel
                        next_av_ch = min(find((cht > sst) & (chf == 0))); 
                        if isempty(next_av_ch)  % extend the ch. avail process
                            chf = [chf, chf(2:end)];
                            cht = [cht, max(cht)+ cht(2:end)];
                            next_av_ch = min(find((cht > sst) & (chf == 0))); 
                        end
                        sst = cht(next_av_ch); 
if sim.debug_active, fprintf('  Retry at time(%d):%1.4f  \n',next_av_ch, sst); end
                        if isempty(next_av_ch), queue.fail=1; end
                    end
                end  %end while
            end %end for r = 1 : rv(p)
            if ret_cont > 100
                queue.fail=1;
                break; %FAIL
            else
                sv2(p) = sst - sst0;
                inter_pack_time = tp(p+1)-tp(p);
                wnv2(p+1) = max(0, wnv2(p) + sv2(p) - inter_pack_time);  %inter-arrival time
if sim.debug_active, fprintf('  st(p):%1.4f   T(p+1):%1.4f wt(p):%1.4f    wt(p+1):%1.4f \n',sv2(p), inter_pack_time,  wnv2(p), wnv2(p+1)); end
            end
        end  %for p = 1: npackets
        wnv2=wnv2(1:end-1);
        queue.wnv2 = wnv2;  queue.sv2 = sv2;  queue.cht = cht;  queue.chf = chf;  

