function [queue] = f_q_evolution(queue, sim, TKind, SU)
tp = queue.tp; kp = queue.kp; queue.fail=0;
sim.control.debug_active=0;
lambda = sim.lambda;
debug_f_q = 1;
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
        PER = 1 - (1 - sim.ch.BER) .^ (sim.Plenv);
    end
else
    %packet lengths are constant in this mode
    npackets = SU*length(kp);
    Plen = sim.Plenv(TKind);
    PER = sim.PERv(TKind);
end

rv = 1+geornd(1-PER,1,npackets);     %retransmission in ARQ
%rv =ones(1,npackets);
sv = (Plen .* rv)/sim.ch.Rch;           %service time
sv_1tranmition = (Plen /sim.ch.Rch);           %service time for one transmission
if sim.Framing_mode==2, sv1=sv_1tranmition; end  %service time for one transmission

wnv = zeros(1,npackets);
wnv(1)=0;
for p = 1: npackets-1
    inter_pack_time = tp(p+1)-tp(p);
    wnv(p+1) = max(0, wnv(p) + sv(p) - inter_pack_time);  %inter-arrival time
end
queue.rv = rv; queue.sv = sv; queue.sv1 = sv_1tranmition;  queue.wnv = wnv;

CRNwnv = []; CRNsv =[];
if sim.cognitive
    if isfield(sim, 'CRN')
        if isfield(sim.CRN, 'Ch_AvailMeanTime') && isfield(sim.CRN, 'Ch_NotAvailMeanTime')
            %  a = sim.CRN.Ch_AvailMeanTime(1); b=sim.CRN.Ch_NotAvailMeanTime(1);
             a = sim.CRN.Ch_AvailMeanTime(TKind); b=sim.CRN.Ch_NotAvailMeanTime(TKind);
       %     b = sim.CRN.chNotAvRatio;        a = sim.CRN.chAvPlenRatio;
            nch = 100* ceil(max(tp) / (a+b));
            nch = min(1e8, max(1e5, nch));
            
            %generate channel availability durations
            avt  = exprnd(a , [1, nch]);     %duration with available channel
            avnt = exprnd(b, [1, nch]);  %duration with no available channel
            avnvt = reshape([avt; avnt], [1, 2*nch]);
            
            %channel availability flag
            chf = reshape([ones(1,nch); zeros(1,nch)], [1, 2*nch]);
            
            %channel availability process
            cht = nan(size(avnvt)); cht(1) = 0;
            for i = 1: length(cht)
                cht(i+1)=cht(i)+avnvt(i);
            end
            chf = [0, chf];
            chcont = [1:length(chf)];
            chs = [avnvt,0];
            CRNwnv = zeros(1,npackets); CRNsv = zeros(1,npackets);
            CRNwnv(1) = 0;
            
            %save temp1, disp('temp1 saved'); pause(2);
            %%%Debug part%%%
            if(debug_f_q == 1)
                P = zeros(1,npackets); % # of insufficient available when arrive at insufficient available
                P3 = 0;% # of  available  arrive
                P1 = 0;% # of insufficient available arrive
                P2 = 0;% # of busy arrive
                P4 = zeros(1,npackets);% # of insufficient available when arrive at busy
                mean_first_u = zeros(1,npackets);%just for check mean of first busy slot in insufficient available time slot
                mean_u_v_available = zeros(1,npackets);%just for check mean of first busy slot in insufficient available time slot
                mean_first_v = [];%zeros(1,npackets); % man of starting at sav (phi)
                mean_u_v_busy = zeros(1,npackets);
                temp = -ones(1,npackets);
                first_slot = zeros(1,npackets);
                first_slot2 = zeros(1,npackets); %sav
                first_slot4 = zeros(1,npackets); %busy
                status_flag = zeros(1,npackets);
                last_slot = zeros(1,npackets);
            end
            %%%EnD of Debug part%%%
            
            tp(npackets+1)=tp(npackets)+1e-10; %just to avoid out of bound error
            
            for p = 1: npackets
                if sim.Framing_mode==1, sv1=sv_1tranmition(p); end
                
                %calculate service time
                sst0 = tp(p) + CRNwnv(p);   %service start time for all retransmissions
                sst = sst0;               %service start time
                %%%Debug part%%%
                if(debug_f_q == 1)
                    first_flag = 0;
                    long_sav_flag = 0;
                    long_sav_flag2 = 0;
                end
                %%%EnD of Debug part%%%
                for r = 1 : rv(p)   %try for all retransmission
                    
                    send_flag = 0; ret_cont=0;
                    while (~send_flag && ret_cont<100)
                        ret_cont=ret_cont+1;
                        cht_s1 = find(cht <= sst,1,'last');
                        cht_s2 = find(cht > sst,1,'first');
                        cht_e1 = find(cht <= sst+sv1,1,'last');
                        cht_e2 = find(cht > sst+sv1,1,'first');
                        sstf = sst;
                        if sim.control.debug_active, fprintf('v:%f u:%f s:%f p:%d r:%d retry:%d  arrival:%1.4f wt:%1.4f   startp:%1.4f  pt[%1.4f-%1.4f], start(%d-%d) end(%d-%d) ',a,b,sv1, p,r,ret_cont, tp(p), CRNwnv(p), sst0, sst, sst+sv1, cht_s1, cht_s2, cht_e1, cht_e2); end
                        
                        if isempty(cht_e2)  %extend channel availabity process and retry
                            chf = [chf, chf(2:end)];
                            cht = [cht, max(cht)+ cht(2:end)];
                            chs = [chs,chs(2:end)];
                            if sim.control.debug_active, fprintf('  ch extended \n'); end
                            
                        elseif  ((cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 == cht_e1) && (chf(cht_s1) == 0))  %packet transmission fits into one availale time slo
                            %%%Debug part%%%
                            if(debug_f_q == 1)
                                if(~first_flag)
                                    status_flag(p) = 1;
                                    first_flag     = 1;
                                    P3= P3+1;
                                end
                                
                                if(status_flag(p)== 2)%sav
                                    last_slot(p)= cht_s1;
                                    %   if(P4(p)~=0)
                                    temp(p) = (last_slot(p)-(first_slot2(p)+1));%->(id of start slot[when packet sent ])-id of end slot[when packet arrived]+1
                                    if(temp(p)~=0)
                                        mean_u_v_available(p) = (cht(last_slot(p)) - cht(first_slot2(p) + 1));%/(temp(p)/2);%->(cht_s1[when packet sent ])-(cht_s2+1)[when packet arrived]
                                    end
                                    %end
                                end
                                if(status_flag(p)== 3)%busy
                                    last_slot(p)= cht_s1;
                                    
                                    %  temp(p) = (last_slot(p)-first_slot4(p));%(id of start slot[when packet sent ])-id of end slot[when packet arrived]
                                    %  if(P(p)~=0)
                                    if(temp(p)~=0)
                                        mean_u_v_busy(p) = (cht(last_slot(p))-cht(first_slot4(p)));%->(cht_s1[when packet sent ])-(cht_s2)[when packet arrived]
                                        mean_first_u_busy(p) = (cht(first_slot4(p))-sst0);%->small part of busy begining of next slot-arrival moment
                                    end
                                    % end
                                end
                            end
                            %%%EnD of Debug part%%%
                            sst = sst+ sv1;  %set start time for next retransmission
                            send_flag = 1;   %mark the packet as transmitted
                            if sim.control.debug_active, fprintf('  OK \n'); end
                        elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (chf(cht_s1) == 1)  %packet transmission arrive at busy time slot
                            
                            next_av_ch = find((cht > sst) & (chf == 0)& (chs>sv1),1,'first');
                            if isempty(next_av_ch)  % extend the ch. avail process
                                chf = [chf, chf(2:end)];
                                cht = [cht, max(cht)+ cht(2:end)];
                                chs = [chs,chs(2:end)];
                                next_av_ch = find((cht > sst) & (chf == 0) & (chs>sv1),1,'first');
                            end
                            sst = cht(next_av_ch);
                            %%%Debug part%%%
                            if(debug_f_q == 1)
                                j=1;
                                if(~first_flag)
                                    status_flag(p)= 3 ;
                                    first_flag = 1;
                                    first_slot(p)= cht_s2;
                                    first_slot4(p)= cht_s2;
                                    first_vb(p)= cht(cht_s2+1)-cht(cht_s2);
                                    P2= P2+1;
                                    % first_slot(j) = slot_s2+1;
                                end
                            end
                            %%%EnD of Debug part%%%
                            if isempty(next_av_ch), queue.fail=1; end
                            if sim.control.debug_active, fprintf('  Retry at time(%d):%1.4f  \n',next_av_ch, sst); end
                        elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 ~= cht_e1)  && (chf(cht_s1) == 0)  %packet transmission does not fits into one availale time slot
                            next_av_ch = find((cht > sst) & (chf == 0)& (chs>sv1),1,'first');
                            if isempty(next_av_ch)  % extend the ch. avail process
                                chf = [chf, chf(2:end)];
                                cht = [cht, max(cht)+ cht(2:end)];
                                chs = [chs,chs(2:end)];
                                next_av_ch = find((cht > sst) & (chf == 0)&(chs>sv1),1,'first');
                            end
                            sst = cht(next_av_ch);
                            %%%Debug part%%%
                            if(debug_f_q == 1)
                                j=1;
                                if(~first_flag)
                                    status_flag(p) = 2;
                                    first_flag     = 1;
                                    first_slot2(p)= cht_s2;
                                    mean_first_u(p) = cht(cht_s2+1)-cht(cht_s2);
                                    first_ub(p)= cht(cht_s2+3)-cht(cht_s2+2);
                                    mean_first_v(end+1) = cht(cht_s2)-sstf;
                                    P1= P1+1;
                                else
                                    if(  status_flag(p)== 3)
                                        if(long_sav_flag==0)
                                            long_sav_flag = 1;
                                        end
                                        P(p)= P(p)+1;
                                    else
                                        
                                        if(long_sav_flag2==0)
                                            first_slot3(p)= cht_s1;
                                            long_sav_flag2 = 1;
                                        end
                                        P4(p)= P4(p)+1;
                                    end
                                end
                            end
                            %%%EnD of Debug part%%%
                            if isempty(next_av_ch), queue.fail=1;pause; end
                            if sim.control.debug_active, fprintf('  Retry at time(%d):%1.4f  \n',next_av_ch, sst); end
                        else %for all other conditions check the first upcoming available channel
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
                %                 if(mod(p,100)==0)
                %                  fprintf('%f %f \n',p,toc);
                %                 end
            end  %for p = 1: npackets
            
            a111=1;
            %%markov approach%%%
            %             for qw=1:length(p)
            %      pp = p(qw);
            % transition_matrix = [pp,0,(1-pp);pp,0,(1-pp);0,1,0];
            % state_zero = [(v(qw)/(v(qw)+u(qw)))*pp;(u(qw)/(v(qw)+u(qw)));(v(qw)/(v(qw)+u(qw)))*(1-pp)];
            % current_state(:,qw) = (transition_matrix^(1))*state_zero
            % end
            % Scenario1 = (current_state(1,:)) .* s1; %arrive at available time slot and it has enough time to do service
            %             Scenario2 =(current_state(2,:)).* (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
            %             Scenario3 = (current_state(3,:)).* (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
            %               state1 = exp(-sv1/a)*a/(a+b);
            %               state3 = (1-exp(-sv1/a))*a/(a+b);
            %               state2 = b/(a+b);
            %            transition = [exp(-sv1/a) exp(-sv1/a) exp(-sv1/a);0 0 0;1-exp(-sv1/a) 1-exp(-sv1/a) 1-exp(-sv1/a)]^1e4;
            % %
            % %               transition = [exp(-sv1/a) 0 1-exp(-sv1/a);exp(-sv1/a) 0 1-exp(-sv1/a);0 1 0]^1e4
            %               state = transition*[state1;state2;state3];
            %%end markov approach%%%
            %%%Debug part%%%
            if(debug_f_q == 1)
                s1=sv1;
                u = b;
                v = 1/a;
                p              = exp(-s1*v);
                v1             = ((((-v*s1-1)*exp(-v*s1))+1)/v)/(1-exp(-s1*v));
                E_uv_2         = 2.*u.^2+((-(s1.^2+2.*s1.*v+2.*v.^2).*exp(-s1./v)+2.*v.^2)./(1-exp(-s1./v)))+2.*u.*v1;
                
                ES_BusyArrive  = (((1./p)-1).*((u)+(v1)))/mean(1-PER);
                ES_AvArrive    = ((1./p).*(u))+((((1./p)-1)).*(v1));
                mean_first_u=mean_first_u(mean_first_u~=0);
                mean_first_v = mean_first_v(mean_first_v~=0);
                %                 first_vb=first_vb(first_vb~=0);
                %                 first_ub=first_ub(first_ub~=0);
                %                 fprintf('fu:%f fvb:%f fub:%f, %f %f \n',mean(mean_first_u),mean(first_vb),mean(first_ub),a,b);
                
                % pfirsttbigger = exp(-sv1*lambda);
                % %psecondttoobig = exp((2*sv1-((((-lambda*sv1-1)*exp(-lambda*sv1))+1)/lambda)/(1-exp(-sv1*lambda)))/(-lambda))*(1-exp(-sv1/lambda));
                % psecondttoobig = (1-exp(-2*lambda*sv1))/(2*lambda);
                % zerostat = [1-pfirsttbigger;pfirsttbigger];
                % transmat = [1-(psecondttoobig) psecondttoobig;1-pfirsttbigger pfirsttbigger];
                % stat = transmat*zerostat;
                % pnonmarkov = pfirsttbigger + psecondttoobig;
                % % stat(1)*(state(2)) + (stat(2))*(b/(a+b))
                % % stat(1)*(state(1)) + (stat(2))*((exp(-sv1/a)* a/(a+b)))
                % % stat(1)*(state(3)) + (stat(2))*(((1-exp(-sv1/a))* a/(a+b)))
                % pnonmarkov = pfirsttbigger + psecondttoobig;
                % pnonmarkov*(state(2)) + (1-pnonmarkov)*(b/(a+b));
                % pnonmarkov*(state(1)) + (1-pnonmarkov)*((exp(-sv1/a)* a/(a+b)));
                % pnonmarkov*(state(3)) + (1-pnonmarkov)*(((1-exp(-sv1/a))* a/(a+b)));
                a1 = -exp(-sv1/a)*(a*(a+sv1)/(a*(1-exp(-sv1/a))))+a^2/(a*(1-exp(-sv1/a)));
                if(sim.trafficmode == 1)
                fprintf('portion of arrive at busy:%f %f arrive at available:%f %f arrive at insufficient available:%f %f\n',P2/npackets,0,P3/npackets,exp(-sv1/a),P1/npackets,((1-exp(-sv1/a))));%,b/(a+b),(exp(-sv1/a)* a/(a+b)),(exp(-sv1/a)* a/(a+b)));
                else
                    fprintf('portion of arrive at busy:%f %f arrive at available:%f %f arrive at insufficient available:%f %f\n',P2/npackets,b/(a+b),P3/npackets,a/(a+b)*exp(-sv1/a),P1/npackets,a/(a+b)*((1-exp(-sv1/a))));%,b/(a+b),(exp(-sv1/a)* a/(a+b)),(exp(-sv1/a)* a/(a+b)));
                end
                fprintf('average number of slots when arrive at busy:%f %f arrive at insufficient available:%f %f\n',mean(temp(first_slot4~=0)),(1./p)-1,mean(temp(first_slot2~=0)),(2*((exp(sv1/a)-1)))/mean(1-PER));
                %    fprintf('portion of arrive at busy:%f %f arrive at available:%f %f arrive at insufficient available:%f %f\n',P2/npackets,state(3),P3/npackets,state(2),P1/npackets,state(3));%,b/(a+b),(exp(-sv1/a)* a/(a+b)),(exp(-sv1/a)* a/(a+b)));
                mean_u_v_available=mean_u_v_available(first_slot2~=0);
                temp=temp(temp~=0);
                fprintf('mean u+v available: %f %f %f\n',mean(mean_u_v_available),ES_BusyArrive);
                mean_u_v_busy=mean_u_v_busy(mean_u_v_busy~=0);
                fprintf('mean u+v busy: %f %f\n',mean(mean_u_v_busy),ES_BusyArrive);
                opt1 = find(status_flag==2);
                fprintf('Number of(u+v)when arrive at busy time slot:%f  Analytical:%f \n',mean(P4(opt1)),(1/exp(-sv1/a))-1);
                opt2 = find(status_flag==3);
                fprintf('Number of (u+v)when arrive at small available time slot:%f  Analytical:%f\n',mean(P(opt2)),(1/exp(-sv1/a))-1);
                last_slot2 = last_slot;
                ww = find(first_slot~=0);
                last_slot2 = last_slot2(ww);
                qq = find(first_slot2~=0);
                last_slot = last_slot(qq);
                first_slot = first_slot(find(first_slot~=0));
                first_slot2 = first_slot2(find(first_slot2~=0));
                ES_AvArrive = ((1/exp(-sv1/a))-1)*(a1+b);
                fprintf('mean of service time of busy time   :%f  Analytical:%f\n', mean(cht(last_slot2)-cht(first_slot)),ES_AvArrive)
                % fprintf('mean of service time of small available time slots from end of end of first busy to last moment:%f  Analytical:%f\n', mean(mean_u_v_available/temp),ES_AvArrive+b)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%New Approximation%%%%%%%%%%%%%%%%%%%%%%%%%%
                % s1 = sv1;
                % v = a;
                % u =b;
                % p = exp(-s1./v);
                % v1             = -exp(-s1./v).*(v.*(v+s1)./(v.*(1-exp(-s1./v))))+v.^2./(v.*(1-exp(-s1./v)));
                % E_uv_2         = 2.*u.^2+((-(s1.^2+2.*s1.*v+2.*v.^2).*exp(-s1./v)+2.*v.^2)./(1-exp(-s1./v)))+2.*u.*v1;
                %
                % ES_BusyArrive  = ((1./p)-1).*((u)+(v1));
                % ES_AvArrive    = ((1./p).*(u))+((((1./p)-1)).*(v1));
                % Scenario1 = s1; %arrive at available time slot and it has enough time to do service
                % Scenario2 = (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
                % Scenario3 = (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
                % Pa = p;
                % Pc = (1-p);
                % P1 = Pa*(v/(u+v));
                % P2 = (u/(u+v));
                % P3 = Pc*(v/(u+v));
                % T = 1/lambda;
                % transition = [(1-exp(-lambda*Scenario1))*Pa,0,(1-exp(-lambda*Scenario1))*Pc,(exp(-lambda*Scenario1))*P1,(exp(-lambda*Scenario1))*P2,(exp(-lambda*Scenario1))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*Scenario2))*Pa,0,(1-exp(-lambda*Scenario2))*Pc,0,0,0,(exp(-lambda*Scenario2))*P1,(exp(-lambda*Scenario2))*P2,(exp(-lambda*Scenario2))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*Scenario3))*Pa,0,(1-exp(-lambda*Scenario3))*Pc,0,0,0,0,0,0,(exp(-lambda*Scenario3))*P1,(exp(-lambda*Scenario3))*P2,(exp(-lambda*Scenario3))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*2*Scenario1))*Pa,0,(1-exp(-lambda*2*Scenario1))*Pc,0,0,0,0,0,0,0,0,0,(exp(-lambda*2*Scenario1))*P1,(exp(-lambda*2*Scenario1))*P2,(exp(-lambda*2*Scenario1))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*(Scenario1+Scenario2)))*Pa,0,(1-exp(-lambda*(Scenario1+Scenario2)))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*(Scenario1+Scenario2)))*P1,(exp(-lambda*(Scenario1+Scenario2)))*P2,(exp(-lambda*(Scenario1+Scenario2)))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*(Scenario1+Scenario3)))*Pa,0,(1-exp(-lambda*(Scenario1+Scenario3)))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*(Scenario1+Scenario3)))*P1,(exp(-lambda*(Scenario1+Scenario3)))*P2,(exp(-lambda*(Scenario1+Scenario3)))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*(Scenario2+Scenario1)))*Pa,0,(1-exp(-lambda*(Scenario2+Scenario1)))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*(Scenario2+Scenario1)))*P1,(exp(-lambda*(Scenario2+Scenario1)))*P2,(exp(-lambda*(Scenario2+Scenario1)))*P3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*2*Scenario2))*Pa,0,(1-exp(-lambda*2*Scenario2))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*2*Scenario2))*P1,(exp(-lambda*2*Scenario2))*P2,(exp(-lambda*2*Scenario2))*P3,0,0,0,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*(Scenario2+Scenario3)))*Pa,0,(1-exp(-lambda*(Scenario2+Scenario3)))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*(Scenario2+Scenario3)))*P1,(exp(-lambda*(Scenario2+Scenario3)))*P2,(exp(-lambda*(Scenario2+Scenario3)))*P3,0,0,0,0,0,0,0,0,0;
                %                (1-exp(-lambda*(Scenario3+Scenario1)))*Pa,0,(1-exp(-lambda*(Scenario3+Scenario1)))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*(Scenario3+Scenario1)))*P1,(exp(-lambda*(Scenario3+Scenario1)))*P2,(exp(-lambda*(Scenario3+Scenario1)))*P3,0,0,0,0,0,0;
                %                (1-exp(-lambda*(Scenario3+Scenario2)))*Pa,0,(1-exp(-lambda*(Scenario3+Scenario2)))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*(Scenario3+Scenario2)))*P1,(exp(-lambda*(Scenario3+Scenario2)))*P2,(exp(-lambda*(Scenario3+Scenario2)))*P3,0,0,0;
                %                (1-exp(-lambda*2*Scenario3))*Pa,0,(1-exp(-lambda*2*Scenario3))*Pc,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(exp(-lambda*2*Scenario3))*P1,(exp(-lambda*2*Scenario3))*P2,(exp(-lambda*2*Scenario3))*P3
                %                ;zeros(27,39)];
                
                % transition = [(1-exp(-lambda*Scenario1))*Pa,0,(1-exp(-lambda*Scenario1))*Pc,(exp(-lambda*Scenario1))*P1,(exp(-lambda*Scenario1))*P2,(exp(-lambda*Scenario1))*P3,0,0,0,0,0,0;
                %                (1-exp(-lambda*Scenario2))*Pa,0,(1-exp(-lambda*Scenario2))*Pc,0,0,0,(exp(-lambda*Scenario2))*P1,(exp(-lambda*Scenario2))*P2,(exp(-lambda*Scenario2))*P3,0,0,0;
                %                (1-exp(-lambda*Scenario3))*Pa,0,(1-exp(-lambda*Scenario3))*Pc,0,0,0,0,0,0,(exp(-lambda*Scenario3))*P1,(exp(-lambda*Scenario3))*P2,(exp(-lambda*Scenario3))*P3;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %                Pa,0,Pc,0,0,0,0,0,0,0,0,0;
                %              ];
                % states_zero = [P1,P2,P3,zeros(1,9)];
                % states = states_zero*(transition^1e3);
                % waiting_matrix = [0;0;0;Scenario1-T;Scenario1-T;Scenario1-T;Scenario2-T;Scenario2-T;Scenario2-T;Scenario3-T;Scenario3-T;Scenario3-T;zeros(27,1)];
                % Total_wait = waiting_matrix *states*1000;
                % Average_wait = sum(Total_wait(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%EnD of Debug part%%%
            fprintf('mean of phi   :%f  Analytical:%f\n', mean(mean_first_v),a1);
            CRNwnv=CRNwnv(1:end-1);
            queue.CRNwnv = CRNwnv;  queue.CRNsv = CRNsv;  queue.cht = cht;  queue.chf = chf;
        end
    end
end

