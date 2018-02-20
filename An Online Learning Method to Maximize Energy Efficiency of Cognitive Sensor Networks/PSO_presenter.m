function [aa,NN1,NN2,NN3,EE_real,t_t,sim] = EE_simulator(k,sim,t_t,n)
kv= k;
ali = [];
P1_total= [];
P2_total= [];
P3_total= [];
slot_tracker =[];
sim.startslot = 0;
    sim.endslot = 0;
%sim = channel_perform(sim);
length(kv);
%tic
for j=1: length(kv)
 %   if(mod(j,100)==0)
       % j
      %  toc
  %  end
    k= kv(j);
    tp = t_t(k);
    if(isfield(sim,'lastsent'))
        if (tp<sim.lastsent)
            tp = sim.lastsent;
        end
    end
    t_t = t_t(k+1:end);
    if(length(t_t)<k)
        aa(j) = 0;
        NN(j) = 0;
        EE_real(j)=0;
        fprintf('return type1 from simuklation');
        return
    end
    sim.control.debug_active=0;
    lambda = sim.lambda;
    debug_f_q = 1;
    chf = sim.chf;
    cht = sim.cht;
   vtrack=sim.vtrack;
utrack=sim.utrack;
    BERch = sim.BERch;
    CRNwnv = []; CRNsv =[];
    npackets=1;
    p=1;
    CRNwnv = zeros(1,npackets); CRNsv = zeros(1,npackets);
    CRNwnv(1) = 0;
    %save temp1, disp('temp1 saved'); pause(2);
    %%%Debug part%%%
    if(debug_f_q == 1)
        N1 = 0;% # of insufficient available arrive
        N2 = 0;% # of busy arrive
        N3 = 0;% # of  available  arrive
        P3 = 0;% # of  available  arrive
        P1 = 0;% # of insufficient available arrive
        P2 = 0;% # of busy arrive
        
    end
    %%%EnD of Debug part%%%
    
    
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
    rv_flag= 1;
    r_flag=0;
    rv= 1e5;
    Plen = k * sim.N + sim.H;
    sv1 = Plen/sim.ch.Rch;
    r=1;
    timespan1(1)=0;
    timespan2(1)=0;
    timespan3(1)=0;
    timespan4(1)=0;
    timespan5(1)=0;
    
    while (~r_flag)   %try for all retransmission
        
        send_flag = 0; ret_cont=0;first_flag=0;
        while (~send_flag )
            ret_cont=ret_cont+1;
            cht_s1 = max(find(cht <= sst)); cht_s2 = min(find(cht > sst));
            cht_e1 = max(find(cht <= sst+sv1)); cht_e2 = min(find(cht > sst+sv1));
           % fprintf(' u:%f v:%f simulation \n',vtrack(cht_s1),utrack(cht_s1));
           if(sim.startslot == 0)
               sim.startslot = cht_s1;
           end
           if(rv_flag)
                
                PER = 1 - (1 - BERch(cht_e1)) .^ (Plen);
                rv = 1+geornd(1-PER,1);     %retransmission in ARQ
                if(isnan(rv))
                    rv=1;
                end
                rv_flag = 0;
                
            end
            if(rv>100)
            fprintf('more than 100 retransmission \n');
            pause;
            end    
            if sim.control.debug_active, fprintf('p:%d r:%d retry:%d  arrival:%1.4f wt:%1.4f   startp:%1.4f  pt[%1.4f-%1.4f], start(%d-%d) end(%d-%d) \n', p,r,ret_cont, tp(p), CRNwnv(p), sst0, sst, sst+sv1, cht_s1, cht_s2, cht_e1, cht_e2); end
            
            if isempty(cht_e2)  %extend channel availabity process and retry
                chf = [chf, chf(2:end)];
                cht = [cht, max(cht)+ cht(2:end)];
    utrack=[utrack,utrack(2:end)];
    vtrack=[vtrack,vtrack(2:end)];
    sim.slotuv = [sim.slotuv, reshape([vtrack; utrack], [1, 2*length(vtrack)])];
fprintf(' Packet number %d  \n',j);
fprintf(' retransmit number %d  \n',r);
                fprintf('  ch extended \n');

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %packet transmission fits into one availale time slot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 == cht_e1) && (chf(cht_e2) == 1)
                sst = sst+ sv1;  %set start time for next retransmission
                sim.lastsent = sst;
                N3(r) = N3(r)+1;
                send_flag = 1;   %mark the packet as transmitted
                if sim.control.debug_active, fprintf('  OK \n'); end
                if(~first_flag)
                    first_flag = 1;
                    P3(r)= P3(r)+1;
                end
                if(~r_flag)
                    timespan3(end+1) = cht(cht_e2)-sst;
                end
                
                timespan3(end+1) = cht(cht_e2)-cht(cht_e1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %packet transmission arrive at busy time slot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (chf(cht_s1) == 1)
                N2(r) = N2(r)+1;
                next_av_ch = min(find((cht > sst) & (chf == 0)));  %chf==0, means the time slot right after unavailable channel
                if isempty(next_av_ch)  % extend the ch. avail process
                    chf = [chf, chf(2:end)];
                    cht = [cht, max(cht)+ cht(2:end)];
                    next_av_ch = min(find((cht > sst) & (chf == 0)));
                end
                sst = cht(next_av_ch);
                if sim.control.debug_active, fprintf('  Busy Retry at time(%d):%1.4f  \n',next_av_ch, sst); end
                if isempty(next_av_ch), queue.fail=1; end
                
                if(~first_flag)
                    first_flag = 1;
                    P2(r)= P2(r)+1;
                end
                
                timespan2(end+1) = cht(cht_s2)-cht(cht_s1);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %packet transmission does not fits into one availale time slot
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif (cht_s1+1 == cht_s2) && (cht_e1+1 == cht_e2) && (cht_s1 ~= cht_e1)  && (chf(cht_s1) == 0)
                N1(r) = N1(r)+1;
                timespan1(end+1) = cht(cht_s2)-sst;
                next_av_ch = min(find((cht > sst) & (chf == 0)));  %chf==0, means the time slot right after unavailable channel
                if isempty(next_av_ch)  % extend the ch. avail process
                    chf = [chf, chf(2:end)];
                    cht = [cht, max(cht)+ cht(2:end)];
                    next_av_ch = min(find((cht > sst) & (chf == 0)));
                end
                sst = cht(next_av_ch);
                if sim.control.debug_active, fprintf('  s av Retry at time(%d):%1.4f  \n',next_av_ch, sst); end
                if isempty(next_av_ch), queue.fail=1; end
                
                if(~first_flag)
                    first_flag = 1;
                    P1(r)= P1(r)+1;
                end
                
                timespan4(end+1) = cht(cht_s2)-cht(cht_s1);
                
                
            else %for all other conditions check the first upcoming available channel
                next_av_ch = min(find((cht > sst) & (chf == 0)));  %chf==0, means the time slot right after unavailable channel
                if isempty(next_av_ch)  % extend the ch. avail process
                    chf = [chf, chf(2:end)];
                    cht = [cht, max(cht)+ cht(2:end)];
                    next_av_ch = min(find((cht > sst) & (chf == 0)));
                end
                sst = cht(next_av_ch);
                if sim.control.debug_active, fprintf(' unkbown  Retry at time(%d):%1.4f  \n',next_av_ch, sst); end
                if isempty(next_av_ch), queue.fail=1; end
            end
            
        end  %end while
        
        
        if(r==rv)
            r_flag=1;
            
        end
        r = r+1;
        N1(r)=0;
        N2(r)=0;
        N3(r)=0;
        P1(r)=0;
        P2(r)=0;
        P3(r)=0;
        
    end %end for r = 1 : rv(p)
    slot_tracker(end+1) = cht_e1;
    sim.endslot = cht_e1;
%     if(cht_s1>=sim.Nparamset*sim.nch*2)
%         sim.Nparamset =sim.Nparamset+1;
%         sim.channel_change_flag =1;
%         aa(j) = 0;
%         NN1(j) = 0;
%         NN2(j) = 0;
%         NN3(j) = 0;
%         EE_real(j)=0;
%         return
%     end
    cht_s1;
    if(isnan(mean(P1(P1~=0))))
        mean_P1 = 0;
    else
        mean_P1 =  mean(P1(P1~=0));
    end
    
    if(isnan(mean(P2(P2~=0))))
        mean_P2 = 0;
    else
        mean_P2 =  mean(P2(P2~=0));
    end
    
    if(isnan(mean(P3(P3~=0))))
        mean_P3 = 0;
    else
        mean_P3 =  mean(P3(P3~=0));
    end
    
    P1_total=  [P1_total sum(P1)];
    P2_total=  [P2_total sum(P2)];
    P3_total=  [P3_total sum(P3)];
    
    %     P1_total=  [P1_total mean_P1];
    %     P2_total=  [P2_total mean_P2];
    %     P3_total=  [P3_total mean_P3];
    N1_mean = mean(N1(N1~=0));
    N1_sum = sum(N1);
    NNN1(j) = mean(N1);
    N2_mean = mean(N2(N2~=0));
    N2_sum = sum(N2);
    
    N3_mean = mean(N3(N3~=0));
    N3_sum = sum(N3);
    
    N_mean = mean(N1_sum);
    %     ali1= N_sum;
    %     ali2 =sum(N2);
    %     ali = [ali ali1];
    if(isnan(N1_mean))
        N1_mean = 0;
    end
    if(length(kv)<=1)
        %  EE_real  = (((1+rv)*Plen + ((Plen/2*N1_mean)))*sim.power_bit+sim.power_sense)/(k*sim.N);
        %EE_real  = ((Plen*(N3_sum+(N1_sum/2)))*sim.power_bit+sim.power_sense)/(k*sim.N);
        EE_real  = (sim.power_bit*(N1_sum*(Plen/2)+Plen*N3_sum)+sim.power_sense)/(k*sim.N);% changed on 09/16/2017
        %aa= N_mean;
        aa(j)= rv;
        NN = mean(ali);
        NN1(j) = N1_sum; %small available
        NN2(j) = N2_sum; %busy
        NN3(j) = N3_sum;%available
    else
        %EE_real(j)  = ((Plen*(N3_sum+(N1_sum/2)))*sim.power_bit+sim.power_sense)/(k*sim.N);
        % EE_real(j)  = (((N1_mean/2)+rv).*Plen.*sim.power_bit+sim.power_sense)/(k*sim.N);
       EE_real(j)  = (sim.power_bit*(N1_sum*(Plen/2)+Plen*N3_sum)+sim.power_sense)/(k*sim.N);% changed on 09/16/2017
        aa(j)= rv;
        NN1(j) = N1_sum; %small available
        NN2(j) = N2_sum; %busy
        NN3(j) = N3_sum;%available
    end
    N1 = 0;% # of insufficient available arrive
    N2 = 0;% # of busy arrive
    N3 = 0;% # of  available  arrive
end
%fprintf('Av : %f SAV : %f B : %f\n',mean(NN3),mean(NN2),mean(NN1));
%aa;
NN = mean(ali);
sim.slottrack  = [sim.slottrack slot_tracker];
  sim.vtrack=vtrack;
sim.utrack=utrack;
mean(NNN1);
coeff11s = (N1_sum*(Plen/2)+Plen*N3_sum);
end
