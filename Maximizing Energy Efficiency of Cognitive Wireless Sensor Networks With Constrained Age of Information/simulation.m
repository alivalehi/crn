clear variables;

check_each_part = 1;
check_p         = 1;

Av_mean   = 3; 
Bu_mean   = 6;
s1        = 3;
n         = 1e4;
k         = 100;
L         = 1e3;
av_time   = exprnd(Av_mean,n,L);
bu_time   = exprnd(Bu_mean,n,L);
channel   = reshape([av_time; bu_time], [n, 2*L]);
chf       = [ zeros(n,1),reshape([ones(n,L); zeros(n,L)], [n, 2*L])]; %0 means next slot is available  1 means next slot is busy 
slot      = nan(n,2*L);
sim       = nan(n,k);
slot(:,1) = 0;
symbols   = nan(1,L);
cht       = nan(1,L);
flag      = nan(1,L);
chs       = nan(1,L);
w         = nan(n,2*L); 
for i=1:2*L
    slot(:,i+1) = slot(:,i) + channel(:,i);
end
                
alpha        = rand(n,2*L);
beta         = rand(n,2*L);
slot_tiles   = repmat(slot(:,end),1,k);
rand_slot    = alpha.* (2*L);
for g=1:n  
    for h=1:2*L
        yyy   = channel(g,1+floor(rand_slot(g,h)));
        yy    = slot(g,1+floor(rand_slot(g,h)));
        yyyy  = beta(g,h);
        w(g,h)= yy+yyy*yyyy;
     end
end
status_flag  = zeros(n,k);
P            = zeros(1,k);
mean_av      = zeros(n,k);
mean_bu      = zeros(n,k);
mean_sig_bu  = zeros(n,k);
mean_sav     = zeros(n,k);
mean_sig_sav = zeros(n,k);
mean_P_sav   = zeros(n,k);
mean_P_bu    = zeros(n,k);
mean_sig_sav_total = zeros(1,k*n);
mean_sig_bu_total = zeros(1,k*n);
first_slot = zeros(1,k);
tic
for i=1:n 
    symbols = sort(w(i,:));
    cht     = slot(i,:);
    flag    = chf(i,:);
    chs     = [channel(i,:),0];
    
    for j=1:k    
        sst        = symbols(1,j);
        first_slot(j) = find(cht > sst,1,'first');
        send_flag  = 0;
        first_flag = 0;
        while (~send_flag )
            slot_s1 = find(cht <= sst,1,'last');
            slot_s2 = find(cht > sst,1,'first');
            slot_e1 = find(cht <= sst+s1,1,'last');
            slot_e2 = find(cht > sst+s1,1,'first');
            
            if isempty(slot_e2)  %extend channel availabity process and retry
                flag = [flag, flag(2:end)];
                cht  = [cht, max(cht)+ cht(2:end)];
                chs  = [channel(i,:),channel(i,:),0];

            elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (slot_s1 == slot_e1) && (flag(slot_s1) == 0)  %packet transmission fits into one availale time slot  
                if(~first_flag)
                    status_flag(i,j) = 1;
                    first_flag     = 1;
                end

                sst       = sst+ s1;  %set start time for next retransmission
                send_flag = 1;        %mark the packet as transmitted
                            
            elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (flag(slot_s1) == 1)  %packet transmission arrive at busy time slot 
                 if(~first_flag)
                    status_flag(i,j) = 2;
                    first_flag     = 1;                               
                 end
                 
                 next_av_ch = find((cht > sst) & (flag == 0)& (chs >= s1),1,'first');   %chf==0, means the time slot right after unavailable channel
                 if isempty(next_av_ch)  % extend the ch. avail process
                    flag = [flag, flag(2:end)];
                    cht  = [cht, max(cht)+ cht(2:end)];
                    chs  = [channel(i,:),channel(i,:),0];
                    next_av_ch = find((cht > sst) & (flag == 0)& (chs >= s1),1,'first'); 
                 end
                 P(i,j) = (next_av_ch-first_slot(j))/2;
                 sst = cht(next_av_ch);
                for o=1:P(i,j)
                    mean_uv_b(i,o) = cht(first_slot(j)+(2*o))-cht(first_slot(j)+(2*(o-1))); 
                    mean_vb(i,o) = cht(first_slot(j)+((2*o)-1))-cht(first_slot(j)+(2*(o-1)));
                    mean_ub(i,o) = cht(first_slot(j)+((2*o)))-cht(first_slot(j)+(2*(o-1))+1);
                    maen_first_u(i,o) = cht(first_slot(j))-symbols(1,j);
                end    
                 
                 if isempty(next_av_ch), queue.fail=1; end

            elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (slot_s1 ~= slot_e1)  && (flag(slot_s1) == 0)  %packet transmission does not fits into one availale time slot 
                 if(~first_flag)
                     status_flag(i,j)= 3 ;
                     first_flag = 1;
                    % first_slot(j) = slot_s2+1;
                 end
                 
                 next_av_ch = find((cht > sst) & (flag == 0)& (chs >= s1),1,'first');
                 if isempty(next_av_ch)  % extend the ch. avail process
                    flag = [flag, flag(2:end)];
                    cht  = [cht, max(cht)+ cht(2:end)];
                    chs  = [channel(i,:),channel(i,:),0];
                    next_av_ch = find((cht > sst) & (flag == 0)& (chs >= s1),1,'first'); 
                 end
                 P(i,j) = (next_av_ch-first_slot(j))/2;
                 sst = cht(next_av_ch);
                 for o=1:P(i,j)
                    mean_uv(i,o) = cht(first_slot(j)+(2*o))-cht(first_slot(j)+(2*(o-1))); 
                    mean_v(i,o) = cht(first_slot(j)+((2*o)-1))-cht(first_slot(j)+(2*(o-1)));
                    mean_u(i,o) = cht(first_slot(j)+((2*o)))-cht(first_slot(j)+(2*(o-1))+1);
                    maen_first_v(i,o) = cht(first_slot(j)-1)-symbols(1,j);
                 end 
                 if isempty(next_av_ch), queue.fail=1; end
            else
                pause;
            end
        end  %end while          while (~send_flag && ret_cont<100)
        sim(i,j) = sst-symbols(1,j);
        sim1(i,j) = sst-s1-cht(first_slot(j));
    end
    symbols = nan(1,n);
    cht     = nan(1,n);
    flag    = nan(1,n);
    chs     = nan(1,n);
    
    if(check_each_part)
        opt1 = find(status_flag(i,:)==1);
        if (~isempty(opt1))
            mean_av(i)      = mean(sim(i,opt1));
        end
        opt2 = find(status_flag(i,:)==2);
        if (~isempty(opt2))
            mean_bu(i)         = mean(sim(i,opt2));
            mean_bu2(i)        = mean(sim(i,opt2).^2);
            mean_sig_bu(i)     = mean(sim1(i,opt2));
            mean_sig_bu2(i)    = mean((sim1(i,opt2)).^2);
            var_sig_bu(i)      = var(sim1(i,opt2));
        end
        opt3 = find(status_flag(i,:)==3);
        if (~isempty(opt3))
            mean_sav(i)     = mean(sim(i,opt3));
            mean_sav2(i)     = mean(sim(i,opt3).^2);
            mean_sig_sav(i) = mean(sim1(i,opt3));
            mean_sig_sav2(i)   = mean((sim1(i,opt3)).^2);
            var_sig_sav(i)     = var(sim1(i,opt3));
            
        end
    end
%     sav_length = find(mean_sig_sav_total~=0,1,'last');
%     b_length = find(mean_sig_bu_total~=0,1,'last');
%     if(isempty(sav_length))
%         sav_length = 0;
%     end    
%     if(isempty(b_length))
%         b_length = 0;
%     end    
%     mean_sig_sav_total(1,sav_length+1:sav_length+length(mean_sig_sav)) = mean_sig_sav(:);
%     mean_sig_bu_total(1,b_length+1:b_length+length(mean_sig_bu)) = mean_sig_bu(:);
   if(mod(i,1000)==0)
    fprintf('%f %f \n',i,toc);
   end 
    % clearvars symbols cht flag chs
end
mean_sig_sav_total = mean_sig_sav_total(2:end);
mean_sig_bu_total  = mean_sig_bu_total(2:end);

    if(check_p)
        opt4 = find(status_flag(:,1)==3);
        if (~isempty(opt4))
            mean_P_sav(i) = mean(P(opt4));
        end
        opt5 = find(status_flag(:,1)==2);
        if (~isempty(opt5))
            mean_P_bu(i)  = mean(P(opt5));
        end
    end
            u              = Bu_mean;
            v              = Av_mean;
            p              = exp(-s1/v);
            v1             = -exp(-s1/v)*(v*(v+s1)/(v*(1-exp(-s1/v))))+v^2/(v*(1-exp(-s1/v)));
            E_uv_2         = 2*u^2+((-(s1^2+2*s1*v+2*v^2)*exp(-s1/v)+2*v^2)/(1-exp(-s1/v)))+2*u*v1;
            a              = 1;
            ES_BusyArrive  = ((1/p)-a)*((u)+(v1));
            ES_AvArrive    = ((1/p)*(u))+((((1/p)-a))*(v1));
            E_bu2          = E_uv_2 *((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v1)^2;
            ES_BusyArrive2 = 2*u.^2/3 + E_bu2+s1^2 +2*s1*(u/2+ES_BusyArrive)+ 2*ES_BusyArrive*u/2;
            E_sav2         = E_bu2+(2*u^2)+2*u*((1/p)-a)*(u+v1);
            var_bu         = E_bu2-ES_BusyArrive^2;
            var_sav        = E_sav2-ES_BusyArrive^2; % it should be (E_sav2- ES_AvArrive^2);
            fprintf('mean_sav :%f %f ',mean(mean_sav(mean_sav~=0)),ES_AvArrive+s1+v/2);
            fprintf('mean_bu :%f %f \n',mean(mean_bu(mean_bu~=0)),ES_BusyArrive+s1+u/2);
            fprintf('mean_sav2 :%f %f ',mean(mean_sav2(mean_sav2~=0)),ES_AvArrive+s1+v/2);
            fprintf('mean_bu2 :%f %f \n',mean(mean_bu2(mean_bu2~=0)),ES_BusyArrive2);
            fprintf('mean_sig_sav2 :%f %f ',mean((mean_sig_sav2(mean_sig_sav2~=0))),E_sav2);
            fprintf('mean _sig_bu2 :%f %f \n',mean((mean_sig_bu2(mean_sig_bu2~=0))),E_bu2);
            fprintf('var_sig_sav :%f %f ',mean(var_sig_sav(var_sig_sav~=0)),var_sav);
            fprintf('var_sig_bu :%f %f \n',mean(var_sig_bu(var_sig_bu~=0)),var_bu);
            fprintf('mean_sig_sav :%f %f ',mean(mean_sig_sav(mean_sig_sav~=0)),ES_AvArrive);
            fprintf('mean_sig_bu :%f %f \n',mean(mean_sig_bu(mean_sig_bu~=0)),ES_BusyArrive);
            fprintf('mean_P_sav :%f %f ',mean(mean_P_sav(mean_P_sav~=0)),((1/p)-a));
            fprintf('mean_P_bu :%f %f \n',mean(mean_P_bu(mean_P_bu~=0)),((1/p)-a));
            fprintf('mean_uv_b :%f mean_uv: %f \n',mean(mean_uv_b(mean_uv_b~=0)),mean(mean_uv(mean_uv~=0)));
            fprintf('mean_u :%f mean_v: %f \n',mean(mean_u(mean_u~=0)),mean(mean_v(mean_v~=0)));
            fprintf('mean_ub :%f mean_vb: %f \n',mean(mean_ub(mean_ub~=0)),mean(mean_vb(mean_vb~=0)));
            fprintf('first_mean_b :%f first_mean_v: %f \n',mean(maen_first_u(maen_first_u~=0)),mean(maen_first_v(maen_first_v~=0)));
            FileName = ['ali''s result\proj',datestr(now, 'mmmm dd yyyy HH-MM-SS'),'.mat'];
       %     save(FileName); 
        %    pause;