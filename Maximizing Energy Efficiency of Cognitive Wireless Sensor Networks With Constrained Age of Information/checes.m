
if 1 %this part randomly generate one  packet and verify service time of that packet
            clear all;
            Av_mean = 2; 
            Bu_mean = 4;
            s1= 3;
            n=1e3;
            lambda = 30;
            u = Bu_mean;
            v = Av_mean;
            check_each_part = 1;
            check_p = 1;
            p  = exp(-s1/v);
            a= 1;
            ES_BusyArrive = ((1/p)-a)*((u)+(v));
            ES_AvArrive = ((1/p)*(u))+((((1/p)-a))*(v));
            ES_BusyArrive2 = 2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2;
            ES_AvArrive2 = (2*(v^2+u^2+v*u)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(v+u)^2)-(u*(1/p-a)+v*(1/p-a))^2+2*u^2+ES_AvArrive^2;
            var_BusyArrive = ((1/p)^2)*((u^2)+(v^2));
            var_AvArrive = ((1/p^2)*(u^2))+(((1/p^2)+1)*(v^2));
            Scenario1 = (p* v/(v+u)) .* s1;%available 
            Scenario2 = (u/(v+u)).* (s1 + u/2+ ES_BusyArrive);%busy
            Scenario3 = ((v/(v+u)).*(1-p)).* (ES_AvArrive+s1+s1/2);%insufficient available
            ESk = Scenario1 +Scenario2+Scenario3;
            ESk2 = ((p* v/(v+u)) .* s1^2)+((u./(v+u)).* (s1.^2 + 2*u.^2/3 +ES_BusyArrive2 + s1.*u + 2*s1*ES_BusyArrive +u*ES_BusyArrive))+((v./(v+u)).*(1-p)).* (s1.^2 + 2*v.^2/3 + s1.*v+ES_AvArrive2+2*s1*ES_AvArrive+v*ES_AvArrive);
            varsk2 = ESk2- ESk.^2;  
  for j=1:n 
            av_time = exprnd(Av_mean,1,n);
            bu_time = exprnd(Bu_mean,1,n);
            channel = reshape([av_time; bu_time], [1, 2*n]);
            chf = reshape([ones(1,n); zeros(1,n)], [1, 2*n]); %0 means next slot is available  1 means next slot is busy 
            slot = nan(size(channel));
          %  w = nan(n);
            slot(1) = 0;
           
            for i=1:length(slot)
                slot(j,i+1) = slot(i) + channel(i);
            end
               
            chf = [0,chf];
            alpha = rand(1,n);
            finalslot = slot(j,end);
            w(j,:) = alpha.* finalslot;
            w = sort(w);
   
            status_flag = zeros(1,n);
            sim1 = zeros(1,n);
            P = zeros(1,n); % # of insufficient available when arrive at insufficient available
            P3 = 0;% # of  available  arrive
            P1 = 0;% # of insufficient available arrive 
            P2 = 0;% # of busy arrive
            P4 = zeros(1,n);% # of insufficient available when arrive at busy
            mean_first_u = zeros(1,n);%just for check mean of first busy slot in insufficient available time slot
            mean_u_v_available = zeros(1,n);%just for check mean of first busy slot in insufficient available time slot
            mean_u_v_busy = zeros(1,n);
            temp = zeros(1,n);
            tic
            for i=1:n 
                sst = w(j,i);
                send_flag = 0;
                tag = 'Available'; 
                first_flag = 0;
                long_sav_flag = 0;
                long_sav_flag2 = 0;
                while (~send_flag )
                        slot_s1 = find(slot(j,:) <= sst,1,'last'); slot_s2 = find(slot(j,:) > sst,1,'first');
                        slot_e1 = find(slot(j,:) < sst+s1,1,'last'); slot_e2 = find(slot(j,:) >= sst+s1,1,'first');
                      

                        if isempty(slot_e2)  %extend channel availabity process and retry
                             chf = [chf, chf(2:end)];
                            slot = [slot, max(slot)+ slot(2:end)];

%                              chf = [chf,1,0];
%                              bu_new_exp = exprnd(Bu_mean);
%                              av_new_exp = exprnd(Av_mean);
%                              slot = [slot, max(slot)+bu_new_exp,max(slot)+bu_new_exp+av_new_exp];
                         elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (slot_s1 == slot_e1) && (chf(slot_e2) == 1)  %packet transmission fits into one availale time slot  
                            if(~first_flag)
                                status_flag(i)= 1 ;
                                first_flag = 1;
                                P3= P3+1;
                            end
                            sst = sst+ s1;   %set start time for next retransmission
                            send_flag = 1;   %mark the packet as transmitted
%                             last_slot(i)= slot_s1;
%                             if(status_flag(i)== 3)
%                                 if(P(i)~=0)
%                                     mean_u_v_available(i) = (slot(last_slot(i))-slot(first_slot2(i)))/(P(i));
%                                     temp(i) = slot(last_slot(i))-slot(first_slot2(i));
%                                 end
%                             end
%                              if(status_flag(i)== 2)
%                                 if(P4(i)~=0)
%                                     mean_u_v_busy(i) = (slot(last_slot(i))-slot(first_slot3(i)))/(P4(i));
%                                 end
%                              end
                             
                        elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (chf(slot_s1) == 1)  %packet transmission arrive at busy time slot 
                             if(~first_flag)
                                status_flag(i)= 2 ;
                                first_flag = 1;                               
                                tag = 'Busy'; 
                                P2= P2+1;
                             end   
                            next_av_ch = find((slot(j,:) > sst) & (chf == 0),1,'first');  %chf==0, means the time slot right after unavailable channel
                          
                            if isempty(next_av_ch)  % extend the ch. avail process
                                chf = [chf, chf(2:end)];
                                slot = [slot, max(slot)+ slot(2:end)];
                                next_av_ch = find((slot > sst) & (chf == 0),1,'first'); 
                            end
                            sst = slot(next_av_ch);
                            if isempty(next_av_ch), queue.fail=1; end

                        elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (slot_s1 ~= slot_e1)  && (chf(slot_s1) == 0)  %packet transmission does not fits into one availale time slot 
                             if(~first_flag)
                                status_flag(i)= 3 ;
                                first_flag = 1;
                                sim3(i) = slot(slot_s2)-w(i);
                                tag = 'Small Available';                                                          
                                first_slot2(i)= slot_s2+1;
                                mean_first_u(i) = slot(slot_s2+1)-slot(slot_s2);
                                P1= P1+1;
                             else 
                                 if(  status_flag(i)== 3)
                                     if(long_sav_flag==0)
                                        first_slot(i)= slot_s1-1;
                                        long_sav_flag = 1;
                                     end   
                                     P(i)= P(i)+1;
                                 else
                                     
                                     if(long_sav_flag2==0)
                                        first_slot3(i)= slot_s1;
                                        long_sav_flag2 = 1;
                                     end
                                     P4(i)= P4(i)+1;
                                 end
                             end
                            next_av_ch = find((slot > sst) & (chf == 0),1,'first');  %chf==0, means the time slot right after unavailable channel
                            if isempty(next_av_ch)  % extend the ch. avail process
                                chf = [chf, chf(2:end)];
                                slot = [slot, max(slot)+ slot(2:end)];
                                next_av_ch = find((slot > sst) & (chf == 0),1,'first'); 
                            end
                            sst = slot(next_av_ch); 
                            if isempty(next_av_ch), queue.fail=1; end
                    
                        else
                            disp('undefined');          
                        end
                 end  %end while
                    
                    sim1(i) = sst-w(i);

                  if(mod(i,200)==0)
                    i
                    toc
                    if(mod(i,20000)==0)
                       FileName=['ali''s result\proj',datestr(now, 'mmmm dd, yyyy HH-MM-SS'),'.mat'];
                       save(FileName);
                    end    
                  end             
            end
  end
                FileName=['ali''s result\proj',datestr(now, 'mmmm dd, yyyy HH-MM-SS'),'.mat'];
                save(FileName);
                fprintf('portion of arrive at busy:%f arrive at available:%f arrive at insufficient available:%f analytical : %f %f %f \n',P2/n,P3/n,P1/n,u/(v+u),(p* v/(v+u)),((v/(v+u)).*(1-p)));
                L2=sim3(sim3~=0);
                opt1 = find(status_flag==2);
                fprintf('Number of(u+v)when arrive at busy time slot:%f  Analytical:%f \n',mean(P4(opt1)),(1/p)-a);
                opt2 = find(status_flag==3);
                fprintf('Number of (u+v)when arrive at small available time slot:%f  Analytical:%f\n',mean(P(opt2)),(1/p)-a);    
                fprintf('scenario %f %f %f\n',Scenario1,Scenario2,Scenario3);
                mean_u_v_available=mean_u_v_available(mean_u_v_available~=0);
                fprintf('mean u+v available: %f %f\n',mean(mean_u_v_available),p*(u+v));
                mean_u_v_busy=mean_u_v_busy(mean_u_v_busy~=0);
                fprintf('mean u+v busy: %f %f\n',mean(mean_u_v_busy),p*(u+v));
                mean_first_u=mean_first_u(mean_first_u~=0);
                fprintf('mean first u: %f %f\n',mean(mean_first_u),u);
                fprintf('length of fisrt part when arrive at insufficient available %f\n',mean(L2));
                opt3 = find(status_flag==1);
                last_slot2 = last_slot;
                ww = find(first_slot~=0);
                last_slot2 = last_slot2(ww);
                qq = find(first_slot2~=0);
                last_slot = last_slot(qq);
                first_slot = first_slot(find(first_slot~=0));
                first_slot2 = first_slot2(find(first_slot2~=0));
                fprintf('mean of service time of small available time includes 1/p=1  :%f  Analytical:%f\n', mean(slot(last_slot2)-slot(first_slot)),ES_AvArrive)
                fprintf('mean of service time of small available time slots from end of end of first busy to last moment:%f  Analytical:%f\n', mean(slot(last_slot)-slot(first_slot2)),ES_AvArrive-u )
                fprintf('mean of service time of available time slots:%f  Analytical:%f\n', mean(sim1(opt3)),s1);   
                opt4 = find(status_flag==2);
                fprintf('mean of service time of busy time slots:%f  Analytical:%f\n', mean(sim1(opt4)),ES_BusyArrive+s1 + u/2); 
                opt5 = find(status_flag==3);
                fprintf('mean of service time of small available time slots:%f  Analytical:%f\n', mean(sim1(opt5)),ES_AvArrive+s1+mean(L2) );    
                fprintf('Si analysis: v:%f u:%f   s1:%f   ES:[%f  %f],  ES2::[%f  %f],var(s):[%f  %f] \n', v,u,s1, mean(sim1),ESk, ESk2, mean(sim1.^2), varsk2,var(sim1));
end
if 0
    clear all;
    k=1e6;
    u=6;v=2;
    p=0.1;
    N = geornd(p,[k,1])+1;
    X = zeros(k,100);
    Y = zeros(k,1);
    Z = zeros(1,k);

     if 1 %check exponential

         
            X  =   exprnd(u, [1,k]);   %arrive at busy
         fprintf('E:%f  mean:%f \n',u,mean(X));
%        histogram(X);
%         hold;
%         Z = exppdf(1:10000,u);
%         plot(Z);

     end
     if 0 %check theorem for exponential R.V


         for l=1:k
            X(l,1:1000)  =   exprnd(u, [1,1000]);   %arrive at busy
        end

        for l=1:k
            S=N(l,1);
            S=1/p;
            if(S>1000)
               fprintf('N exceeds max'); pause;
            end
            Y(l)  =   sum(X(l,1:S));   %arrive at busy
        end

         Z  =   gamrnd(1/p,u,[1,1000]);   %arrive at busy
        close all;
        figure;
        histogram(Y,'Normalization','pdf');
        hold all
        histogram(Z,'Normalization','pdf');
        legend('Simulation', 'Analytical')
        mean(Y)
        mean(Z)
        var(Y)
        var(Z)
     end

    if 0 %check theorem for gamma R.V using Welch?Satterthwaite equation
        N = geornd(p,[k,1])+1;
        if 1
              X  =   exprnd(u, [1,k]);   %arrive at busy
              
              X2  =   exprnd(v, [1,k]);   %arrive at busy
              
                    Yt = X+X2;   %arrive at busy
              
         else
            load('matlab.mat');
        end
      K =1;
      Egamma=(K)*((u)+(v));
      vargamma=(K^2)*(u^2+v^2);
      mean(Yt);
      var(Yt);
    
      Egamma=(K)*((u)+(v));
      vargamma=(K^2)*(u^2+v^2);
     fprintf('Egamma:%f  mean:%f vargamma:%f var:%f \n',Egamma,mean(Yt),vargamma,var(Yt));
     end
 end
if 0 %verify summation of   
    %a=mean channel avail time,   b=mean channel not avail time
    %s = s1 (with prob. a/a+b)   s=s1+x * ux (with prob b/a+b) and x=expo(b), ux=uniform(0,x)
    clear; clc;
    a=1;
    v=3;  %mean ch. available interval time
    u=6;  %mean ch. busy interval time
    n=243300;
    s1 = 4;
    p  = exp(-s1/v);
    nv=floor(n *p* v/(v+u));nsv=floor(n *(1-p)* v/(v+u)); nu=floor(n *u/(v+u));
    sv = s1 * ones(1, nv);          %avail interval
    ssv = zeros(1,nsv);
    su = zeros(1,nu);
    N = 1e6;
    r = geornd(p, [1, N])+1;

         vex_b = exprnd(v,max(r), N);
         uex_b = exprnd(u,max(r), N);
    for k=1:nu
        su(k) = sum(uex_b([1:r(k)], k))+ sum(vex_b([1:r(k)], k));   %arrive at busy
    end
  
         vex_a = exprnd(v,max(r)-1, N);
         uex_a = exprnd(u,max(r), N);

    for k=1:nsv
        ssv(k) =  sum(vex_a([1:r(k)-1], k))+sum(uex_a([1:r(k)], k));   %arrive at small available
    end


    fprintf('E[ssv]: %1.3f  = %1.3f   \n', mean(ssv),(((1/p))*(u))+(((1/p)-1)*(v)));
    a=1
    fprintf('E[ssv^2]: %1.3f  = %1.3f   \n', mean(ssv.^2),2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2);


    fprintf('E[su]: %1.3f  = %1.3f   \n', mean(su), ((1/p)-1)*((u)+(v)));
    fprintf('E[su^2]: %1.3f  = %1.3f   \n',mean(su.^2),2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2);
%  
     fprintf('var[su]: %1.3f  = %1.3f   \n',var(su),(2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2)-(v*(1/p-a)+u*(1/p-a))^2);
     fprintf('var[ssv]: %1.3f  = %1.3f   \n',var(ssv),(2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2+2*(u+v)*u+2*u^2)-(v*(1/p)+u*(1/p-a))^2);
%     %split with a uniform splitting
%     ux = rand([1, nu]);
%     wu = exprnd(u, [1, nu]).*ux; suu=s1 + wu + su;
%     vx = rand([1, nsv]);
%     wv = exprnd(u, [1, nsv]).*vx; svv=s1 + wv + ssv;
% 
%     sk = [sv, suu,svv];
%     ES_Ep= (1/p)-1;
% 
%     ES_BusyArrive = (1/p)*((u)+(v));
%     ES_AvArrive = ((1/p)*(u))+(((1/p)+1)*(v));
%     ES_BusyArrive2 = 2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2;
%     ES_AvArrive2 =(2*(u^2+v^2+u*v)*((1/p)-a)+(((2-p)/p^2)-((1+2*a)/p)+a^2+a)*(u+v)^2)-(v*(1/p-a)+u*(1/p-a))^2+2*v^2+ES_AvArrive^2;
%     var_BusyArrive = ((1/p)^2)*((u^2)+(v^2));
%     var_AvArrive = ((1/p^2)*(u^2))+(((1/p^2)+1)*(v^2));
% 
%     fprintf('S A :sigmau:%f %f sigmav:%f %f moment2u:%f %f moment2v:%f %fvaru:%f %f varv:%f%f\n',mean(su),ES_BusyArrive,mean(ssv),ES_AvArrive,mean(su.^2),ES_BusyArrive2,mean(ssv.^2),ES_AvArrive2,var(su),var_BusyArrive,var(ssv),var_AvArrive);
% 
%     Scenario1 = (p* v/(v+u)) .* s1
%     Scenario2 =(u./(v+u)).* (s1 + u/2+ ES_BusyArrive);
%     Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+v/2);
%     ESk = Scenario1 +Scenario2+Scenario3;
%     ESk2 =((p* v/(v+u)) .* s1^2)+((u./(v+u)).* (s1.^2 + 2*u.^2/3 +ES_BusyArrive2 + s1.*u + 2*s1*ES_BusyArrive +u*ES_BusyArrive))+((v./(v+u)).*(1-p)).* (s1.^2 + 2*v.^2/3 + s1.*v+ES_AvArrive2+2*s1*ES_AvArrive+v*ES_AvArrive);
%     varsk2 = ESk2- ESk.^2;   %ok
% 
%    if abs(var(sk)-varsk2)/abs(var(sk)+varsk2) > 0.01, warning('Potentialerror in calculating var(CRN:S) Err:%f<>%f \n', var(sk), varsk2); end
% 
%     fprintf('Si analysis: v:%f u:%f   s1:%f   ES:[%f  %f],  ES2::[%f  %f],var(s):[%f  %f] \n', v,u,s1, ESk, mean(sk), ESk2, mean(sk.^2), varsk2,var(sk));
end