
if 1 %this part randomly generate one  packet and verify service time of that packet
              clear variables;
              n=5;

                 P = zeros(1,n);
                
                P3 = zeros(1,n);
                P1 = zeros(1,n);
                P2 = zeros(1,n);
w = [3,.5,3.5,10,1];
chf= [0,1,0,1,0,1,0];
slot = [0,1,3,6,10,15,21];
s1=5;
                for i=1:n 
                sst = w(i);
                send_flag = 0;
                tag = 'Available'; 
                first_flag = 0;
                long_sav_flag = 0;
                
                while (~send_flag )
                        slot_s1 = max(find(slot <= sst)); slot_s2 = min(find(slot > sst));
                        slot_e1 = max(find(slot < sst+s1)); slot_e2 = min(find(slot >= sst+s1));
                      

                        if isempty(slot_e2)  %extend channel availabity process and retry
                            chf = [chf, chf(2:end)];
                            slot = [slot, max(slot)+ slot(2:end)];

                        elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (slot_s1 == slot_e1) && (chf(slot_e2) == 1)  %packet transmission fits into one availale time slot  
                            if(~first_flag)
                                status_flag(i)= 1 ;
                                first_flag = 1;
                                
                                P3(i)= P3(i)+1;
                       %         first_slot(i)= slot_e1;
                            end
                            sst = sst+ s1;  %set start time for next retransmission
                            send_flag = 1;   %mark the packet as transmitted
                            if(  status_flag(i)== 3)
                                last_slot(i)= slot_s1;
                            end
                            %      P(i)= P(i)+1;
                        elseif (slot_s1+1 == slot_s2) && (slot_e1+1 == slot_e2) && (chf(slot_s1) == 1)  %packet transmission arrive at busy time slot 
                             if(~first_flag)
                                status_flag(i)= 2 ;
                                first_flag = 1;                               
                                tag = 'Busy'; 
                                P2(i)= P2(i)+1;
                             else
                             %   P(i)= P(i)+1;
                             end   
                            next_av_ch = min(find((slot > sst) & (chf == 0)));  %chf==0, means the time slot right after unavailable channel
                          
                            if isempty(next_av_ch)  % extend the ch. avail process
                                chf = [chf, chf(2:end)];
                                slot = [slot, max(slot)+ slot(2:end)];
                                next_av_ch = min(find((slot > sst) & (chf == 0))); 
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
                                P1(i)= P1(i)+1;
                             else 
                                 sim2(i) = slot(slot_s2+1)-slot(slot_s2);
                                 if(  status_flag(i)== 3)
                                     if(long_sav_flag==0)
                                        first_slot(i)= slot_s1;  
                                        long_sav_flag = 1;
                                     end   
                                     P(i)= P(i)+1;
                                 end
                             end
                            next_av_ch = min(find((slot > sst) & (chf == 0)));  %chf==0, means the time slot right after unavailable channel
                            if isempty(next_av_ch)  % extend the ch. avail process
                                chf = [chf, chf(2:end)];
                                slot = [slot, max(slot)+ slot(2:end)];
                                next_av_ch = min(find((slot > sst) & (chf == 0))); 
                            end
                            sst = slot(next_av_ch); 
                            if isempty(next_av_ch), queue.fail=1; end
                        else
                            uncondition=1;
                        end
                    end  %end while          while (~send_flag && ret_cont<100)
                    
                    sim1(i) = sst-w(i);
                
                end

               
                  end
if 0
    clear all;
    k=10000;
    u=4;v=1;
    p=0.1;
    N = geornd(p,[k,1])+1;
    X = zeros(k,100);
    Y = zeros(k,1);
    Z = zeros(1,k);

     if 0 %check exponential

          for l=1:k
            X  =   exprnd(u, [1,10000]);   %arrive at busy
          end
       histogram(X);
        hold;
        Z = exppdf(1:10000,u);
        plot(Z);

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

    if 1 %check theorem for gamma R.V using Welch?Satterthwaite equation
        N = geornd(p,[k,1])+1;
        if 1
              X(1:k,1:k)  =   exprnd(u, [k,k]);   %arrive at busy
              for l=1:k
                  S=N(l,1);
                  if(S>k)
                    S = k;
                  end
                  Y(l,L) = sum(X(l,1:S));   %arrive at busy
              end
              X2(1:k,1:k)  =   exprnd(v, [k,k]);   %arrive at busy
              for l=1:k
                  S=N(l,1);
                  if(S>k)
                    S = k;
                  end
                  Y2(l,L) = sum(X2(l,1:S));   %arrive at busy
              end
          
        Yt= Y+Y;
        else
            load('matlab.mat');
        end
      K =1/p;
      Egamma=(K)*((u)+(v))
      vargamma=(K^2)*(u^2+v^2)
      mean(Yt)
      var(Yt)
     end
 end
if 0 %verify summation of   
    %a=mean channel avail time,   b=mean channel not avail time
    %s = s1 (with prob. a/a+b)   s=s1+x * ux (with prob b/a+b) and x=expo(b), ux=uniform(0,x)
    clear; clc;
    a=1;
    v=2;  %mean ch. available interval time
    u=5;  %mean ch. busy interval time
    n=1e6;
    s1 = 3;
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


    fprintf('E[su]: %1.3f  = %1.3f   \n', mean(su), (1/p)*((u)+(v)));
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