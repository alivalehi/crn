function [tp, kp, fiv] = f_perform_framing(Framing_mode, x, TK, debug_active)  %excluded lambda and NinP
%This function receives input symbols, framing mode, and framing time/number 
% of symbols and generates packet process
% Inputs: 
%    Framing mode: 0:Frame based packetization, 1:Time based packetization 
%    x: inter-arrival times
%    TK: Framing time for time-based based packetization policy
%        or Number of symbols at each packet for number-based packetization 
%    lambda: input symbol rates

%a=1/lambda;
Nt = length(x);
t = zeros(1,Nt);   %arrivals of Poisson process
fiv = zeros(1,Nt);  %distance of symbols from end of a packet period


t(1) = x(1);
for i = 1:Nt-1,  
    t(i+1) = t(i) + x(i+1); %calc arrival times from interarrival times
end


%%
if Framing_mode == 1     %Time based             
    T=TK;
    Kp = ceil(t(end)/T);         %Number of packets for Time based framing
    t = t(t<=Kp*T); %truncate, extra check, is not necessary
    tp = zeros(1, Kp);   kp = zeros(1, Kp); 
elseif Framing_mode == 2     %Number based             
    K=TK; if abs(mod(K,1)) > 1e-6, error('K should be integer, while is : %d\n', K); end
    Kp = floor(Nt/K);   %Number of packets for Number-based framing: Number of symbols in a packet is fixed to NinP
    t=t(1:K*Kp);
    kp = K* ones(1, Kp); %each packet includes exactly k symbols, only time is mnot defined
else
    error('Invalid packetization policy: %d\n', TK);
end

%Packet process params
Nt=length(t);   fiv = zeros(1, Nt);

if Framing_mode == 1
    t0=t;  %initial t
    nPackets = 0; Tlen=0; Tf=[];
    for i = 1 : Kp
        tp(i) = i * T;     %end of current packetization interval
        ts = t(t <= i*T); %packets at this interval
        t = t(t > i*T);   %remove packetized symbols from input process
        if isempty(ts) 
            Tlen=Tlen+1;  %packet is not formed, the next packet interarrival is added by one time slot
        else   %form a new packet
            nPackets=nPackets+1;
            kp(i) = length(ts);   %number of symbols in the current interval
            fi = tp(i) - ts ;     %time of chosen symbols to end of interval
            fiv (sum(kp)-kp(i)+1:1:sum(kp)) = fi; 
            Tf(nPackets) = (Tlen+1)*T;  %packet length, geometrically distributed
            Tlen=0;   %reste
        end
    end
    Tf = mean(Tf);  %Was Tf = T;  %should include zero length packets
    %Theoretic value
    %Tf = T/(1-exp(-T/a));  %P(G=k)=p(1-p)^(k-1) => E[G]=1/p
    if debug_active,    fprintf('T/2:%1.6f      Framing delay:%1.6f  Deviation:%1.6f percent \n', Tf/2, mean (fiv), 200*( mean(fiv)-Tf/2)/Tf ); end

else   %Number based framing
        ts = reshape(t, K, Kp);
        tp = ts(K, :);
        fi = ones(K,1)*tp - ts; 
        fiv = reshape(fi ,1, K * Kp);
end
return

