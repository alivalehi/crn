function []= fplot_service_time_dist(lambda, T, N, H, BER, kmax, rmax)
%This program plots service time distribution for a packet with length L=kN+H bits
% ARQ system with PER=1-(1-BER)^L
if nargin < 2
    %default vars
    lambda = 5; 
    T = 2; N = 10;  H = 50;%250;
    BER = 1e-3; %3.0098e-04;
    kmax =  100; rmax = 10;
end
fprintf('Plot service time distribution using packet length and retransmission distribution.\n');

Pkr = zeros(kmax,rmax);
St = zeros(1, 1+kmax * rmax);
Ps = zeros(1, 1+kmax * rmax);


Rch = 1; %channel bit rate
S(1) = 0;  %service time
Ps(1) = exp (-lambda * T);

i = 1;

for k = 1 : kmax
    Pk = exp (-lambda * T) * ((lambda * T)^k / factorial(k)); %Prob of encapsulating k symbols
    Lp = N*k + H; %length of packet
    PER = 1 - (1-BER)^Lp;

    for rr = 1: rmax
        Pr = (1-PER) * (PER ^(rr-1)); %Prob of r retransmission        
        Pkr(k,rr) = Pk * Pr;

        i = i+1;
        St(i) = rr * (N*k + H) / Rch;
        Ps (i) = Pkr(k,rr); 
    end
end
[St, ind] = sort(St);
Ps = Ps(ind);


a = sum(Ps); b = Ps(1) + sum(sum(Pkr(:)));
fprintf('Check integration to 1 for p[s]:%d   p[r|k]:%d\n', a,b);

figure;
kv = 1:kmax; rrv = 1: rmax;
stem3(rrv, kv, Pkr);
%surf(rrv, kv, Pkr); 
xlabel('r'); ylabel('k'); zlabel('P[r,k]'); xlim([0, 10]); ylim([0, 40]); pause;

figure;
stem(St, Ps, 'r.'); xlabel('Service Time'); ylabel('P[s]'); xlim([0,600]); pause;
