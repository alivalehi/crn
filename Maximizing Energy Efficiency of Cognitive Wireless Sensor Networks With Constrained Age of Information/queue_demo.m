clc
clear all
close all
disp('This program simulates a single queue M/M/1 with deadlines [only for demo: not required for resuls]');

STE_ACTIVE = 0;
DEADLINE_ACTIVE = 0;


%%
use_default = 0;
if use_default

    tend = 25;
    fi = 1; %deterministic
    fo = 1; %deterministic
    fg = 1; %0: no deadline 1:deterministic, 2:poisson
    

%generate input process
    N = 6;
    Tint = [0, 1, 5, 1, 1, 1, 1];
    gi = [4, 4, 2, 1, 1, 100];
    si = [5,3, 5 , 4, 1 , 1];
else

    tend = input('please enter end time (0..100, def:5): ');
    if isempty(tend)
        tend = 5;
    end

    fi = input('Please enter input distribution (1:deterministic, 2:poisson), [def:2]: ');
    if isempty(fi)
        fi = 2;
    end

    fo = input('Please enter service time distribution (1:deterministic, 2:exponential), [def:2]: ');
    if isempty(fo)
        fo = 2;
    end

    fg = input('Please enter deadline distribution (0: no deadline 1:deterministic, 2:poisson), [def:0]: ');
    if isempty(fg)
        fg = 0;
    end

    x = input('Enter input rate , [def:5]: ');
    if isempty(x) 
        x = 2;
    end

    mu = input('Enter service rate, , [def:8]: ');
    if isempty(mu) 
        mu = 8;
    end

    if DEADLINE_ACTIVE && (fg > 0 )
        mu_g = input('Enter deadlines: (def = 2*x): ');
        if isempty(mu_g) 
            mu_g = 2*x;
        end
    end


    %generate input process
    %dt= (1/10)*min(1/x,1/mu);
    N = floor(tend*x); %number of arrivals
    if fi == 2
        Tint = exprnd((1/x), 1, N); 
    else
        Tint = 1/x * ones(1,N);
    end
    
    if DEADLINE_ACTIVE
        %generate deadlines 
        if fg == 2
            gi = exprnd(1/mu_g, [1,N]); 
        elseif fg == 1
            gi = (1/mu_g)* ones(1,N);
        end
        zi = zeros(1,N); %drop rate
        
    end
    
    %generate output process, service times
    if fo == 2
        si = exprnd(1/mu, [1,N]); 
    else
        si = (1/mu)* ones(1,N);
    end
end
dt= (1/x)/50;
ti = zeros(1,N);  %arrival times
wi = zeros(1,N);  %waiting times
zi = zeros(1,N);  %discard flag
to = zeros(1,N);  %departure time





%%
%generate input arrival times
ti(1) = 0;
to(1) = si(1)
pause;
for i=1:N
    ti(i+1) = ti(i) + Tint(i);  %next arrivaltime
end


wi(1) = 0;
zi(1) = 0;  %Packet 1 is not discarded
for i=1:N-1

    %calculatewaiting time
    if DEADLINE_ACTIVE
        if STE_ACTIVE
            %put the packet in appropriate position
            %save the new order 
            [ei, sort_index] = sort(ei);
            error('Later check carefully');
        end

        if zi(i) == 1 %previous packet is discarded 
            wi(i+1) = wi(i) - Tint(i); %the time required to have the queue empty, previous packet was discarded
        else
            wi(i+1) = wi(i) + si(i) - Tint(i); %the time required to have the queue empty
        end

    else  %no deadline
        wi(i+1) = wi(i) + si(i) - Tint(i); %the time required to have the queue empty
    end
    
    if wi(i+1) < 0  %if arrives after queue gets empty => no waiting time
        wi(i+1) = 0;  
    end
        
     if (DEADLINE_ACTIVE) && (fg > 0) && (wi(i+1) >= gi(i+1)) %(wi(i)+si(i)>= gi(i))  %discard the packet
        zi(i+1)= 1;
    else
        zi(i+1)= 0;
    end
    
    %generate departure process
    to(i+1) = ti(i+1) + wi(i+1) + si(i+1);
end
%remove the last extra packet
ti = ti(1:end-1);

si_Ti = si - Tint; %this amount adds to waiting time if positive

%number of arrivals
ai = 1:1:N;
bi = 1:1:N;

%%
if DEADLINE_ACTIVE
    ei = ti + gi;  %extinction time
    if (fg > 0)
        Nf = sum(zi); %Fail
    else
        Nf = 0; %Fail due to deadline
    end
    Ns = N-Nf; %Sucess



    is = 0;
    for i=1:N
        if zi(i) == 1  %discarded packet
        else
            is= is+1;
            tis(is)= ti(i);
            wis(is) = wi(i);
            tos(is)= to(i);
        end
    end

    ais = 1:1:Ns;
    bis = 1:1:Ns;
    [tios nis]= combine_arriv_depart(tis, tos);
end

if use_default
    fprintf('[Tint;ti;gi;ti+gi;wi;si;to;zi]');
    matrix =[Tint;ti;gi;ti+gi;wi;si;to;zi]
    pause;
end

[tio ni]= combine_arriv_depart(ti, to);


fprintf('\n');

figure;
subplot(211);
splot(ti,ai, 'b');
hold on;
pause;
%better plot from origin
splot(to,bi, 'r');
title('Arrivals, blue A(t), Departures, red B(t)');
pause;
splot(tio,ni+0.1, 'g');

subplot(212);
title('N(t)');
splot(tio,ni, 'g');
pause;



figure;
subplot(411);
splot(ti,si, 'b');
ylabel('Service time, S(t)');

subplot(412);
splot(ti,si_Ti, 'b');
ylabel('Service time - Inter arrival time, S(t)');

% subplot(413);
% ylabel('Service time - Inter arrival time, S(t)');
% splot(ti,si_Ti, 'b');

subplot(414);
splot(ti,wi, 'g');
ylabel('Waiting Time, W(t)');
pause;

if DEADLINE_ACTIVE
    figure;
    subplot(211);
    splot(tis,ais, 'b');
    hold on;
    pause;
    %better plot from origin
    tos = [0, tos];
    bis = [0, bis]; 
    splot(tos,bis, 'r');
    pause;
    splot(tios,nis+0.1, 'g');
    subplot(212);
    splot(tios,nis+0.1, 'g');
    Discard_Rate= Ns/N
end





Nav = calcav(tio, ni,dt);
Wav = sum(wi)/length(wi);

if (fi == 2) && (fo ==2)%M/M/1
    ru = min(1,mu/x);
    Navth= x/(mu-x);
    fprintf('Average N = (%1.3f , %1.3f) \n', Nav, Navth);
    
    Wavth = ru/(mu-x);
    fprintf('Average W = (%1.3f , %1.3f) \n', Wav, Wavth);

end
