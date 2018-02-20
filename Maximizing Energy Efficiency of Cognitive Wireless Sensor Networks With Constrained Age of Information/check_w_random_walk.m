close all;
clc;
clear all;
disp('This program shows how wn+1=max(0, wn-1) can be modeled as');
disp('sup of accumulated steps of a random walk');
disp('wn = max(v0,v1,v2,...)');
disp('extra demo, not required for main results');
lambda = 1;
a = 1\lambda;
Nsym = 1000;
T = a * 1.1;

plot_active = 1;

Ew1 = zeros(1,100);
Ew2 = zeros(1,100);
Ew3 = zeros(1,100);
for run = 1: 1000
    if (mod(run,100)==0), fprintf('\nx(%d/%d)',floor(run/100) ,10); end

    x = random('exp', a,[1,Nsym]);  %interarrivals
    u = x - T;

    w1 = zeros(1,Nsym);  %waiting time, calc directly
    w2 = zeros(1,Nsym); %waiting time, calc indirectly using backward random alk process
    w3 = zeros(1,Nsym); %sum of accumulated forward random walk process

    %Expectation: w1=w2  E(w1)= E(w3)

    wsum = zeros(1,Nsym);


    w1(1) = 0;
    w2(1) = 0;
    w3(1) = 0;
    wsum(1) = 0;
    for i = 1: Nsym
        if (mod(i,1000)==0), fprintf('.'); end
        w1(i+1) = max(w1(i) + u(i), 0);


        v = zeros(1,i);

        v(1) = u(i);
        for j = 2 : i
            v(j) = v(j-1) + u(i - j + 1);
        end

        w2(i+1) = max(0, max(v));


        wsum(i+1) = wsum(i) + u(i);

        w3(i+1) = max(wsum);

    %     wcomp=[w;w2];
    %     pause;
    end

    Tv = T * (1:1:Nsym);
    Zv = zeros(1,Nsym);

    if (run == 1) && (plot_active)
        figure;
        subplot(311);
        stem(Tv, u, '.');
        %splot(Tv, u, 'b');
        hold on;
        plot(Tv,Zv,'k:');
        title('U = S - T');


        subplot(312);
        splot(Tv, wsum(2:end), 'b');
        hold on;
        plot(Tv,Zv,'k:');
        title('Waiting Time, no max');

        subplot(313);
        splot(Tv, w1(2:end), 'b');
        hold on;
        plot(Tv,Zv,'k:');
        %splot(Tv, w2(2:end)-0.2, 'r');
        title('Waiting Time');

        pause;

        figure;
        subplot(211);
        splot(Tv, w1(2:end), 'b');
        hold on;
        plot(Tv,Zv,'k:');
        title('Waiting Time');

        subplot(212);
        lastind = min(5+ max(find(w3 ~= w3(end))), length(Tv));
        splot(Tv(1:lastind), w3(1+(1:lastind)), 'b');
        hold on;
        plot(Tv(1:lastind),Zv(1:lastind),'k:');
        title('Waiting Time: Obtained From Random walk proc');
    end

    Ew1(run) = mean(w1);
    Ew2(run) = mean(w2);
    Ew3(run) = mean(w3);
end
E = [mean(Ew1), mean(Ew2), mean(Ew3)]
pause;

