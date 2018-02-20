function [optimum_k] = callNB(u,v,Ber,i,k)

%This program runs simulations for different parameter setups
%see also call_framing_delay for similar test setups

%clear all; close all; clc;
quick_run=1;
sim.method = 1;
%#######################################################################
%     SIMULATION OPTIONS
%#######################################################################
if (sim.method == 1) %method1
    
    sim = init_sim(2);
    sim.variable =1;
    sim.u = u;
    sim.v  = v;
    
    sim.ch.BER = Ber;
    sim.H = 80;%228;
    sim.N = 16;%50;
    sim.SU =1;
    if(sim.variable==1)
%         sim.Nexperts = 1e4;
%         sim.ch.Rch =1000;
%         sim.nch = 1e4;
%         sim.k_min=1 ;
%         sim.k_max= 100;
%         %   sim.u = [.25:.05:4];%[1 1 1 1 1 1 1 1 2 3 4 5 6 7 8];%randi([4,10],1,sim.m);
%         %   sim.v = ones(1,length(sim.u))*2;%[8 7 6 5 4 3 2 1 1 1 1 1 1 1 1];%randi([4,10],1,sim.m);
%         sim.m=length(sim.u);
%         sim.run.Nsym = 1e5;
%         sim.lambda = 10;%randi([1,10],1,sim.m);
%         Berr = [0 .1 .01 .001 .0001];
%         %  sim.ch.BER = ones(1,sim.m).*1e-4;%Berr(randperm(length(Berr),1));%Berr(randperm(length(Berr),sim.m));
        n=1;
        
        sim.v = [1:.5:3 2.5:-.5:1];%repelem([1,2,3],2*n);%1*ones(1,2*n);%repmat([.5,1.5],[1,n]);%randi([1,15],1,n);%[2 5 1 4];
        sim.u = [.5 .5 .5 .5 .5 .5 .5 .5 .5];%repelem([1,2,3],2*n);%1*ones(1,length(sim.v));%randi([1,15],1,n);%[2 4 6 8];
        
        sim.Nexperts = 1e3;
        sim.m=length(sim.u);
        sim.nch = randi([5 10],sim.m,1);%sim.nch = 1e4;
        sim.k_min=1 ;
        sim.k_max= 100;
        %   sim.u = [.25:.05:4];%[1 1 1 1 1 1 1 1 2 3 4 5 6 7 8];%randi([4,10],1,sim.m);
        %   sim.v = ones(1,length(sim.u))*2;%[8 7 6 5 4 3 2 1 1 1 1 1 1 1 1];%randi([4,10],1,sim.m);
        
        sim.run.Nsym = 1e8;
        sim.lambda = 1;%randi([1,10],1,sim.m);
        
        %  sim.ch.BER = ones(1,sim.m).*1e-4;%Berr(randperm(length(Berr),1));%Berr(randperm(length(Berr),sim.m));
        sim.lambda = .1;
        sim.ch.Rch =1e2;
        sim.Nparamset = 1;
        
        sim.Nparamset = 1;
    else
        sim.lambda = 5;
        sim.Nexperts = 1e2;
        sim.ch.Rch =1e4;
        sim.nch = 5e6;
        sim.k_min=1 ;
        sim.k_max= 30;
        sim.run.Nsym = 1e10;
        sim.Nparamset = Inf;
        sim.m=1;
        %sim.u = 4;
        %sim.v = 8;
        %  sim.ch.BER = 7e-4;
    end
    
    sim.power_bit = 10;
    sim.power_sense = 20;
    kappa = 1;
    sim.phi1 = 2.05;
    sim.phi2 = 2.05;
    sim.phi = sim.phi1 + sim.phi2;
    sim.chi = 2*kappa/abs(2-sim.phi-sqrt(sim.phi^2-4*sim.phi));
    %     sim.alpha = .5+(rand/2);%chi;% coeffient of velocity
    %     sim.beta = chi*rand*phi1; %coeffient of Memory Best
    %     sim.gamma = chi*rand*phi2; %coeffient of current best    sim.BER_Mbest = 0;
    sim.u_Mbest = [];
    sim.v_Mbest = [];
    sim.v_past = [];
    sim.u_past = [];
    sim.BER_past = [];
    sim.omega = .65;
    sim.omega1=.33;
    sim.omega2=.33;
    sim.omega3=.33;
    sim.channel_change_flag =0;
    sim.run.NUM_RUNS=1;
    %########
    sim.n =3e2; %number of packets for each segment simulation 
    sim.windowsshift = 30;
    %########
    sim.Itlimit =2000;
    sim.i=i;
    sim.umin =.1;
    sim.umax =10;
    sim.vmin =.1;
    sim.vmax =10;
    sim.kcoeff = log([1:((exp(1)-1)/sim.n):exp(1)]);
    %#######################################################
    % calculate simulation papram
    %#######################################################
    nump = sim.run.Nsym /sim.k_max;
    T = 1/sim.lambda;
    Plen = sim.k_max* sim.N + sim.H;
    ret =(((1 - sim.ch.BER) .^(-Plen)))%retransmisson
    S = Plen/ sim.ch.Rch
    var1 = (S+(T*sim.k_max))*sim.n*sim.Itlimit*ret;
    var2 = sim.nch*(sim.u+sim.v);
    %fprintf('var1:%f var2:%f \n',var1,var2);
    %     sim.alpha = .5+(rand/2);%chi;% coeffient of velocity
    % sim.beta = sim.chi*rand*sim.phi1; %coeffient of Memory Best
    % sim.gamma = sim.chi*rand*sim.phi2; %coeffient of current best
    %     sim.alpha = .333;%chi;% coeffient of velocity
    % sim.beta = .333;%coeffient of Memory Best
    % sim.gamma = .333; %coeffient of current best
    if(k==0)
        sim.k = randi([sim.k_min,sim.k_max],1);
    else
        sim.k = k;
    end
    optimum_k = fqueu_NB1_constant(sim); %fqueu_NB1_v4mix(sim);
    t=clock;
fprintf('\n\n***********************************************************************************\n');
fprintf('End: Date:%d-%d-%d   Time: %d:%d \n', floor(t(1:5)));
fprintf('***********************************************************************************\n');
end
%######################################################################
%                   Method 2
%######################################################################
sim.method = 1;
if (sim.method == 2)
    
    sim = init_sim(2);
    sim.variable =0;
    sim.u = u;
    sim.v  = v;
    n= 100 ;
    sim.ch.BER = Ber
    sim.H = 228;
    sim.N = 50;
    sim.SU =1;
    if(sim.variable==1)
        sim.u = randi([1,15],1,n);%[2 4 6 8];
        sim.v = randi([1,15],1,n);%[2 5 1 4];
        sim.Nexperts = 1e3;
        sim.nch = randi(10,1e4,1);%sim.nch = 1e4;
        sim.k_min=1 ;
        sim.k_max= 20;
        %   sim.u = [.25:.05:4];%[1 1 1 1 1 1 1 1 2 3 4 5 6 7 8];%randi([4,10],1,sim.m);
        %   sim.v = ones(1,length(sim.u))*2;%[8 7 6 5 4 3 2 1 1 1 1 1 1 1 1];%randi([4,10],1,sim.m);
        sim.m=length(sim.u);
        sim.run.Nsym = 1e5;
        sim.lambda = .01;%randi([1,10],1,sim.m);
        
        Berr = [0 .1 .01 .001 .0001];
        %  sim.ch.BER = ones(1,sim.m).*1e-4;%Berr(randperm(length(Berr),1));%Berr(randperm(length(Berr),sim.m));
        sim.lambda = 1;
        sim.ch.Rch =1e3;
        sim.Nparamset = 1;
    else
        sim.lambda = .01;
        sim.Nexperts = 1e6;
        sim.ch.Rch =1e3;
        sim.nch = 5e6;
        sim.k_min=1 ;
        sim.k_max= 20;
        sim.run.Nsym = 1e6;
        sim.Nparamset = Inf;
        sim.m=1;
        %sim.u = 4;
        %sim.v = 8;
        %  sim.ch.BER = 7e-4;
    end
    
    sim.power_bit = 100;
    sim.power_sense = 200;
    kappa = 1;
    sim.phi1 = 2.05;
    sim.phi2 = 2.05;
    sim.phi = sim.phi1 + sim.phi2;
    sim.chi = 2*kappa/abs(2-sim.phi-sqrt(sim.phi^2-4*sim.phi));
    %     sim.alpha = .5+(rand/2);%chi;% coeffient of velocity
    %     sim.beta = chi*rand*phi1; %coeffient of Memory Best
    %     sim.gamma = chi*rand*phi2; %coeffient of current best    sim.BER_Mbest = 0;
    sim.u_Mbest = [];
    sim.v_Mbest = [];
    sim.v_past = [];
    sim.u_past = [];
    sim.BER_past = [];
    sim.omega = .65;
    sim.omega1=.33;
    sim.omega2=.33;
    sim.omega3=.33;
    sim.channel_change_flag =0;
    sim.run.NUM_RUNS=1;
    sim.n = 1e3; %number of packets to simulation result
    sim.Itlimit =100;
    sim.i=i;
    sim.umin =1;
    sim.umax =10;
    sim.vmin =1;
    sim.vmax =10;
    sim.kcoeff = log([1:((exp(1)-1)/sim.n):exp(1)]);
    %#######################################################
    % calculate simulation papram
    %#######################################################
    nump = sim.run.Nsym /sim.k_max;
    T = 1/sim.lambda;
    Plen = sim.k_max* sim.N + sim.H;
    ret =(((1 - sim.ch.BER) .^(-Plen)))%retransmisson
    S = Plen/ sim.ch.Rch
    var1 = (S+(T*sim.k_max))*sim.n*sim.Itlimit*ret;
    var2 = sim.nch*(sim.u+sim.v);
    % fprintf('var1:%f var2:%f \n',var1,var2);
    %     sim.alpha = .5+(rand/2);%chi;% coeffient of velocity
    % sim.beta = sim.chi*rand*sim.phi1; %coeffient of Memory Best
    % sim.gamma = sim.chi*rand*sim.phi2; %coeffient of current best
    %     sim.alpha = .333;%chi;% coeffient of velocity
    % sim.beta = .333;%coeffient of Memory Best
    % sim.gamma = .333; %coeffient of current best
    if(k==0)
        sim.k = randi([sim.k_min,sim.k_max],1);
    else
        sim.k = k;
    end
    fprintf('system number:%d',i);
    optimum_k = fqueu_NB(sim);
    t=clock;
fprintf('\n\n***********************************************************************************\n');
fprintf(' Date:%d-%d-%d   Time: %d:%d \n', floor(t(1:5)));
fprintf('***********************************************************************************\n');
end
ALI=2;
end



