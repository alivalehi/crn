if(1)%new
    clear all
    already_saved = 0;
    Kpredict_array = [];
    BERpredict_array = [];
    upredict_array = [];
    vpredict_array = [];
    EE_real_track_array = [];
    EE_diff_best_track_array = [];
    kAnal_array = [];
    k_array = [];
    EE_prediction_vector1_array = [];
    n= 4; %10 for method 1
    sim.u = ([13 14 2 14]);%randi([1,15],1,n);
    %[2 4 6 8];%[3:.5:7];% 4 5 6 7 2 6 12 1 2 4];
    %[.25:.1:4 fliplr(.25:.1:4)];%[1 1 1 1 1 1 1 1 2 3 4 5 6 7
    %8];%randi([4,10],1,sim.m);[2 6 12 1 2 4]
    sim.v = ([3 5 15 8]);%randi([1,15],1,n);
    %[2 4 6 8];%[4:.5:8];
    rho = (sim.v) ./(sim.u)
    %[fliplr(sim.u) .25:.1:4];%ones(1,length(sim.u))*3;%[8 7 6 5 4 3 2 1 1
    %1 1 1 1 1 1];%randi([4,10],1,sim.m);[1 2 3 4 6 8]
    sim.BER = 1e-4;
    nn = 1; length(sim.u)
    uu = [];
    vv = [];
    k=0;
    if(already_saved ==0)
        for i=1:nn
            u= sim.u;
            v = sim.v;
            Ber = sim.BER;
            [optimum_k] = callNB(u(i),v(i),Ber,i,k);
            k =optimum_k; 
        end
    end
    datestr(now) 

end 