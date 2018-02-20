clear; %clc; 
datapath = [cd, '\PCodeN\NData\CRN_K_Delay\'];


   



% if 1   %uncoded system
%     fname1 = 'NB1_NH8_64_SymRate30_Rch500_BER10000_CS11_CR1_1_impBER1_10_CRN0_0_1_Kv360.mat'; 
%     user_command='plot_delays=1'; user_condition = '';
%     plot_data_fileNB(datapath,fname1, 0, user_command, user_condition);   %active flags:  plot_delays
% end
%%ali
if 0   %uncoded system
    fname1 = 'NB1_NH8_64_SymRate30_Rch1o253804e+03_BER10000_CS11_CR1_1_impBER1_10_CRN1_1_40_Kv2535.mat'; 
    user_command='plot_delays=1'; user_condition = '';
    plot_data_fileNB(datapath,fname1, 0, user_command, user_condition);   %active flags:  plot_delays
end



if 0  %coded system
    fname2 ='NB1_NH8_64_SymRate30_Rch1000_BER100000_CS11_CR2_2_impBER10_10_CRN0_0_1_Kv360.mat'; 
    user_command='plot_delays=1'; user_condition = '';
    %fname2 ='NB1_NH8_64_SymRate30_Rch1500_BER100000_CS11_CR2_3_impBER10_100_CRN0_0_1_Kv360.mat'; 
    plot_data_fileNB(datapath,fname2, 0, user_command, user_condition);   %active flags:  plot_delays
end


 

if 0   %comp_ber
    
datapath = [cd, '\NData (another copy)\'];
    f=dir('NData\NB1_NH8_64_SymRate30_Rch1000_*.mat');  
    fname_ber={};    for i=1:length(f), fname_ber{i}=f(i).name; end;
    
    fname_ber = {'NB1_NH8_128_SymRate2o000000e-01_Rch100_BER100_CS11_CR1_1_impBER1_10_CRN1_2o000000e-01_1_Kv465.mat' ...
             'NB1_NH8_128_SymRate2o000000e-01_Rch100_BER112_CS11_CR1_1_impBER1_10_CRN1_2o000000e-01_1_Kv465.mat', ...
             'NB1_NH8_128_SymRate2o000000e-01_Rch100_BER125_CS11_CR1_1_impBER1_10_CRN1_2o000000e-01_1_Kv465.mat', ...
             'NB1_NH8_128_SymRate2o000000e-01_Rch100_BER143_CS11_CR1_1_impBER1_10_CRN1_2o000000e-01_1_Kv465.mat', ...
             'NB1_NH8_128_SymRate2o000000e-01_Rch100_BER167_CS11_CR1_1_impBER1_10_CRN1_2o000000e-01_1_Kv465.mat', ...
             'NB1_NH8_128_SymRate2o000000e-01_Rch100_BER200_CS11_CR1_1_impBER1_10_CRN1_2o000000e-01_1_Kv465.mat'};
         
    user_command='plot_simulation=0; plot_delays=1'; user_condition = '1';%'ismember(sim.ch.BER, [0 1e-3 2e-3 3e-3 4e-3 5e-3])';
    plot_data_fileNB(datapath,fname_ber, 0, user_command, user_condition);   %active flags:  plot_delays=1   plot_simulation = 0;
    %for i=1:length(fname_ber), plot_data_fileNB(datapath,fname_ber{i}, 0); end;  %active flags:  plot_delays
end


if 0 %ber_k
    load([cd,'\NData\BER_K\data_ber_k.mat']);
    %stem(BER_vect,K_vect)  
    L = min(find(isnan(K_vect))); BER_vect=BER_vect(1:L); K_vect= K_vect(1:L);
    [K,ind]= unique(K_vect); B=BER_vect(ind);
    disp('limits for BER to change K: '); [B; K], pause;
    %Bmid = 0.5 * ([B, 0] + [0, B]); Bmid=Bmid(1:end-1);

      
    %bar(Bmid,K)
    figure; plotdef(0)
    for i=1:length(K)-1
        l=B(i); u=B(i+1); k=K(i);
        area([l, u], [k, k]); hold on;
    end        
    xlabel('Channel Bit Error Probability \beta'); ylabel('Optimal Framing Parameter: k');
    
end


if 0  %CRN, avail time fixed
%     f=dir('NData\NB0_NH8_64_SymRate30_Rch1000000_BER100000_CS11_CR2_2_impBER10_10_CRN1_*.mat');  
%     flist={};    for i=1:length(f), flist{i}=f(i).name; end;
    flist = { 'NB0_NH8_128_SymRate2o000000e-01_Rch1000000_BER100000_CS11_CR2_2_impBER10_10_CRN1_100_0_Kv1275.mat',...
              'NB0_NH8_128_SymRate2o000000e-01_Rch1000000_BER100000_CS11_CR2_2_impBER10_10_CRN1_10_0_Kv1275.mat', ...
              'NB0_NH8_128_SymRate2o000000e-01_Rch1000000_BER100000_CS11_CR2_2_impBER10_10_CRN1_5_0_Kv1275.mat', ...
              'NB0_NH8_128_SymRate2o000000e-01_Rch1000000_BER100000_CS11_CR2_2_impBER10_10_CRN1_1_0_Kv1275.mat', ...
              'NB0_NH8_128_SymRate2o000000e-01_Rch1000000_BER100000_CS11_CR2_2_impBER10_10_CRN1_1o000000e-01_0_Kv1275.mat'};
         
    user_command='plot_simulation=0; plot_delays=1'; user_condition = 'ismember(sim.CRN.chNotAvRatio, [5 10 100])'; %[1 5 10 100])'  %user_condition = 'ismember(sim.CRN.chNotAvRatio, [0.1 1 5 10 100])';
    plot_data_fileNB(datapath,flist, 0, user_command, user_condition);   %active flags:  plot_delays=1   plot_simulation = 0;
    
    %user_command='plot_simulation=0; plot_delays=1'; user_condition = 'ismember(sim.CRN.chNotAvRatio, [0.1 1 5 10])';
    %xlim([0 20]);
    %for i=1:length(flist), plot_data_fileNB(datapath,flist{i}, 0, user_command, user_condition); end;  %active flags:  plot_delays
    
end


if 1
    load NData\CRN_ES\CRD_data
    close all; plotdef(0); %figure;
    for k = [1,5:5:20]
        Esk = ES(:, k)'; loglog(chNotAvRatiov, Esk); hold on; pause;
    end
    legend('k:1', 'k:5', 'k:10','k:15', 'k:20');
    xlabel('Channel Unavailability Factor: \rho_{ch}=u/v'); ylabel('Service Time: E[S]');
end
