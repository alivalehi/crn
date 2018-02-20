function [sim] = fanalyzesimTB(sim)
%This program gets the simulation paramteres and then calculates
%(analythicaly) the expected delay [for Time based system]
% step 1: calculateing E(s) and E(s^2)
% step 2: using kingman formula 


adjust_kingman = 0;
newEf = 1;

if nargin < 1
    %default vars
    lambda = 10; N = 3;    H = 10;
    
    T = 5/lambda;
    Exp_PER = 0.2;
    Tv = (1/lambda) * [0.1:0.1:1,    2:1:10];
else
    H = sim.H;
    N = sim.N;
    lambda = sim.lambda;

    Tv = sim.Tv;
    if isempty(sim.Tv2)
        Tstep = 0.01/lambda; Tv2 = Tv(1): Tstep: Tv(end); sim.Tv2=Tv2; %for analytical results, smaller time resolutions can be used, no need for heavy run as needed for MC simulations
        Tv2 = Tv; sim.Tv2=Tv2; %same resolution for analytical
    else
        Tv2 = sim.Tv2;
    end

    
    if (sim.Framing_mode ~= 1) || (sim.input_process ~= 1)
        error('Err: Invalid params');
    end
    
    BER = sim.ch.BER;    
    Rch = sim.ch.Rch;
    
end

% Plen_av = (T* lambda * N + H);
% Rb_av = Plen_av/T;
% Rch = (1/(1-Exp_PER))* Rb_av;
% %info_bit_rate = lambda * N;
% BER = 1 - (1 - Exp_PER) ^ (1/Plen_av);
fprintf('This program calculates queuing parameters analytically.\n');

LT = lambda * Tv;
BERv = [BER 0];

H0=H; N0=N;
for i = 1:2  %calculate for the given BER and also error-free channel
    B = BERv(i);
    a0 = 1 - B;
    if sim.CodedSystem   %Modifications
        
        consider_data_tail_bits_in_header=1;
        if consider_data_tail_bits_in_header
            sim.coding.nTbitsH = sim.coding.nTbitsH+sim.coding.nTbitsD; sim.coding.nTbitsD=0; %Tails bits in header
        else
            sim.coding.nTbitsD = sim.coding.nTbitsD .* (LT >= 1);
        end
        
        fprintf('uncoded H:%d  N:%d  alpha:%1.4f   Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f]  ', H0,N0,a0, sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD);
        %find equivalent alpha, beta, H and N in order to use the same equations
        
        H = (H0 / sim.coding.CRH + sim.coding.nTbitsH); if length(H)==1, H=repmat(H, [1,length(LT)]); end  %equivalent H to use equations for uncoded system
        N = (N0 / sim.coding.CRD + sim.coding.nTbitsD./ LT); if length(N)==1, N=repmat(N, [1,length(LT)]); end  
        a = ((1 - B / sim.coding.impBERH).^(H)) .*((1-B/sim.coding.impBERD).^(N.*LT));
        a = a .^ (1./(H0 + N0 .* LT));
        fprintf('\n   Coded H:%s  N:%s  alpha:%s\n', num2str(H([1:3, end-1, end])),num2str(N([1:3, end-1, end])),num2str(a([1:3, end-1, end]))); % pause;
    else
        a = a0;
        fprintf('uncoded H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f', H,N,a, B, sim.ch.Rch);
    end
    Es = (N/(Rch * a.^H)) .* ((H./N + LT./(a.^N)).* exp(-LT.*(1-a.^(-N))) - (H./N).* exp(-LT) ) ;
    
    
            if 0   %%only for temp debug [later remove]
                if Tv > 100/lambda   %Approximation for large T
                    Es=Es

                    H = H*3;
                    Es_3H = (N/(Rch * a^H)) * ((H/N + LT/a^N).*exp(-LT*(1-1/a^N)) - (H /N)* exp(-LT) ) 
                    H=H/3;

                    Es_test = (N * LT/(Rch * a^(H+N))).*exp(-LT*(1-1/a^N)) 
                    Es_test = (N * LT/(Rch * a^(H+N)))

                    a=a
                    aN=a^N

                    t1=LT
                    t2 = (1-1/a^N)
                    pause;
                    term1=exp(-LT*(1-1/a^N))
                    term2=exp(-LT)
                    pause;
                end
            end


    %% paper notation: Accurate Result
        AN2_1 = (2*(N.^2)./a.^(2*H)).*(LT.*a.^(-2*N) + LT.^2 .* a.^(-4*N)) .* exp(-LT.*(1-a.^(-2*N))) ;
        AN2_2 = -((N.^2)./a.^(H)).*(LT .* a.^(-N) + LT.^2 .* a.^(-2*N)) .* exp(-LT.*(1-a.^(-1*N))) ;
        AN2 = (AN2_1+AN2_2) / Rch^2;

        AKHN_1 =  (4*H.*N./a.^(2*H)).*(LT.*a.^(-2*N)) .* exp(-LT.*(1-a.^(-2*N))) ;
        AKHN_2 = -(2*H.*N./a.^(H)).*(LT.*a.^(-N)) .* exp(-LT.*(1-a.^(-1*N))) ;
        AKHN = (AKHN_1 + AKHN_2) / Rch^2;

        AH2_1 = (2*H.^2./a.^(2*H)) .* (exp(-LT.*(1-a.^(-2*N)))-exp(-LT)) ;
        AH2_2 = -(1*H.^2./a.^(1*H)) .* (exp(-LT.*(1-a.^(-1*N)))-exp(-LT));
        AH2 = (AH2_1 + AH2_2) / Rch^2;
        Es2 = AN2 + AKHN + AH2;

        sigma_s2 = Es2 - (Es).^2;
        sigma_s = sqrt(sigma_s2);



%%
%**********************************************************************
%  Approach 1:
%  ALLOWING ZERO-SYMBOL PACKET GENERATION 
%  THEN TREATING ZERO LENGTH PACKETS AS REGULAR PACKETS, 
%  DETERMINISTIC INTER-ARRIVAL DISTRIBUTIONS
%**********************************************************************
    
    %%
    % Using Deterministic Packet Arrival
    %one notation of Kingman formula
    ET = Tv;
    ET2 = Tv.^2;
    sigmaT=0;
    
    ru = Es ./ ET;  
    Ka = sigmaT./ET;
    Ks = sigma_s./Es;
    mu = 1 ./Es;
    Ew = ru .* (Ka.^2 + Ks.^2) ./ (2 * mu .* (1 - ru));
    
    
if adjust_kingman
    Ew  = Ew - ET2 /(2*ET);
end
    
    Ew(find(ru >= 1)) = NaN;
    if newEf
        Ed = Ew + Es + Tv ./2;
    else
        Ed = Ew + Es + Tv ./2; %approx for Ed = Ew + Es + Tv .*(1 + exp(LT))/2;
        %Ed = Ew + Es + Tv .*(1 + exp(LT))/2;
    end
    %Another notation of Kingman formula after manipulation, later check
    %Ew2 = sigma_s .* Es./(2*(Tv - Es));
    %Ew2(find(Tv < Es)) = NaN;
    

    if i == 1
        sim.Res.ET_approx = ET;
        sim.Res.ru_approx = ru;
        sim.Res.Es_approx = Es;
        sim.Res.Es2_approx = Es2;
        sim.Res.Ew_approx = Ew;
        sim.Res.Ed_approx = Ed; %typo: was approax !!!
    end

%%
%**********************************************************************
%  Approach 2:
%  NO ZERO-SYMBOL PACKET GENERATION IS ALLOWED
%  MODIFIED ES AND ES2 BY EXCLUDING EMPTY PACKETS
%  GEMOETRICALLY DISTRIBUITED INTER-ARRIVAL DISTRIBUTIONS
%  PROBLEM: INDEPENDENCY VIOLATION BETWEEN S AND T
%**********************************************************************

    %Modify Es, exclude empty packets
    Es= Es ./ (1- exp(-LT));
    Es2 = Es2 ./ (1- exp(-LT));
    sigma_s2 = Es2 - (Es).^2;
    sigma_s = sqrt(sigma_s2);

    %Time interval is a geometric distribution with p(success) = 1-P0, Pfail = P0
    P0 = exp (-lambda * Tv);
    ET = Tv ./(1-P0);  
    ET2 = (Tv.^2) .* (1+P0)./((1-P0).^2);
    sigmaT2= ET2 - ET.^2;
    sigmaT = sqrt(ET2 - ET.^2);
    %sigmaT_2 = Tv.* sqrt(P0)./((1-P0).^2);

    ru = Es ./ ET;
    Ka = sigmaT ./ ET;
    Ks = sigma_s./Es;
    mu = 1 ./Es;
    
    Ew = (ru./(1-ru)) .* ((Ka.^2 + Ks.^2) / 2) .* Es;
    % Wrong in Check Approx: Ew2 = (ru./(1-ru)) .* ((sigmaT2.^2 + sigma_s2.^2) / 2) .* Es; %=(ru./(1-ru)) .* ((Ka.^2 + Ks.^2) ./ (2* mu)
    
    if adjust_kingman, Ew  = Ew - ET2 /(2*ET);  end
    Ew(find(ru >= 1)) = NaN;
    
    if newEf
        Ef = (Tv/2); % razi changed 20150620, was [Ef = (Tv/2) .*(1+P0)] Consider Ef=Tv/2 according to revised version
    else
        Ef = (Tv/2) .*(1+P0);
    end
    Ed = Ew + Es + Ef;% razi changed  20140731, was Ed = Ew + Es + Tv ./2;
   
    if i == 1  %for given BER
        sim.Res.ET = ET;
        sim.Res.ru = ru;
        sim.Res.Es = Es;
        sim.Res.Es2 = Es2;
        sim.Res.Ew = Ew;
        sim.Res.Ed = Ed;

    else     %error free channel
        sim.ResErrorFree.ET = ET;
        sim.ResErrorFree.ru = ru;
        sim.ResErrorFree.Es = Es;
        sim.ResErrorFree.Es2 = Es2;
        sim.ResErrorFree.Ew = Ew;
        sim.ResErrorFree.Ed = Ed;

    end
end


% pause;
if sim.control.plot_active
%     figure;
%     subplot(211);
%     plot(LT, Es*lambda, 'r');
%     title('Es');
% 
%     subplot(212);
%     plot(Tv*lambda, Es2);
%     title('E(s^2)');
%     pause;
% 
    close all; figure;
    semilogy(Tv*lambda, sim.Ew*lambda,'b');
    hold on;
    semilogy(Tv*lambda, sim.Ew_approx*lambda, 'r:');
    legend('Approach 1', 'Approach 2');
    
    if sim.CodedSystem   %Modifications
        tit=sprintf('Coded Heq:%d  Neq:%d  alpha_eq:%1.4f  Rch:%f Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f] ', H,N,a, sim.ch.Rch, sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD);
    else
        tit = sprintf('uncoded H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f', H,N,a, B, sim.ch.Rch);
    end
    title(tit); pause;

    sim.Res.Ew*lambda,    sim.Res.Ew_approx*lambda, pause;
    %ylim([0, 100000]); 
    % hold on;
end
return