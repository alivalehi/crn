%Disp: check why results are totally different than calc_ber_conv for exactly the same setup !!!!


clear; clc; close all; figure; clrs={'b','g','r','m','c','k','y', 'm:', 'g:', 'k:'};
EsNodBv= -3:1:20; 
EsNoVect= 10.^(0.1*EsNodBv);

frameLength = 240;%4000;         % this value must be an integer multiple of 3
targetErrors = 400;%400;
maxNumTransmissions = 5e5;%5e6;

use_QPSK=0;
if use_QPSK
    M=4; EbNoVect= 10.^(0.1*EsNodBv) * log2(M);
else  %BPSK
    M=2; EbNoVect= 10.^(0.1*EsNodBv);
end

figure;
allBERVec = zeros(1,length(EsNodBv));
for scenario =[2]  %1:uncoded    2:conv(1/2)  3:conv(3/4)  5:turbo(1/3)
    fprintf('\nRunning Scenario:%d\n', scenario);
    
    %rng default
    
    if scenario ==5 
        intrlvrIndices = randperm(frameLength);
        hEnc = comm.TurboEncoder('TrellisStructure', poly2trellis(4, [13 15 ],13), 'InterleaverIndices', intrlvrIndices);     %LTE (4, [13 15 17],13)   %Default (4, [13 15],13)
        hDec = comm.TurboDecoder('TrellisStructure', poly2trellis(4, [13 15 ],13), 'InterleaverIndices', intrlvrIndices, 'NumIterations',10);
    elseif scenario == 2 || scenario == 3
        hEnc = comm.ConvolutionalEncoder(poly2trellis(7, [171 133]));
        hDec = comm.ViterbiDecoder(poly2trellis(7, [171 133]), 'InputFormat', 'Unquantized');
        
        if scenario == 3 
            hEnc.PuncturePatternSource = 'Property';  hEnc.PuncturePattern = [1;1;0;1;1;0];  %rate 3/4
            hDec.PuncturePatternSource =  'Property'; hDec.PuncturePattern = hEnc.PuncturePattern;
            hDec.TracebackDepth = 96;
        else
            %hEnc.PuncturePatternSource = 'None';   hDec.PuncturePatternSource = 'None'; %rate 1/2
            hEnc.PuncturePatternSource = 'Property';  hEnc.PuncturePattern = [1;1;1;1;1;1];  %rate 1/2
            hDec.PuncturePatternSource =  'Property'; hDec.PuncturePattern = hEnc.PuncturePattern;
            hDec.TracebackDepth = 40;
        end
    end

    
    hMod = comm.BPSKModulator;
    %hDemod = comm.BPSKDemodulator('DecisionMethod','Log-likelihood ratio');
    hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)', 'SignalPower', 1, 'SamplesPerSymbol', log2(M));
    %hErrorCalc = comm.ErrorRate('ReceiveDelay', hDec.TracebackDepth);


    for snrid = 1:length(EsNodBv)
        if scenario~=1,    reset(hEnc);  reset(hDec); end
        
        EsNo = EsNoVect(snrid); noiseVar = 1/EsNo; EbNo= EsNo*log2(M); fprintf('.');

        hChan.EbNo = EbNoVect(snrid); 
        %hDemod.Variance = noiseVar;
        
        nTransmit=0; nErr=0; 
        while (nErr < targetErrors) && (nTransmit < maxNumTransmissions)
            
            data = randi([0 1], frameLength, 1);

            if scenario == 1,
                encData=data; 
            else
                encData = step(hEnc, data); 
            end
            
            modData = step(hMod, encData);
            channelOutput = step(hChan, modData);
            
            if scenario == 1
                decData =0.5-0.5*sign(real(channelOutput)); 
            else
                decData = step(hDec, real(channelOutput)); 
            end
% fprintf('scenario:%d snrid:%d snr:%f Db EsN0:%f  EbN0:%f  NoiseVar:%f \n', scenario, snrid, EsNodBv(snrid), EsNo, EbNo, noiseVar);
% [encData(1:30)';  modData(1:30)';  channelOutput(1:30)'], pause; 
% [data(1:30)';decData(1:30)'], pause
%             size(data)
%             size(decData), pause;

            nErr= nErr + sum(data~= decData);
            nTransmit =nTransmit + frameLength;
            
        end
        allBERVec(scenario,snrid) = nErr/nTransmit; 
    end
    semilogy(EsNodBv, allBERVec(scenario,:), clrs{scenario}); hold on;
end
legend('Uncoded', 'Conv(1/2)', 'Conv(3/4)', 'Turbo[1/3]');
save result
