% This program is to simulate conv encoder/decoder for various rates to
% plot SNR-BER curves.
% average over multiple runs for more accurate results

clear; close all; clc; 
clrs={'b','g','r','m','c','k','y', 'm:', 'g:', 'k:'};

rv = [1/2, 3/4];
EsNodBv = [-20:0.25:20];
EsNoVect= 10.^(0.1*EsNodBv);
   
%quick run for test
% nRuns = 2;
% frameLength = 2400;        
% targetErrors = 40;
% maxNumTransmissions = 1e4;

nRuns=10; 
frameLength = 4800;  
targetErrors = 400;   
maxNumTransmissions = 1e7; 



if maxNumTransmissions > 1e6, fprintf('This may take very long time ...\n'); pause; end

rnd=floor(100*rand); fname = 'BER_r'; %fname='convBER2_r';
while (~isempty(dir([fname,num2str(rnd),'.mat'])))
    rnd=floor(100*rand);
end    


use_QPSK=0;
if use_QPSK
    M=4; EbNoVect= 10.^(0.1*EsNodBv) * log2(M);error('add demodulator');
else  %BPSK
    M=2; EbNoVect= 10.^(0.1*EsNodBv);
end 



alluncodedBERVec=cell(nRuns,1);
hMod = comm.BPSKModulator;
hChan = comm.AWGNChannel('NoiseMethod', 'Signal to noise ratio (Eb/No)', ...
  'SignalPower', 1, 'SamplesPerSymbol', log2(M));

fprintf('Uncoded System ...\n');
alluncodedBERVec  = zeros(nRuns,length(EbNoVect)); 
for nrun = 1 : nRuns
    fprintf('\nrun %d/%d: ', nrun,nRuns);
    for n=1:length(EbNoVect)
        fprintf('.');
        hChan.EbNo = EbNoVect(n); 
        nErr=0; nTransmit=0;
        while (nErr < targetErrors) && (nTransmit < maxNumTransmissions)
            data = randi([0 1], frameLength, 1);
            modData = step(hMod, data);
            channelOutput = step(hChan, modData);
            decData = 0.5-0.5 * sign(real(channelOutput));
            nErr=nErr + sum(data~=decData);
            nTransmit=nTransmit+frameLength;
        end
        alluncodedBERVec(nrun,n) = nErr/nTransmit;
    end
end
uncodedBERVec=mean(alluncodedBERVec);



%%
allconvBERVec=cell(nRuns, length(rv));
for rvid = 1:length(rv)
    fprintf('\n Encoder Type :%d/%d\n',rvid, length(rv));
    hConvEnc = comm.ConvolutionalEncoder(poly2trellis(7, [171 133]));
    if rvid == 2
        hConvEnc.PuncturePatternSource = 'Property';
        hConvEnc.PuncturePattern = [1;1;0;1;1;0];  %rate 3/4
    elseif rvid==1   
        hConvEnc.PuncturePatternSource = 'Property';
        hConvEnc.PuncturePattern = [1;1;1;1;1;1];  %rate 1/2, no puncturing
    end

    hVitDec = comm.ViterbiDecoder(poly2trellis(7, [171 133]), ...
      'InputFormat', 'Unquantized');
    hVitDec.PuncturePatternSource =  'Property';
    hVitDec.PuncturePattern = hConvEnc.PuncturePattern;

    if rvid == 2
        hVitDec.TracebackDepth = 96;
    else
        hVitDec.TracebackDepth = 40;
    end
    hErrorCalc = comm.ErrorRate('ReceiveDelay', hVitDec.TracebackDepth);

    for nrun = 1 : nRuns
        fprintf('\nrun %d/%d: ', nrun,nRuns);
        convBERVec = zeros(3,length(EbNoVect)); 
        for n=1:length(EbNoVect)
          fprintf('.');
          reset(hErrorCalc)
          reset(hConvEnc)
          reset(hVitDec)
          hChan.EbNo = EbNoVect(n); % Set the channel EbNo value for simulation
          while (convBERVec(2,n) < targetErrors) && (convBERVec(3,n) < maxNumTransmissions)
            data = randi([0 1], frameLength, 1);
            encData = step(hConvEnc, data);
            modData = step(hMod, encData);
            channelOutput = step(hChan, modData);
            decData = step(hVitDec, real(channelOutput));
            convBERVec(:,n) = step(hErrorCalc, data, decData);
          end

        end
        allconvBERVec{nrun, rvid}=convBERVec;
    end
end

clearvars convBERVec
convBERVec(1:length(rv),1)={0};
for nrun=1:nRuns
    for rvid = 1: length(rv)
        convBERVec{rvid} = convBERVec{rvid} + allconvBERVec{nrun, rvid}(1,:);
    end
end
convBERVec{rvid} = convBERVec{rvid}/nRuns;


save ([fname,num2str(rnd),'.mat'], 'rv', 'EsNodBv', 'use_QPSK', 'alluncodedBERVec', 'allconvBERVec', 'uncodedBERVec','convBERVec');

if 1
figure;
for nrun=1:nRuns
    semilogy(EbNoVect, uncodedBERVec, 'b'); hold on;
    semilogy(EbNoVect, convBERVec{1}(1,:), 'g'); 
    semilogy(EbNoVect, convBERVec{2}(1,:), 'r'); 
end
legend('Uncoded', 'Conv[1]', 'Conv[2]');
end