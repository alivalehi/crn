if  0
close all
clear all
plotindesx =3;
%if plotindesx == 1
cd C:\Temp\matlab\
%load alldatafinal2.mat
load alldatafinal2.mat
clearvars range midrange selectrd_predict_k meanvar predict_k err
predict_k = KPRED11;
range = (1:2851);
rangetrack(1) = length(range);
midrange(1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(1) = mean(selectrd_predict_k);
varvar(1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(1) = sum((abs((KPRED11(range)-4))>0.5));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  
%err(1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

range = (2852:5848);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-6))>0.5));%sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));

range = (5849:7339);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-7))>0.5));%sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));

range = (7339:10334);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-8))>0.5));%sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));

range = (10334:12789);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-9))>0.5));%sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));

range = (12789:14296);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-8))>0.5));%sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));

range = (14296:16069);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-7))>0.5));%sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));

range = (16069:18477);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-6))>0.5));%sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  
%err(end+1) =sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));

range = (18477:length(predict_k));
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-4))>0.5));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1))
%err(end+1) =sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));;  

errpso = err./rangetrack;
varvarpso =  varvar;
meanvarpso = meanvar;
terrpso = sum(err)/sum(rangetrack);
%plot(predict_k);
%hold on;

clearvars -except  varvarpso meanvarpso errpso terrpso 
%elseif plotindesx == 3
%##############################################################################################################
%
%##############################################################################################################
cd C:\Temp\matlab\
load('alldatamixde.mat','KPRED11'); 
clearvars range midrange selectrd_predict_k meanvar predict_k err varvar
predict_k = KPRED11;
range = (1:1345);
rangetrack(1) = length(range);
midrange(1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(1) = mean(selectrd_predict_k);
varvar(1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(1) = sum((abs((KPRED11(range)-4))>0.5));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  
%err(1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

range = (1345:4240);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-6))>0.5));%sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));
%err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));

range = (4240:7633);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-7))>0.5));%sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (7633:10309);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-8))>0.5));%sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));

range = (10309:14230);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-9))>0.5));%sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  

range = (14230:15951);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-8))>0.5));%sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));

range = (15951:18110);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-7))>0.5));%sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (18110:20610);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-6))>0.5));%sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (20610:length(predict_k));
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-4))>0.5));%sum((abs((KPRED11(range)-4))>0));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  
%err(end+1) = sum((abs((KPRED11(range)-4))>0));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

errde = err./rangetrack;
terrde = sum(err)/sum(rangetrack);
meanvarde = meanvar;
varvarde = varvar;
%plot(predict_k);

clearvars -except  varvarde meanvarde varvarpso meanvarpso errpso errde terrpso terrde terrsa
%elseif plotindesx == 3
%###################################################################################################################################
%
%###################################################################################################################################


cd C:\Temp\matlab\
load ('alldatamixsa.mat ','KPRED11','Kanalplot');

%plot(Kanalplot);
clearvars range midrange selectrd_predict_k meanvar predict_k err varvar
predict_k = KPRED11;
range = (1:2253);
rangetrack(1) = length(range);
midrange(1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(1) = mean(selectrd_predict_k);
varvar(1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(1) = sum((abs((KPRED11(range)-4))>0.5));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  
%err(1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

range = (2253:5249);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-6))>0.5));%sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  

range = (5249:9611);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-7))>0.5));%sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (9611:12620);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-8))>0.5));%sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  
%(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  

range = (12620:17330);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-9))>0.5));%sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  

range = (17330:19630);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-8))>0.5));%sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  

range = (19630:22860);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-7))>0.5));%sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (22860:25140);
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-6))>0.5));%sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  
%err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  

range = (25140:length(predict_k));
rangetrack(end+1) = length(range);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
varvar(end+1) = var(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sum((abs((KPRED11(range)-4))>0.5));%sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1)); 
%err(end+1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

errsa = err./rangetrack;
terrsa = sum(err)/sum(rangetrack);

meanvarsa = meanvar;
varvarsa = varvar;
%plot(predict_k);

figure;

barmat = [4,meanvarpso(1),meanvarde(1),meanvarsa(1);...
          6,meanvarpso(2), meanvarde(2), meanvarsa(2);...
          7,meanvarpso(3), meanvarde(3), meanvarsa(3);...
          8,meanvarpso(4), meanvarde(4), meanvarsa(4);...
          9,meanvarpso(5), meanvarde(5), meanvarsa(5);...
         8,meanvarpso(6), meanvarde(6), meanvarsa(6);...
        7,meanvarpso(7), meanvarde(7), meanvarsa(7);...
     6,meanvarpso(8), meanvarde(8), meanvarsa(8);...
      4,meanvarpso(9), meanvarde(9), meanvarsa(9);...
          ];
%       barmat = [varvarpso(1), varvarsa(1),varvarde(1),4,meanvarpso(1),meanvarde(1),meanvarsa(1);...
%           varvarpso(2), varvarsa(2), varvarde(2),6,meanvarpso(2), meanvarde(2), meanvarsa(2);...
%           varvarpso(3), varvarsa(3), varvarde(3),7,meanvarpso(3), meanvarde(3), meanvarsa(3);...
%           varvarpso(4), varvarsa(4), varvarde(4),8,meanvarpso(4), meanvarde(4), meanvarsa(4);...
%           varvarpso(5), varvarsa(5), varvarde(5),9,meanvarpso(5), meanvarde(5), meanvarsa(5);...
%           varvarpso(6), varvarsa(6), varvarde(6),8,meanvarpso(6), meanvarde(6), meanvarsa(6);...
%           varvarpso(7), varvarsa(7), varvarde(7),7,meanvarpso(7), meanvarde(7), meanvarsa(7);...
%           varvarpso(8), varvarsa(8), varvarde(8),6,meanvarpso(8), meanvarde(8), meanvarsa(8);...
%           varvarpso(9), varvarsa(9), varvarde(9),4,meanvarpso(9), meanvarde(9), meanvarsa(9);...
%           ];
  %   errpso = errpso/ norm(errpso);
   %  errsa = errsa/ norm(errsa);
    % errde = errde/ norm(errde);
      barmat2 = [ errpso(1), errsa(1), errde(1);...
           errpso(2), errsa(2), errde(2);...
           errpso(3), errsa(3), errde(3);...
          errpso(4), errsa(4), errde(4);...
           errpso(5), errsa(5), errde(5);...
           errpso(6), errsa(6), errde(6);...
           errpso(7), errsa(7), errde(7);...
          errpso(8), errsa(8), errde(8);...
           errpso(9), errsa(9), errde(9);...
          ];
      
bar(barmat);
legend('Optimum K','\mu_{PSO}','\mu_{SA}','\mu_{DE}');
hold on
ali = [varvarpso(1), varvarsa(1),varvarde(1);...
    varvarpso(2), varvarsa(2),varvarde(2);...
    varvarpso(3), varvarsa(3),varvarde(3);...
    varvarpso(4), varvarsa(4),varvarde(4);...
    varvarpso(5), varvarsa(5),varvarde(5);...
    varvarpso(6), varvarsa(6),varvarde(6);...
    varvarpso(7), varvarsa(7),varvarde(7);...
    varvarpso(8), varvarsa(8),varvarde(8);...
    varvarpso(9), varvarsa(9),varvarde(9);...
    ];
barmat3 = barmat;
barmat3(:,1) = [];
a1= transpose(barmat3);
a2= transpose(ali);
% xaxisnum = [0.666,0.833, 1.1666, 1.333, 1.666,1.833, 2.1666, 2.333,2.666,2.833, 3.1666, 3.333, 3.666,3.833, 4.1666, 4.333, 4.666,4.833, 5.1666, 5.333, 5.666,5.833, 6.1666, 6.333,6.666,6.833, 7.1666, 7.333, 7.666,7.833, 8.1666, 8.333...
%     8.666,8.833, 9.1666, 9.333]
% xaxisnum = [0.7499,0.9166, 1.0833, 1.250, 1.7499,1.9166, 2.0833, 2.250,2.7499,2.9166, 3.0833, 3.250, 3.7499,3.9166, 4.0833, 4.250, 4.7499,4.9166, 5.0833, 5.250, 5.7499,5.9166, 6.0833, 6.250,6.7499,6.9166, 7.0833, 7.250, 7.7499,7.9166, 8.0833, 8.250...
%     8.7499,8.9166, 9.0833, 9.250];
xaxisnum = [0.9166, 1.0833, 1.278, 1.9166, 2.0833, 2.278,2.9166, 3.0833, 3.278,3.9166, 4.0833, 4.278,4.9166, 5.0833, 5.278,5.9166, 6.0833, 6.278,6.9166, 7.0833, 7.278,7.9166, 8.0833, 8.278...
    8.9166, 9.0833, 9.278];

errorbar (xaxisnum,a1(:),a2(:),'o');
%legend('\sigma_{PSO}','\sigma_{SA}','\sigma_{DE}','Optimum K','\mu_{PSO}','\mu_{SA}','\mu_{DE}');
figure
set(gca,'XTickLabel',{' '})
hold on
%c = categorical({'errorPSO','errorSA','errorDE'});
%bar([terrpso, terrsa, terrde]);
% bar(1, Data, 'colorcode')
% hold on
% bar(2, Data, 'colorcode') , bar (3, Data, 'colorcode')
bar(1,terrpso,'m');
%set(gca,'xticklabel',{'error_{PSO}'})

bar(2,terrsa,'c');
%set(gca,'xticklabel',{'error_{SA}'})
bar(3,terrde,'k');

%bar(sum(barmat2,1));
%set(gca,'xticklabel',{'error_{DE}'})

legend('error_{PSO}','error_{SA}','error_{DE}');

%################################################################################################################################################################################################
end
if 1
close all
clear all
% cd C:\Temp\matlab\
% load alldatafinal.mat
% EE_real_basic = EE_real_track;
% clearvars -except EE_real_MLE EE_real_basic
% EE_real_basic = EE_real_basic(1:714);
% plot(cumsum(EE_real_MLE))
% hold on
% plot(cumsum(EE_real_basic))
% plot(cumsum(EE_real_MLE))
% hold on
% plot(cumsum(EE_real_basic))
%####################################################################
 cd C:\Temp\matlab\
%load alldatafinal2.mat
load alldatafinal2.mat
clearvars range midrange selectrd_predict_k meanvar predict_k err
predict_k = KPRED11;
range = (1:2851);
midrange(1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

range = (2852:5848);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  

range = (5849:7339);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (7339:10334);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  

range = (10334:12789);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  

range = (12789:14296);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  

range = (14296:16069);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (16069:18477);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  

range = (18477:length(predict_k));
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

%err(end+1) = sqrt(sum((predict_k(range)).^2)/(length(data)-1));  
%plot(predict_k);
%hold on;
%errorbar(midrange,meanvar,err,'*');

errorbar(midrange+300,meanvar,err,'o');
hold on;
clear all;
%#############################################################################################################################################
%#
%#############################################################################################################################################
clear all
load alldataMLE.mat
clearvars range midrange selectrd_predict_k meanvar predict_k err
%load alldataMLE.mat
predict_k = KPRED11;
range = (1:2851);
midrange(1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

range = (2852:5848);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  

range = (5849:7339);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (7339:10334);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  

range = (10334:12789);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-9).^2)/(length(predict_k(range))-1));  

range = (12789:14296);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-8).^2)/(length(predict_k(range))-1));  

range = (14296:16069);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-7).^2)/(length(predict_k(range))-1));  

range = (16069:18477);
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-6).^2)/(length(predict_k(range))-1));  

range = (18477:length(predict_k));
midrange(end+1) = round(median(range));
selectrd_predict_k = predict_k(range);
meanvar(end+1) = mean(selectrd_predict_k);
predict_k(range) = ones(1,length(range)).* mean(selectrd_predict_k);
err(end+1) = sqrt(sum((KPRED11(range)-4).^2)/(length(predict_k(range))-1));  

%err(end+1) = sqrt(sum((predict_k(range)).^2)/(length(data)-1));  
%plot(predict_k);
%hold on;
errorbar(midrange-500,meanvar,err,'*');
hold on;
%errorbar(midrange+1000,meanvar,err,'o');
end