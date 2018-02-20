if(1)%new
    clear all
    close all
    Kpredict_array = [];
    BERpredict_array = [];
    upredict_array = [];
    vpredict_array = [];
    EE_real_track_array = [];
    EE_diff_best_track_array = [];
    kAnal_array = [];
    k_array = [];
    EE_prediction_vector1_array = [];
    sim.u = randi([1,15],1,10);%[2 4 6 8];%[3:.5:7];% 4 5 6 7 2 6 12 1 2 4];
    %[.25:.1:4 fliplr(.25:.1:4)];%[1 1 1 1 1 1 1 1 2 3 4 5 6 7
    %8];%randi([4,10],1,sim.m);[2 6 12 1 2 4]
    sim.v = randi([1,15],1,10);%[2 4 6 8];%[4:.5:8];
    rho = (sim.v) ./(sim.u)
    %[fliplr(sim.u) .25:.1:4];%ones(1,length(sim.u))*3;%[8 7 6 5 4 3 2 1 1
    %1 1 1 1 1 1];%randi([4,10],1,sim.m);[1 2 3 4 6 8]
    sim.BER = 1e-3;
    nn = length(sim.u)
    uu = [];
    vv = [];
    k=0;

    fig = figure;
    
    kkk=[];
    nn = [1:4]%[1:10]
    j = 1;
    for i=nn
        load(['result_for' num2str(i) '.mat'])
                vars = {'Kpredict','BERpredict','upredict','vpredict','EE_real_track','EE_diff_best_track','kAnal','k','EE_prediction_vector1'};
                clear(vars{:});
                 load(['result_for' num2str(i) '.mat']);
                %  load(['result_for' num2str(i) '.mat'],'Kpredict','BERpredict','upredict','vpredictend','EE_real_trackend','EE_diff_best_trackend','kAnal','k','EE_prediction_vector1','uuu','vvv');
              %  load(['result_for' num2str(i) '.mat'],'Kpredictend','BERpredictend','upredictend','vpredictend','EE_real_trackend','EE_diff_best_trackend','kAnal','k','EE_prediction_vector1','uuu','vvv');
                Kpredict_array = [Kpredict_array Kpredict];
%                 BERpredict_array = [BERpredict_array BERpredictend];
%                 upredict_array = [upredict_array upredictend];
%                 vpredict_array = [vpredict_array vpredictend];

                EE_real_track_array = [EE_real_track_array EE_real_trackend];
                EE_diff_best_track_array = [EE_diff_best_track_array EE_diff_best_trackend];
                kAnal_array = [kAnal_array kAnal*ones(1,length(Kpredict)-1)];
                k_array = [k_array k];
                EE_prediction_vector1_array = [EE_prediction_vector1_array EE_prediction_vector1];
                uu = [uu uuu(end)];%.*ones(1,length(EE_real_trackend))];
                vv = [vv vvv(end)];%.*ones(1,length(EE_real_trackend))];
                figure
                plot(EE_prediction_vector1),hold on,plot(EE_real_track+EE_diff_best_track)
    title('EE-rho');
                figure;
plot(Kpredict),hold on,plot(.05+(kAnal*ones(1,length(Kpredict)-1)));
title(['optimum k prediction for u=' num2str(uuu(end)) 'v=' num2str(vvv(end))])
legend('predicted optimum k','calculated optimum k');
xlabel('Packet number')
ylabel('Packet framing factor k')
%saveas(gcf,['good result\png\optimum k prediction for u-' num2str(uuu(end)) 'v-' num2str(vvv(end)) '.png'])
%saveas(gcf,['good result\fig\optimum k prediction for u-' num2str(uuu(end)) 'v-' num2str(vvv(end)) '.fig'])
          %      fprintf('u:%f v:%f BER:%f kcalculated:%f kanal:%f\n',upredictend,vpredictend,BERpredictend,k,kAnal);
 %       uuu
%       vvv
%         if(i~=nn(1)&& i ~= nn(end))
%             Kpredict = [Kpredictend Kpredict];
%             x_values = [((j-1)*100):100*j]  ;
%             kkk =[kkk Kpredict];
% %             plot(x_values,Kpredict-1)
%         elseif(i == nn(end))
%         x_values = [((j-1)*100)+1:((j)*100)-15]  ;
%         kkk =[kkk Kpredict(1:(end-25))];
%         else %first
%             x_values = [((j-1)*100)+1:100*j]  ;
%             Kpredict = Kpredict +9;
%             kkk =[kkk Kpredict];
% %             h11 =    plot(x_values,Kpredict-1)
% %             a11 = ['predicted optimum k'];
%         end
%         a = ['\rho = ' num2str(rho(i))];
%         plot(x_values,ones(1,length(x_values))*kAnal,'DisplayName',a);
%         hold on;
%         
%         Kpredictend = Kpredict(end);
%         
% j= j+1;
    end
    uu
    vv
  %  plot(kkk-1,'DisplayName','predicted optimum k');
   % legend
  %  legend('\rho = 3','\rho = .5','\rho = .25','\rho = .33','\rho = 2','\rho = 4','predicted optimum k');
    %legend([h11 h1 h2 h3 h4 h5],[a11 a1 a2 a3 a4 a5])
  %  title('Packet number','Packet framing factor k')
  %  close all
  figure;
  nn = length(EE_real_track_array);
    %     uu = uuu(1:nn);
    %     vv = vvv(1:nn);
    rho= vv./uu;
    B = repmat((rho)',1,length(upredict_array));
    C = reshape(B',length(upredict_array)*length(uu),1);
    %rho = vpredict_array./upredict_array;
  %  plot(vpredict_array./upredict_array);
  %  title('rho');
    figure
    plot(EE_real_track_array),hold on,plot(EE_real_track_array+EE_diff_best_track_array)
    title('EE-rho');
    figure
    plot(Kpredict_array),hold on, plot(kAnal_array);%,plot(rho)
  title(['optimum k prediction for all u-v'])
legend('predicted optimum k','calculated optimum k');
xlabel('Packet number')
ylabel('Packet framing factor k')
saveas(gcf,'good result\png\optimum k prediction for all.png')
saveas(gcf,'good result\fig\optimum k prediction for all.fig')
    %legend('Kpredicted','kAnalytical','rho')
 %   legend('Kpredicted','kAnalytical')
 %   title('K tracking');
end
% if(0)%old
%     clear all
%     Kpredict_array = [];
%     BERpredict_array = [];
%     upredict_array = [];
%     vpredict_array = [];
%     EE_real_track_array = [];
%     EE_diff_best_track_array = [];
%     kAnal_array = [];
%     k_array = [];
%     EE_prediction_vector1_array = [];
%     nn = 6
%     for i=1:nn
%         vars = {'Kpredict','BERpredict','upredict','vpredict','EE_real_track','EE_diff_best_track','kAnal','k','EE_prediction_vector1'};
%         clear(vars{:});
%         load(['result_for' num2str(i) '.mat'],'Kpredict','BERpredict','upredict','vpredict','EE_real_track','EE_diff_best_track','kAnal','k','EE_prediction_vector1','uuu','vvv');
%         Kpredict_array = [Kpredict_array Kpredict];
%         BERpredict_array = [BERpredict_array BERpredict];
%         upredict_array = [upredict_array upredict];
%         vpredict_array = [vpredict_array vpredict];
%         EE_real_track_array = [EE_real_track_array EE_real_track];
%         EE_diff_best_track_array = [EE_diff_best_track_array EE_diff_best_track];
%         kAnal_array = [kAnal_array kAnal];
%         k_array = [k_array k];
%         EE_prediction_vector1_array = [EE_prediction_vector1_array EE_prediction_vector1];
%         fprintf('u:%f v:%f BER:%f kcalculated:%f kanal:%f\n',upredict(end),vpredict(end),BERpredict(end),k,kAnal);
%     end
%     close all
%     nn = length(EE_real_track_array);
%     uu = uuu(1:nn);
%     vv = vvv(1:nn);
%     rho= vv./uu;
%     B = repmat((rho)',1,length(upredict_array));
%     C = reshape(B',length(upredict_array)*length(uu),1);
%     %rho = vpredict_array./upredict_array;
%     plot(vpredict_array./upredict_array);
%     title('rho');
%     figure
%     plot(rho,EE_real_track_array),hold on,plot(rho,EE_real_track_array+EE_diff_best_track_array)
%     title('EE-rho');
%     figure
%     plot(Kpredict_array),hold on, plot(kAnal_array),plot(rho)
%     legend('Kpredicted','kAnalytical','rho')
%     title('K tracking');
% end