function [] = plot_data_fileNB(dname, fname, save_flag, user_command  ,user_condition)

%default seetings
plot_simulation = 1;
plot_analysis = 1;
adjust_kingman = 0; %might be not necessary or might be done in fanalysis
plot_cognitive=1;

show_data=0;
plot_delays = 1;
plot_all = 0;
plot_analyticals = 0;


%plot more details
plot_Es = 0;
plot_Es2 = 0;
plot_sigmaEs = 0;
plot_coeffEs = 0;
multiSU = 0;
plot_ET = 0;
plot_ET2 = 0;
plot_sigmaT = 0;
plot_coeffT = 0;
plot_Energy = 1 ;
if ~isempty(user_command), eval(user_command); end



% clrs = {'b','g','r','y','c','k','m','b-.','g-.','r-.','y-.','c-.','k-.','m-.','b--','g--','r--','y--','c--','k--','m--',...
%     'b-*','g-*','r-*','y-*','c-*','k-*','m-*','b-x','g-x','r-x','y-x','c-x','k-x','m-x'};
clrs = {'y','c', 'b','g','r','k','m','b-.','g-.','r-.','y-.','c-.','k-.','m-.','b--','g--','r--','y--','c--','k--','m--',...
    'y-*','c-*','b-*','g-*','r-*','k-*','m-*','b-x','g-x','r-x','y-x','c-x','k-x','m-x'};
clrs2 = {'y*','c*','b*','g*','r*','k*','m*','bx','gx','rx','yx','cx','kx','mx'};
markers = {'*', 'o', 'x', 'S', 'd', '^', '<', '>', '.'};
dfnames = strcat(dname,fname); %ali commented
%dfnames = fname;
if iscell(dfnames), nF = length (dfnames); else nF=1; end
close all;

pi = 0;
for i = 1: nF
    if iscell(dfnames), filename = dfnames{i}; else filename = dfnames; end;    %was filename = strcat(dname, readfile(i).name);    
    if isempty(strfind(filename, '.mat')), filename = strcat(filename,'.mat'); end %append.mat
    clearvars sim
    if ~isempty(dir(filename)), load(filename, 'sim'); end  %load file
    if multiSU == 1
     dfnames = fname;
        lpi = 3;
         figure;
%                     plotdef(0);
%                     hold all;
             %       load(filename, 'MultiSU'); 
                for aa=1:4    
             
             S = MultiSU.CRN{1, aa};
             Kv = sim.Kv; if isfield(sim, 'KvA'), KvA=sim.KvA; else KvA=sim.Kv; end
             lambda=sim.lambda; N=sim.N; H=sim.H; a=1/lambda;
             
                   
                    Kinds = [1:length(Kv)];KindsA = [1:length(KvA)];
                     semilogy(Kv(Kinds), S.Davsym(lpi,Kinds));
                     hold on
             
                end
        pause;
    end
    
    if ~exist('sim','var'),    fprintf('Plot Error: File:%s does not exist or corrupted\n', filename); 
        
    else    
        Kv = sim.Kv; if isfield(sim, 'KvA'), KvA=sim.KvA; else KvA=sim.Kv; end
        lambda=sim.lambda; N=sim.N; H=sim.H; a=1/lambda;
        if plot_cognitive && sim.cognitive
            A = sim.Res.A_CRN1;
            if plot_simulation, S = sim.Res.S_CRN; end
        else
            A = sim.Res.A_SN;
            if plot_simulation, S = sim.Res.S_SN; end
        end

        if adjust_kingman
            A.Ew  = A.Ew -KvA./2;
            A.Ed  = A.Ed -KvA./2;
        end
    %    save tp; disp('tp saved'); delay = [A.Ed; S.Davsym], pause;

        plot_simulation = plot_simulation * sim.control.simulation;
        plot_analysis = plot_analysis * sim.control.analysis;

        if show_data

            compEs = [A.Es; S.sav]
            compEs2 = [A.Es2; S.s2av]

            compET = [A.ET; S.Tav]
            compET2 = [A.ET2; S.T2av]


            compW = log([A.Ew; S.wav])
            compD = log([A.Ed; S.Davsym])

    %         compW = log(repmat(Kv/2, [1+size(S.wav,1),1]) + [A.Ew-KvA/2; S.wav])   %kingman adjustment
    %         compD = log(repmat(Kv/2, [1+size(S.wav,1),1]) + [A.Ed-KvA/2; S.Davsym])%kingman adjustment

            pause;
        end


        %fprintf('File:%s  N:%d H:%d BER:%1.3f  Coding \n', filename, sim.N, sim.H, sim.ch.BER);  
        fprintf('File:%s \n Params [sample rate:%f Info rate:%1.3f]  [N:%d] [H:%d]  [Coding Rate D:%1.3f  %1.3f] [Tail Bits D:%1.3f H:%1.3f] \n [BER:%d impN:%d  impH:%d] [Channel Bit rate:%1.3f]\n', ...
        filename, lambda, lambda*N, N, H, sim.coding.CRD, sim.coding.CRH, sim.coding.nTbitsD, sim.coding.nTbitsH, sim.ch.BER, sim.coding.impBERD, sim.coding.impBERH, sim.ch.Rch);
        if sim.cognitive
        fprintf(' Cognitive Params chNotAvRatio:%f  chAvPlenRatio:%f   Ch_AvailMeanTime:%f   Ch_NotAvailMeanTime:%f\n', sim.CRN.chNotAvRatio, sim.CRN.chAvPlenRatio, sim.CRN.Ch_AvailMeanTime, sim.CRN.Ch_NotAvailMeanTime); end



        if sim.CodedSystem   %Modifications
            tit=sprintf('Coded H:%d  N:%d  alpha:%1.4f  Rch:%f Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f] ', sim.H,sim.N,sim.a, sim.ch.Rch, sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD);
        else
            tit = sprintf('Uncoded H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f', sim.H,sim.N,sim.a, sim.ch.BER, sim.ch.Rch);
        end

        if (sim.Framing_mode== 2)
            lpi = 3; %last 50% of simulation

            if nF > 1 %many files, compare
                plotdef(0);  figure(1);  

                if eval (user_condition) 
                    %fprintf('Plot File:%s [%d/%d]\n', filename, pi, nF);%:%s  N:%d H:%d PER:%1.3f  Nsym:%d \n', filename, sim.N, sim.H, sim.ch.Exp_PER, sim.Nsym);  
                    pi = pi+1;
                    if plot_analysis
                        p(pi) = semilogy(KvA, A.Ed, clrs{1+mod(pi,length(clrs))}); hold on; pause;
                        ind = find(floor(KvA)==KvA); if ~isempty(ind), p(pi) = semilogy(KvA(ind), A.Ed(ind), [clrs{1+mod(pi,length(clrs))},markers{1+mod(pi,length(markers))}]); pause; end; 
                        %p(pi) = semilogy(KvA, A.Es, clrs{1+mod(pi,length(clrs))}); hold on; pause;
                        %ind = find(floor(KvA)==KvA); if ~isempty(ind), p(pi) = semilogy(KvA(ind), A.Es(ind), [clrs{1+mod(pi,length(clrs))},markers{1+mod(pi,length(markers))}]); pause; end; 
                    end

                    if plot_simulation
                        p(pi) = semilogy(Kv, S.Davsym(lpi,:), clrs2{1+mod(pi,length(clrs))}); hold on;
                    end
                    %lgnd{pi} = sprintf('N:%d H:%d BER:%d', N, H, sim.ch.BER);   legend(p, lgnd); pause;
                     lgnd{pi} = sprintf('BER:%1.0d', sim.ch.BER);   legend(p, lgnd);%pause;
                   % lgnd{pi} = sprintf('Ch Busy/Av:%1.3f', sim.CRN.chNotAvRatio);   legend(p, lgnd);%pause;
                    xlabel('Framing Parameter: k'); ylabel('End-to-End Delay: E[D]');
                    %xlabel('Framing Parameter: k'); ylabel('Service Time: E[S]');
                    hold on;
                end

            else  %plot a single file
                if plot_all
                    figure;
                    plot(Kv, S.sav(lpi,:), 'b');
                    hold on;
                    plot(Kv, S.wav(lpi,:), 'g');
                    plot(Kv, S.Dav(lpi,:), 'm');
                    plot(Kv, S.Davsym(lpi,:), 'r');
                    legend('Serive Time', 'Waiting', 'Packet Latency', 'Symbol Latency');
                    xlim([.2 5]);
                    if save_flag
                        saveas(gcf, 'TB2', 'fig'); saveas(gcf, 'TB2', 'png');
                    end
                end
                if plot_Es
                    figure;
                    plotdef(2);
                    if plot_analysis
                        plot(KvA, A.Es, 'm-o'); hold on;
                    %   plot(KvA, A.Esll);
                    end
                    if plot_simulation
                        plot(Kv, S.sav(lpi,:), 'c-x'); hold on
                    end

                    xlabel('Packetization Time: T');
                    ylabel('Expected Service Time: E(s)');
                    legend('Analytical','Analytical(PER:0)','Simulation');
                end

                if plot_Es2
                    figure;
                    plotdef(2);
                    if plot_analysis
                        plot(KvA, A.Es2, 'm-o'); hold on;
                     %   plot(KvA, A.Es2ll, 'y'); 
                    end
                    if plot_simulation
                        plot(Kv, S.s2av(lpi,:), 'c-x'); hold on
                    end

                    xlabel('Packetization Time: T');
                    ylabel('E(s^2)');
                    legend('Analytical','Analytical(PER:0)','Simulation');
                end


                if plot_sigmaEs
                    figure;
                    plotdef(2);
                    if plot_analysis
                        plot(KvA, sqrt((A.Es2)-((A.Es).^2)), 'm-o'); hold on;
                    end
                    if plot_simulation
                        sigEs = sqrt(S.s2av(lpi,:)-S.sav(lpi,:).^2);
                        plot(Kv, sigEs, 'c-x'); hold on
                    end

                    xlabel('Packetization Time: T');
                    ylabel('\sigma_s(s)');
                    legend('Analytical','Analytical(PER:0)','Simulation');
                end


                if plot_coeffEs
                    figure;
                    plotdef(0);
                    if plot_analysis
                        plot(KvA, sqrt((A.Es2) - ((A.Es).^2))  ./ A.Es, 'm'); hold on;
                    end
                    if plot_simulation
                        sigEs = sqrt(S.s2av(lpi,:)-S.sav(lpi,:).^2);
                        coeffEs = sigEs./ S.sav(lpi,:);
                        plot(Kv,  coeffEs, 'm:o', 'markersize', 15); hold on
                    end


                    disp('Check PER = 0.2');
                    pause;
                    xlabel('Packetization Time: T');
                    legend('Analytical: PER=0.2','Analytical: PER=0','Simulation:PER=0.2','Simulation:PER=0');
                    ylabel('C(s) = \sigma(S) / E[S]');
                end


                if plot_ET
                    figure;
                    plotdef(2);
                    if plot_analysis
                        plot(KvA, A.ET, 'm-o'); hold on;
                    end
                    if plot_simulation
                        plot(Kv, S.Tav, 'c-x'); hold on
                    end

                    xlabel('Packetization Time: T');
                    ylabel('Expected Inter-packet Arrivals: E(T)');
                    legend('Analytical','Simulation');
                end

                if plot_ET2
                    figure;
                    plotdef(2);
                    if plot_analysis
                        plot(KvA,  A.ET2, 'm-o'); hold on;
                    end
                    if plot_simulation
                        plot(Kv, S.T2av, 'c-x'); hold on
                    end

                    xlabel('Packetization Time: T'); ylabel('E(T^2)');
                    legend('Analytical','Analytical(PER:0)','Simulation');
                end


                if plot_sigmaT
                    figure;
                    plotdef(2);
                    if plot_analysis
                        plot(KvA, sqrt((A.ET2)-((A.ET).^2)), 'm-o'); hold on;
                    end
                    if plot_simulation
                        sigT = sqrt(S.T2av-S.Tav.^2);
                        plot(Kv, sigT, 'c-x'); hold on
                    end

                    xlabel('Packetization Time: T'); ylabel('\sigma(T)');
                    legend('Analytical','Simulation');
                end


                if plot_coeffT
                    figure;
                    plotdef(0);
                    if plot_analysis
                        plot(KvA, sqrt((A.ET2)-((A.ET).^2)) ./ A.ET, 'm'); hold on;
                    end
                    if plot_simulation
                        coeffT = sqrt(S.T2av-S.Tav.^2)./ S.Tav;
                        plot(Kv, coeffT, 'm:o', 'markersize', 15); hold on
                    end


                    xlabel('Packetization Time: T');
                    legend('Analytical: PER=xx','Simulation:PER=xx');
                    ylabel('C(T) = \sigma(T) / E[T]');
                end


                if plot_delays
 
                    figure;
                    plotdef(0);
                    Kinds = [1:length(Kv)];KindsA = [1:length(KvA)];
                    if plot_analysis
                        semilogy(KvA, A.Ef, 'b'); hold on; 
                        semilogy(KvA, A.Es, 'm'); 
                        semilogy(KvA, A.Ew, 'g'); 
                        semilogy(KvA, A.Ed, 'r'); 
                        KindsA = find(isfinite(A.Ew)); minK=KvA(min(KindsA)); maxK=KvA(max(KindsA)); 
                        Kinds = find((Kv >= minK) & (Kv <= maxK));
                    end
                    if plot_simulation
                        semilogy(Kv, S.fav, 'bo', 'linewidth', 1, 'markersize', 8); hold on; 
                        semilogy(Kv, S.sav(lpi,:), 'mx', 'linewidth', 1, 'markersize', 8);  
                        semilogy(Kv(Kinds), S.wav(lpi,Kinds), 'g*', 'linewidth', 1, 'markersize', 8);  
                        semilogy(Kv(Kinds), S.Davsym(lpi,Kinds), 'r.', 'linewidth', 1, 'markersize', 20); 
                    end

                    %xlim([0 5]);    ylim([0 20])
                    xlabel('Framing Parameter: k');
                    ylabel('Expected End to End Delay');
                    %legend('Packetization Delay (Analytical)', 'Waiting Time (Analytical)', 'Average Delay (Analytical)', 'Packetization Delay (Simulation)', 'Waiting Time (Simulation)', 'Average Delay (Simulation)');
                    legend('E[F]:(An)','E[S]:(An)','E[W]:(An)','E[D]:(An)', 'E[F]:(Sim)', 'E[S]:(Sim)', 'E[W]:(Sim)','E[D]: (Sim)');
                    %ylim([0 12]);   xlim([0 4]);
                    %legend('Ew(Analytical)','Ed','Ed/sym');
                    pause;

                    if save_flag
                        saveas(gcf, 'NB_Edsym', 'fig'); saveas(gcf, 'NB_Edsym', 'png');
                    end
                end
               if plot_Energy
                     figure;
                    plotdef(0);
                    Kinds = [1:length(Kv)];KindsA = [1:length(KvA)];
                    if plot_analysis
                    %    subplot(2,1,1);
                        semilogy(KvA, A.Ed, 'r'); 
                        KindsA = find(isfinite(A.Ew)); minK=KvA(min(KindsA)); maxK=KvA(max(KindsA)); 
                        Kinds = find((Kv >= minK) & (Kv <= maxK));
                        hold all
                    end
                    if plot_simulation

                        plot(Kv(Kinds), S.Davsym(lpi,Kinds), 'r.', 'linewidth', 1, 'markersize', 20); 
                    end

                    %xlim([0 5]);    ylim([0 20])
                    xlabel('Framing Parameter: k');
                    ylabel('Expected End to End Delay');
                    %legend('Packetization Delay (Analytical)', 'Waiting Time (Analytical)', 'Average Delay (Analytical)', 'Packetization Delay (Simulation)', 'Waiting Time (Simulation)', 'Average Delay (Simulation)');
                    legend('E[F]:(An)','E[S]:(An)','E[W]:(An)','E[D]:(An)', 'E[F]:(Sim)', 'E[S]:(Sim)', 'E[W]:(Sim)','E[D]: (Sim)');
                    %ylim([0 12]);   xlim([0 4]);
                    %legend('Ew(Analytical)','Ed','Ed/sym');
                 %  subplot(2,1,2);
                   figure; 
                 plot(KvA, A.EE, 'm-o'); hold on;
                    xlabel('Framing Parameter: k');
                    ylabel('Energy');
                    pause
                end
                if plot_analyticals   
                    if sim.CodedSystem   %Modifications
                        tit=sprintf('Coded Heq:%d  Neq:%d  alpha_eq:%1.4f  Rch:%f Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f] ', sim.H,sim.N,sim.a, sim.ch.Rch, sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD);
                    else
                        tit = sprintf('Uncoded H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f', sim.H,sim.N,sim.a, sim.ch.BER, sim.ch.Rch);
                    end

                    if plot_cognitive && sim.cognitive %cognitive network
                            tit =[' CRN ', tit];
                    else %single network
                            tit =[' SN ', tit];
                    end

                    figure; lgnd={'BER:0', ['BER:',num2str(sim.ch.BER)]};
                    subplot(231); semilogy(Kv, A.Es*lambda, 'r-x');  title('E[S]'); legend(lgnd);  %was h=stem(), hb = get(h,'Baseline');set(hb,'Visible','off');set(gca,'yscal','log')
                    subplot(232); semilogy(Kv, A.var_s*lambda, 'r-x');  title('Var[S]');  legend(lgnd);
                    subplot(233); semilogy(Kv, A.var_s ./ (A.Es.^2), 'r-x');  title('Var[S]/E[S]^2'); legend(lgnd);

                    subplot(234); stem(Kv, A.ET*lambda, 'r-x');  title('E[T]'); legend(lgnd);
                    subplot(235); stem(Kv, A.varT*lambda, 'r-x');  title('Var[T]'); legend(lgnd);
                    subplot(236); stem(Kv, A.varT./ (A.ET.^2), 'r-x');  title('Var[T]/E[T]^2'); legend(lgnd);
                    pause;


                    figure;
                    semilogy(KvA, A.Ef*lambda,'b-x'); hold on;
                    semilogy(KvA, A.Es*lambda,'g-x');
                    semilogy(KvA, A.Ew*lambda,'c-x');
                    semilogy(KvA, A.Ed*lambda,'r-x');
                    legend({'E[F]', 'E[S]', 'E[W]', 'E[D]'});
                    title(['Error Ch.: ', tit]); pause;
                    title(tit); pause;
                end

            end  %if nF>1

        end
    end  %file loaded        
end


