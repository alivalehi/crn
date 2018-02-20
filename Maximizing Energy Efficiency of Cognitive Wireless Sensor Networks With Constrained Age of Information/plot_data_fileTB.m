function [] = plot_data_fileTB(dname, fname, save_flag)

plot_simulation = 1;
plot_analysis = 1;

plot_Es = 0;
plot_Es2 = 0;
plot_sigmaEs = 0;
plot_coeffEs = 0;
plot_Ew = 1;
plot_all = 0;

clrs = {'b','g','r','y','c','k','m','b-.','g-.','r-.','y-.','c-.','k-.','m-.','b--','g--','r--','y--','c--','k--','m--',...
    'b-*','g-*','r-*','y-*','c-*','k-*','m-*','b-x','g-x','r-x','y-x','c-x','k-x','m-x'};

if ~isempty(dir),   search_name = strcat(dname,fname); end

if isempty(strfind(search_name, '.mat')), search_name = strcat(search_name,'.mat'); end
readfile = dir(search_name);

nF = length (readfile);
if nF < 1, fprintf('Data file does not exist: %s\n', search_name); end
close all;

pi = 0;
for i = 1: nF
    filename = strcat(dname, readfile(i).name);    load(filename);
    
    plot_simulation = plot_simulation * sim.simulation;
    plot_analysis = plot_analysis * sim.analysis;
    
    fprintf('File:%s  N:%d H:%d PER:%1.3f  \n', filename, sim.N, sim.H, sim.Exp_PER);  
    
    if (sim.Framing_mode== 1) %&& ((sim.N + sim.H== 16+20)||(sim.N+sim.H==8+30)) && (sim.Exp_PER == 0.2)%&& ((sim.H == 20)|| (sim.H == 30)) && (sim.Exp_PER == 0.2)   %time based , (sim.H == 30) && && (sim.N == 16) 
    %if (sim.Framing_mode== 1) && ((sim.N + sim.H== 16+20)||(sim.N+sim.H==8+30)) && (sim.Exp_PER == 0.2)%&& ((sim.H == 20)|| (sim.H == 30)) && (sim.Exp_PER == 0.2)   %time based , (sim.H == 30) && && (sim.N == 16) 
    %if (~isempty(strfind(fname, 'data_TB_D1_17_T7_52_5'))) ||(~isempty(strfind(fname, 'data_TB_D1_17_T3_46_31')))

        lpi = 1; %last 50% of simulation
        if nF > 1 %many files, compare
            plotdef(2);  figure(1);  
            fprintf('Plot File \n');%:%s  N:%d H:%d PER:%1.3f  Nsym:%d \n', filename, sim.N, sim.H, sim.Exp_PER, sim.Nsym);  
            
            lpi = 1; pi = pi+1;
            if plot_analysis
                p(pi) = plot(sim.Tv2 * sim.lambda, sim.Ed* sim.lambda, clrs{1+mod(pi,50)}); hold on;
            end
            
            if plot_simulation
                p(pi) = plot(Tv* sim.lambda, sim.Davsym(lpi,:)* sim.lambda, clrs{1+mod(pi,50)}); hold on;
            end
            lgnd{pi} = sprintf('N:%d H:%d PER:%1.2f', sim.N, sim.H, sim.Exp_PER);  
            hold on;
        else
            if plot_all
                figure;
                plot(Tv, sav(lpi,:), 'b');
                hold on;
                plot(Tv, wav(lpi,:), 'g');
                plot(Tv, Dav(lpi,:), 'm');
                plot(Tv, Davsym(lpi,:), 'r');
                legend('Serive Time', 'Waiting', 'Packet Latency', 'Symbol Latency');
                xlim([.2 5]);
                if save_flag
                    saveas(gcf, 'TB2', 'fig');
                    saveas(gcf, 'TB2', 'pdf');
                end
            end
            if plot_Es
                figure;
                plotdef(2);
                if plot_analysis
                    plot(sim.Tv2, sim.Es, 'm-o'); hold on;
                    plot(sim.Tv2, sim.Es_0, 'g-s'); hold on;
                end
                if plot_simulation
                    plot(Tv, sim.sav(lpi,:), 'm-x'); hold on
                end
                
                xlabel('Packetization Time: T');
                ylabel('Expected Service Time: E(s)');
                legend('Analytical','Analytical(PER:0)','Simulation');
            end
            
            if plot_Es2
                figure;
                plotdef(2);
                if plot_analysis
                    plot(sim.Tv2, sim.Es2, 'm-o'); hold on;
                    plot(sim.Tv2, sim.Es2_0, 'g-s'); hold on;
                end
                if plot_simulation
                    simEs2 = sim.s2av(lpi,:);
                    plot(Tv, simEs2, 'm-x'); hold on
                end
                
                xlabel('Packetization Time: T');
                ylabel('E(s^2)');
                legend('Analytical','Analytical(PER:0)','Simulation');
            end
          
          
            if plot_sigmaEs
                figure;
                plotdef(2);
                if plot_analysis
                    plot(sim.Tv2, sqrt((sim.Es2)-((sim.Es).^2)), 'm-o'); hold on;
                    plot(sim.Tv2, sqrt((sim.Es2_0)-((sim.Es_0).^2)), 'g-s'); hold on;
                end
                if plot_simulation
                    simsigEs = sqrt(sim.s2av(lpi,:)-sim.sav(lpi,:).^2);
                    plot(Tv, simsigEs, 'm-x'); hold on
                end
                
                xlabel('Packetization Time: T');
                ylabel('\sigma_s(s)');
                legend('Analytical','Analytical(PER:0)','Simulation');
            end
            
          
            if plot_coeffEs
                figure;
                plotdef(0);
                if plot_analysis
                    plot(sim.Tv2, sqrt((sim.Es2)-((sim.Es).^2)) ./ sim.Es, 'm'); hold on;
                    plot(sim.Tv2, sqrt((sim.Es2_0)-((sim.Es_0).^2)) ./ sim.Es, 'g'); hold on;
                end
                if plot_simulation
                    simsigEs = sqrt(sim.s2av(lpi,:)-sim.sav(lpi,:).^2);
                    simcoeffEs = simsigEs./ sim.sav(lpi,:);
                    simcoeffEs_0 = sqrt((sim.Es2_0)-((sim.Es_0).^2)) ./ sim.Es;
                    
                    
                    %downsample
                    ind = [1,3, 5, 25, find(mod(sim.Tv * 2, 2) == 0)];
                    Tvds = sim.Tv(ind);
                    simcoeffEs = simcoeffEs(ind);
                    simcoeffEs_0 = simcoeffEs_0(ind);
                    
                    plot(Tvds, simcoeffEs, 'm:o', 'markersize', 15); hold on
                    plot(Tvds, simcoeffEs_0, 'g:*', 'markersize', 15); hold on; %we are sure that there is a perfect match between analythical and theory
                end

                
                disp('Check PER = 0.2');
                pause;
                xlabel('Packetization Time: T');
                legend('Analytical: PER=0.2','Analytical: PER=0','Simulation:PER=0.2','Simulation:PER=0');
                ylabel('K_s = \sigma_s / E[s]');
            end
            
            if plot_Ew
                
                adjust_kingman = 0; %might be not necessary or might be done in fanalysis
                if adjust_kingman
                    sim.Ew  = sim.Ew -sim.Tv2./2;
                    sim.Ed  = sim.Ed -sim.Tv2./2;
                    sim.Ew_approx  = sim.Ew_approx-sim.Tv2./2;
                    
                end
                
                figure;
                plotdef(0);
                if plot_analysis

                    plot(sim.Tv2, sim.Ed, 'r'); hold on;
                    plot(sim.Tv2, sim.Ew_approx, 'g'); 
                    plot(sim.Tv2, sim.Tv2/2, 'b'); 
                end
                if plot_simulation

                    %ind = 1:1:length(Tv);
                    ind  = find((mod(Tv * 20,1)==0)|(Tv < .6));
                    
                    Tvds = Tv(ind);
                    wav = sim.wav(lpi,ind);
                    Dav = sim.Dav(lpi,ind);
                    Davsym = sim.Davsym(lpi,ind);
                    
                    plot(Tvds, Davsym, 'r.-', 'linewidth', 1, 'markersize', 20); hold on
                    plot(Tvds, wav, 'g-*', 'linewidth', 1, 'markersize', 8); hold on

                    ind  = find(mod(Tv * 2,1)==0);
                    Tvds = Tv(ind);
                    Tav = sim.Tav(lpi,ind);
                    plot(Tvds, Tav/2, 'b-o', 'linewidth', 1, 'markersize', 8); hold on
                end
                
                xlim([0 5]);    ylim([0 20])
                xlabel('Packetization Time: T');
                ylabel('Expected End to End Delay');
                legend('Average Delay (Analytical)','Waiting Time (Analytical)','Packetization Delay (Analytical)', 'Average Delay (Simulation)','Waiting Time (Simulation)','Packetization Delay (Simulation)');
                
                
                
                legend('E[f]: (Analytical)','E[w]: (Analytical)','E[d]: (Analytical)', 'E[f]: (Simulation)','E[w]: (Simulation)','E[d]: (Simulation)');
                ylim([0 12]);   xlim([0 4]);
                %legend('Ew(Analytical)','Ed','Ed/sym');
                pause;

                if save_flag
                    saveas(gcf, 'Edsym', 'fig');
                    saveas(gcf, 'Edsym', 'pdf');
                    saveas(gcf, 'Edsym', 'bmp');
                    saveas(gcf, 'Edsym', 'png');
                end
                
            
            end
            

        end

    end
end

if nF > 1 %many files, compare
    figure(1);
    legend(p, lgnd);
    xlim([0 15])
    ylim([0 3e4]);    
    xlabel('Packetization Time: T\lambda');
    ylabel('Average Delay: D\lambda');
end

