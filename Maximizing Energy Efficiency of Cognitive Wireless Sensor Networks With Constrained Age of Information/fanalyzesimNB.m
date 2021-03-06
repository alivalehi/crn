function [sim] = fanalyzesimNB(sim)
%This program gets the simulation paramteres and then calculates
%(analythicaly) the expected delay [for number based system]
% step 1: calculateing E(s) and E(s^2)
% step 2: using kingman formula

adjust_kingman = 1;

if nargin < 1
    %default vars
    lambda = 10; N = 3;  H = 10;  BER=1e-5; Rch=1;
    Kv = [1:1:10, 20:5: 50];
    Plenv = Kv * N + H;
    PERv = 1 - (1 - BER) .^ Plenv;
else
    H = sim.H;
    N = sim.N;
    lambda = sim.lambda;
    Kv = sim.KvA;
    
    %these values are not used
    if sim.CodedSystem
        sim.Plenv = Kv * (sim.N / sim.coding.CRD) + (sim.H / sim.coding.CRH) + sim.coding.nTbitsD + sim.coding.nTbitsH; %Expected average packet length for each T
        sim.PERv = 1 - ((1 - sim.ch.BER/sim.coding.impBERD) .^ (Kv * (sim.N / sim.coding.CRD))+ sim.coding.nTbitsD)*((1 - sim.ch.BER/sim.coding.impBERH) .^ ((sim.H / sim.coding.CRH) + sim.coding.nTbitsH));
    else
        sim.Plenv = Kv * sim.N + sim.H;
        sim.PERv = 1 - (1 - sim.ch.BER) .^ (sim.Plenv);
    end
    
    if (sim.Framing_mode ~= 2) || (sim.input_process ~= 1)
        error('Err: Invalid params');
    end
    
    BER = sim.ch.BER;
end
%sim.KvA, sim.KvA, sim.PERv,sim.Plenv,disp('del1'); pause;

H0=H; N0=N;
B0 = BER;
a0 = 1 - B0;
if sim.CodedSystem   %Modifications
    %find equivalent alpha, beta, H and N in order to use the same equations
    
    H = (H0 / sim.coding.CRH + sim.coding.nTbitsD + sim.coding.nTbitsH); if length(H)==1, H=repmat(H, [1,length(Kv)]); end  %equivalent H to use equations for uncoded system
    N = (N0 / sim.coding.CRD); if length(N)==1, N=repmat(N, [1,length(Kv)]); end
    Lenv = Kv.*N+H;
    PERv = 1 - ((1 - B0/sim.coding.impBERD) .^ (Kv * sim.N / sim.coding.CRD + sim.coding.nTbitsD)) .* ((1 - B0/sim.coding.impBERH) .^ (sim.H / sim.coding.CRH + sim.coding.nTbitsH));  %(1-berData)^nDatabits * (1-berHeader)^nHeaderbits
    a = (1-PERv) .^ (1./(H0 + N0 .* Kv)); B=1-a;
    
    %PERv,Lenv,disp('del2'); pause;
    fprintf('Analytical Results: Parameters N:[%d==>%1.1f]  H:[%d==>%1.1f] Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f]   BER:[%1.4f==>%s]   \n', N0, N(1), H0, H(1), sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD, B0, num2str(B([1:min(3,end), max(1,end-1), end])));
    
else
    a = a0; B=B0;
    Lenv = Kv.*N+H;
    PERv = 1-(1 - B).^ Lenv;
    %PERv,Lenv,disp('del2'); pause;
    fprintf('Analytical Results: H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f\n', H,N,a, B, sim.ch.Rch);
end


SN=[];
SN.Lenv_anal = Lenv;
SN.PERv_anal=PERv;
SN.Ef= 0.5*(Kv-1)./ sim.lambda; %packet formation time for number based policy

%moments of inter-arrival time for Gamma(shape=k, rate=lambda, scale=1/lambda)
SN.ET = Kv / lambda;
SN.ET2 = (Kv .*(Kv+1)) ./ (lambda.^2);
SN.varT=Kv / lambda^2;

%moments of service time :L/Rch * gemoetric(Ps = 1-PER = (1-b)^L= (1-b)^(KN+H)
Ps= 1-PERv;
Rch = sim.ch.Rch;

s1=Lenv./Rch;
SN.Es = (s1) ./ Ps;  %Geometric distributed
SN.Es2 = ((s1).^2) .*  ((2-Ps)./Ps.^2);
SN.var_s = ((s1).^2) .*  (PERv./Ps.^2);


SN.ru = SN.Es ./ SN.ET;
SN.CT2 = SN.varT./(SN.ET.^2);
SN.Cs2 = SN.var_s./(SN.Es.^2);
%one notation of Kingman formula, later check
%save temp, pause;

SN.Ew = 0.5 * (SN.ru ./ (1-SN.ru)) .* (SN.CT2 + SN.Cs2) .* SN.Es;  % % % % %     if adjust_kingman, Ew  = Ew - ET2 /(2*ET);  end
SN.Ew(find(SN.ru >= 1)) = NaN;

%Total Delay
SN.Ed = SN.Ew + SN.Es + SN.Ef;

%SN.Es, SN.ET, SN.ru, SN.Ed ,disp('del3'); pause;

CRN1 = [];
if sim.cognitive   %calculate service time and waiting time for partially available channel
    if sim.CRN.Approx == 1   %use approximation of chAvailTime,chNotAvailTime >> Plen/Rch
        %Test provided in test_ES_CRN
        
        v = sim.CRN.Ch_AvailMeanTime; u = sim.CRN.Ch_NotAvailMeanTime; s=s1;
        
        %moments for one copy transmission time (PER=0)
        
        ESk = (v./(v+u)) .* s1   + (u./(v+u)) .* (s1 + u/2);  %ok
        ESk2 = (v./(v+u)) .* s1.^2 + (u./(v+u)) .* (s1.^2 + 2*u.^2/3 + s1.*u);  %ok
        varsk =(u.^3) .* (8*v+5*u)./(12 *((v+u).^2));   %ok
        varsk2 = ESk2-ESk.^2;
        if abs(varsk-varsk2)/abs(varsk+varsk2) > 0.01, warning('Potential error in calculating var(CRN:S) Err:%f<>%f \n', varsk, varsk2); end
        
        CRN1.Es = ESk ./ Ps;
        
        
        CRN1.Es2   = ESk2./Ps  + (ESk.^2) .* (2*(1-Ps)./(Ps.^2));
        
        CRN1.var_s = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));
        
        if sim.control.debug_active
            %paper notations
            ES_2 =(1 ./ Ps) * (s1  + (u.^2 ./(2*(v+u))));  %another notation
            ES2_2 = (1./Ps).*((v./(v+u)) .* s1.^2 + (u./(v+u)) .* (s1.^2 + 2*u.^2/3  + u.* s1)) +  ...
                (2*(1-Ps)./(Ps.^2)) .* (s1 + (u.^2 ./(2*(u+v)))).^2;
            varS_2 = CRN1.Es2 - CRN1.Es.^2;   %for extr check
            %closed form
            varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
                (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
            %different notation
            varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
                ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));
            
            %coefficient of variation
            CS = sqrt(CRN1.var_s) ./ CRN1.Es;
            CS_2 = sqrt( 1 - 2*Ps + Ps *(s1.^2 + (u./(v+u)) .* (2 * u.^2 /3 + u.*s1))/((s1 + u.^2/(2*(u+v))).^2));
            
            
            %check various notations
            if mean(abs(CRN1.Es-ES_2)./abs(CRN1.Es+ES_2)) > 0.01
                warning('Potential error in calculating E(CRN:S) Err:%f<>%f\n', mean(CRN1.Es), mean(ES_2));
            end
            
            if mean(abs(CRN1.Es2-ES2_2)./abs(CRN1.Es2+ES2_2)) > 0.01
                warning('Potential error in calculating E(CRN:S2) Err:%f<>%f\n', mean(CRN1.Es2), mean(ES2_2));
            end
            
            if mean(abs(CRN1.var_s-varS_2)./abs(CRN1.var_s+varS_2)) > 0.01 || ...
                    mean(abs(CRN1.var_s-varS_3)./abs(CRN1.var_s+var_3)) > 0.01 || ...
                    mean(abs(CRN1.var_s-varS_4)./abs(CRN1.var_s+varS_4)) > 0.01
                warning('Potential error in calculating var(CRN:S) Err:%f<>%f<>%f<>%f \n', mean(CRN1.var_s), mean(CRN1.var_s_2), mean(CRN1.var_s3), mean(CRN1.var_s4));
            end
            
        end
        
        CRN1.ET = SN.ET; CRN1.ET2 = SN.ET2; CRN1.varT = SN.varT; CRN1.CT2=SN.CT2;  CRN1.Ef=SN.Ef; %remaind unchanged
        CRN1.ru = CRN1.Es ./ CRN1.ET;
        CRN1.Cs2 = CRN1.var_s./(CRN1.Es.^2);
        CRN1.Ew = 0.5 * (CRN1.ru ./ (1-CRN1.ru)) .* (CRN1.CT2 + CRN1.Cs2) .* CRN1.Es;  % % % % %     if adjust_kingman, Ew  = Ew - ET2 /(2*ET);  end
        CRN1.Ew(find(CRN1.ru >= 1)) = NaN;
        CRN1.Ed = CRN1.Ew + CRN1.Es + CRN1.Ef;
        
        %save tp3; disp('tp3 saved in fanalyzesimNB'); pause;
    elseif sim.CRN.Approx == 0   %use approximation of chAvailTime,chNotAvailTime >> Plen/Rch
        %%%% ali added
        P_s = 10;
        P_p = 60;
        Plenv = Kv * N + H;
        CRN1.EE = (((1 - BER) .^( -Plenv).*(Kv.*N+H).*(P_s))+(P_p))./Kv.*N;
        %plot(CRN1.EE);
        %##################################################################################################################################################################################################################
        if(sim.trafficmode ==1) %%heavy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Test provided in test_ES_CRN
            
            %    v = ones(1,length(sim.CRN.Ch_AvailMeanTime)).*sim.CRN.Ch_AvailMeanTime(1);
            %     u = ones(1,length(sim.CRN.Ch_NotAvailMeanTime)).*sim.CRN.Ch_NotAvailMeanTime(1);
            v = sim.CRN.Ch_AvailMeanTime; u = sim.CRN.Ch_NotAvailMeanTime;
            %  u = sim.CRN.chNotAvRatio;
            %  v = sim.CRN.chAvPlenRatio;
            s=s1;
            bb = 1;
            %moments for one copy transmission time (PER=0)
            
            p              = exp(-s1./v);
            v1             = ((((-v.*s1-1).*exp(-v.*s1))+1)./v)./(1-exp(-s1.*v));
            %  v1             = ((((-v.*s1-1).*exp(-s1./v))+1)./v)./(1-exp(-s1./v));
            E_uv_2         = 2.*u.^2+((-(s1.^2+2.*s1.*v+2.*v.^2).*exp(-s1./v)+2.*v.^2)./(1-exp(-s1./v)))+2.*u.*v1;
            a              = 1;
            ES_BusyArrive  = ((1./p)-a).*((u)+(v1));
            ES_AvArrive    = ((1./p).*(u))+((((1./p)-a)).*(v1));
            E_bu2          = E_uv_2 .*((1./p)-a)+(((2-p)./p.^2)-((1+2.*a)./p)+a.^2+a).*(u+v1).^2;
            ES_BusyArrive2 = 2.*u.^2./3 + E_bu2+s1.^2 +2.*s1.*(u./2+ES_BusyArrive)+ 2.*ES_BusyArrive.*u./2;
            ES_AvArrive2   = ES_BusyArrive2+2.*u.^2+2.*u.*ES_BusyArrive;
            E_sav2         = E_bu2+(2.*u.^2)+2.*u.*((1./p)-a).*(u+v1);
            var_bu         = E_bu2-ES_BusyArrive.^2;
            var_sav        = E_sav2-ES_BusyArrive.^2; % it should be (E_sav2- ES_AvArrive^2);
            %               state1 = exp(-sv1/v)*v/(u+v);
            %               state3 = (1-exp(-sv1/v))*v/(u+v);
            %               state2 = b/(a+b);
            %            transition = [exp(-sv1/v) exp(-sv1/v) exp(-sv1/v);0 0 0;1-exp(-sv1/v) 1-exp(-sv1/v) 1-exp(-sv1/v)]^1e4;
            %
            %               transition = [exp(-sv1/a) 0 1-exp(-sv1/a);exp(-sv1/a) 0 1-exp(-sv1/a);0 1 0]^1e4
            %              state = transition*[state1;state2;state3]
            %               if()
            %                     Scenario1 = (p.* v./(v+u)) .* s1; %arrive at available time slot and it has enough time to do service
            %                     Scenario2 =(u./(v+u)).* (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
            %                     Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
            %               else
            Scenario1 = (p) .* s1; %arrive at available time slot and it has enough time to do service
            Scenario2 =.3*(u./(v+u)).* (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
            alpha = u./v;
            ali = (s1 + u./2+ ES_BusyArrive);
            % Scen2approx = ((((alpha.*v.^2)./2)+1).*exp(-s1./v))+ ((((alpha.*v.^2))+1).*exp(s1./v))+(-1.5.*(v.^2)-2);
            Scen2approx = 1.5 .*exp(-s1)+2.*exp(s1)-3.5;
            %######################################after revision
            %    Scenario3 = ((1-p)).* (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
            alpha = ((1-exp(-s1./v)+(s1./v).*expint(s1./v)));
            phi = v1;%(v./2-(1./2).*(s1+v).*exp(-s1./v)+((s1.^2)./(2.*v)).*expint(s1./v))./alpha;
            Scenario3 = ((1-p)).* (ES_AvArrive+s1+phi);
            %#######################################
            
            %  end
            %  ESk = Scenario1+Scenario2+Scenario3;
            v12 = ((2.*(1./v).^2+(-s1.^2-2.*(1./v).*s1-2.*(1./v).^2).*exp(-s1.*v)))./(1-exp(-s1.*v));
            ESk = Scenario1+Scenario3;
            Scenario1_2 = (((p).^2) .* s1.^2);
            Scenario2_2 = .3*((u./(v+u)).^2.* (s1.^2 + 2.*u.^2/3 + ES_BusyArrive2 + 2.*s1.*(u/2 +ES_BusyArrive)+2.*u/2.*ES_BusyArrive));
            %Scenario3_2 = ((1-p)).^2.* (s1.^2 + 2.*v1.^2/3 + ES_AvArrive2 + 2.*s1.*(v1/2+ES_AvArrive)+2.*v1/2.*ES_AvArrive);
            Scenario3_2 = ((1-p).^2).* (s1.^2 + v12 + ES_AvArrive2 + 2.*s1.*(v1+ES_AvArrive)+2.*v1.*ES_AvArrive);
            ESk2 =Scenario1_2+Scenario3_2+2.*Scenario1.*(Scenario3);
            % ESk2 =Scenario1_2+Scenario2_2+Scenario3_2+2.*Scenario1.*(Scenario2+Scenario3)+2.*Scenario2.*Scenario3;
            varsk2 = ESk2- ESk.^2;   %ok
%                       v12 = ((-2.*v.^2+(-s1.^2+-2.*v.*s1-2.*v.^2).*exp(-s1./v))./v)./(1-exp(-s1.*v));
%             ESk = Scenario1+Scenario3;
%             Scenario1_2 = ((p).^2 .* s1.^2);
%             Scenario2_2 = .3*((u./(v+u)).^2.* (s1.^2 + 2.*u.^2/3 + ES_BusyArrive2 + 2.*s1.*(u/2 +ES_BusyArrive)+2.*u/2.*ES_BusyArrive));
%             %Scenario3_2 = ((1-p)).^2.* (s1.^2 + 2.*v1.^2/3 + ES_AvArrive2 + 2.*s1.*(v1/2+ES_AvArrive)+2.*v1/2.*ES_AvArrive);
%             Scenario3_2 = ((1-p)).^2.* (s1.^2 + v12 + ES_AvArrive2 + 2.*s1.*(v1+ES_AvArrive)+2.*v1.*ES_AvArrive);
%             ESk2 =Scenario1_2+Scenario3_2+2.*Scenario1.*(Scenario3);
%             % ESk2 =Scenario1_2+Scenario2_2+Scenario3_2+2.*Scenario1.*(Scenario2+Scenario3)+2.*Scenario2.*Scenario3;
%             varsk2 = ESk2- ESk.^2;   %ok  
            
            % %%%%light
            %  ESkll = Scenario1 +Scenario2+Scenario3;
            %             Scenario1_2ll = ((p.* v./(v+u)).^2 .* s1.^2);
            %             Scenario2_2ll = ((u./(v+u)).^2.* (s1.^2 + 2.*u.^2/3 + ES_BusyArrive2 + 2.*s1.*(u/2 +ES_BusyArrive)+2.*u/2.*ES_BusyArrive));
            %             Scenario3_2ll = ((v./(v+u)).*(1-p)).^2.* (s1.^2 + 2.*v1.^2/3 + ES_AvArrive2 + 2.*s1.*(v1/2+ES_AvArrive)+2.*v1/2.*ES_AvArrive);
            %             ESk2ll =Scenario1_2+Scenario2_2+Scenario3_2+2.*Scenario1.*(Scenario2+Scenario3)+2.*Scenario2.*Scenario3;
            %  CRN1.Esll = ESkll ./ Ps;
            %
            %
            %                 CRN1.Es2ll   = ESk2ll./Ps  + (ESkll.^2) .* (2*(1-Ps)./(Ps.^2));
            %
            %                 CRN1.var_sll = ESk2ll./Ps  + (ESkll.^2) .* ((1-2*Ps)./(Ps.^2));
            %                                 CRN1.ETll = SN.ET; CRN1.ET2ll = SN.ET2; CRN1.varTll = SN.varT; CRN1.CT2ll=SN.CT2;  CRN1.Efll=SN.Ef; %remaind unchanged
            %                 CRN1.rull = CRN1.Esll ./ CRN1.ETll;
            %                 CRN1.Cs2ll = CRN1.var_sll./(CRN1.Esll.^2);
            %                 CRN1.Ewll = 0.5 * (CRN1.rull ./ (1-CRN1.rull)) .* (CRN1.CT2ll + CRN1.Cs2ll) .* CRN1.Esll;  % % % % %     if adjust_kingman, Ew  = Ew - ET2 /(2*ET);  end
            %                 CRN1.Ewll(find(CRN1.rull >= 1)) = NaN;
            %                 CRN1.Edll = CRN1.Ewll + CRN1.Esll + CRN1.Efll;
            % %%%%
            
            CRN1.Es = ESk ./ Ps;
            
            
            CRN1.Es2   = ESk2./Ps  + (ESk.^2) .* (2*(1-Ps)./(Ps.^2));
            
            CRN1.var_s = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));
            
            if sim.control.debug_active
                %paper notations
                ES_2 =(1 ./ Ps) .* (s1  + (u.^2 ./(2.*(v+u))));  %another notation
                ES2_2 = (1./Ps).*((v./(v+u)) .* s1.^2 + (u./(v+u)) .* (s1.^2 + 2*u.^2/3  + u.* s1)) +  ...
                    (2*(1-Ps)./(Ps.^2)) .* (s1 + (u.^2 ./(2*(u+v)))).^2;
                varS_2 = CRN1.Es2 - CRN1.Es.^2;   %for extr check
                %closed form
                varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
                    (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
                %different notation
                varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                    (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
                    ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));
                
                %coefficient of variation
                CS = sqrt(CRN1.var_s) ./ CRN1.Es;
                CS_2 = sqrt( 1 - 2.*Ps + Ps .*(s1.^2 + (u./(v+u)) .* (2 .* u.^2 /3 + u.*s1))/((s1 + u.^2/(2.*(u+v))).^2));
                
                
                %check various notations
                if mean(abs(CRN1.Es-ES_2)./abs(CRN1.Es+ES_2)) > 0.01
                    warning('Potential error in calculating E(CRN:S) Err:%f<>%f\n', mean(CRN1.Es), mean(ES_2));
                end
                
                if mean(abs(CRN1.Es2-ES2_2)./abs(CRN1.Es2+ES2_2)) > 0.01
                    warning('Potential error in calculating E(CRN:S2) Err:%f<>%f\n', mean(CRN1.Es2), mean(ES2_2));
                end
                
                if mean(abs(CRN1.var_s-varS_2)./abs(CRN1.var_s+varS_2)) > 0.01 || ...
                        mean(abs(CRN1.var_s-varS_3)./abs(CRN1.var_s+var_3)) > 0.01 || ...
                        mean(abs(CRN1.var_s-varS_4)./abs(CRN1.var_s+varS_4)) > 0.01
                    warning('Potential error in calculating var(CRN:S) Err:%f<>%f<>%f<>%f \n', mean(CRN1.var_s), mean(CRN1.var_s_2), mean(CRN1.var_s3), mean(CRN1.var_s4));
                end
                
            end
            
            CRN1.ET = SN.ET; CRN1.ET2 = SN.ET2; CRN1.varT = SN.varT; CRN1.CT2=SN.CT2;  CRN1.Ef=SN.Ef; %remaind unchanged
            CRN1.ru = CRN1.Es ./ CRN1.ET ;
            CRN1.ru(find(CRN1.ru >= 1))=.99; %test after revision
            %                figure;
            %                 plot(CRN1.ru);
            %                 title(['H is ',num2str(sim.H),'N is ',num2str(sim.N),])
            
            CRN1.Cs2 = CRN1.var_s./(CRN1.Es.^2);
            CRN1.Ew = 0.5 * (CRN1.ru ./ (1-CRN1.ru)) .* (CRN1.CT2 + CRN1.Cs2) .* CRN1.Es;  % % % % %     if 1, Ew  = Ew - ET2 /(2*ET);  end
          
            % CRN1.Ew(find(CRN1.ru >= 1)) = 0;
            CRN1.Ed = CRN1.Ew + CRN1.Es + CRN1.Ef;
            
            %save tp3; disp('tp3 saved in fanalyzesimNB'); pause;
            %############light#################################################################################################################################################################################################
        else
            %Test provided in test_ES_CRN
            
            v = sim.CRN.Ch_AvailMeanTime; u = sim.CRN.Ch_NotAvailMeanTime; s=s1;
            bb = 1;
            %moments for one copy transmission time (PER=0)
            
            p              = exp(-s1./v);
            v1             = -exp(-s1./v).*(v.*(v+s1)./(v.*(1-exp(-s1./v))))+v.^2./(v.*(1-exp(-s1./v)));
            E_uv_2         = 2.*u.^2+((-(s1.^2+2.*s1.*v+2.*v.^2).*exp(-s1./v)+2.*v.^2)./(1-exp(-s1./v)))+2.*u.*v1;
            a              = 1;
            ES_BusyArrive  = ((1./p)-a).*((u)+(v1));
            ES_AvArrive    = ((1./p).*(u))+((((1./p)-a)).*(v1));
            E_bu2          = E_uv_2 .*((1./p)-a)+(((2-p)./p.^2)-((1+2.*a)./p)+a.^2+a).*(u+v1).^2;
            ES_BusyArrive2 = 2.*u.^2./3 + E_bu2+s1.^2 +2.*s1.*(u./2+ES_BusyArrive)+ 2.*ES_BusyArrive.*u./2;
            ES_AvArrive2   = ES_BusyArrive2+2.*u.^2+2.*u.*ES_BusyArrive;
            E_sav2         = E_bu2+(2.*u.^2)+2.*u.*((1./p)-a).*(u+v1);
            var_bu         = E_bu2-ES_BusyArrive.^2;
            var_sav        = E_sav2-ES_BusyArrive.^2; % it should be (E_sav2- ES_AvArrive^2);
            %               state1 = exp(-sv1/v)*v/(u+v);
            %               state3 = (1-exp(-sv1/v))*v/(u+v);
            %               state2 = b/(a+b);
            %            transition = [exp(-sv1/v) exp(-sv1/v) exp(-sv1/v);0 0 0;1-exp(-sv1/v) 1-exp(-sv1/v) 1-exp(-sv1/v)]^1e4;
            %
            %               transition = [exp(-sv1/a) 0 1-exp(-sv1/a);exp(-sv1/a) 0 1-exp(-sv1/a);0 1 0]^1e4
            %              state = transition*[state1;state2;state3]
            %               if()
            %                     Scenario1 = (p.* v./(v+u)) .* s1; %arrive at available time slot and it has enough time to do service
            %                     Scenario2 =(u./(v+u)).* (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
            %                     Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
            %               else
            Scenario1 = (p.* v./(v+u)) .* s1; %arrive at available time slot and it has enough time to do service
            Scenario2 =(u./(v+u)).* (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
            %######################################after revision
            %Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
            alpha = ((1-exp(-s1./v)+(s1./v).*expint(s1./v)));
            phi = (v./2-(1./2).*(s1+v).*exp(-s1./v)+((s1.^2)./(2.*v)).*expint(s1./v))./alpha;
            Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+phi./2);%arrive at available time slot but it doesnt have  enough time to do service
            %#######################################
            Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+v1./2);%arrive at available time slot but it doesnt have  enough time to do service
            %  end
            ESk = Scenario1 +Scenario2+Scenario3;
            Scenario1_2 = ((p.* v./(v+u)).^2 .* s1.^2);
            Scenario2_2 = ((u./(v+u)).^2.* (s1.^2 + 2.*u.^2/3 + ES_BusyArrive2 + 2.*s1.*(u/2 +ES_BusyArrive)+2.*u/2.*ES_BusyArrive));
            Scenario3_2 = ((v./(v+u)).*(1-p)).^2.* (s1.^2 + 2.*v1.^2/3 + ES_AvArrive2 + 2.*s1.*(v1/2+ES_AvArrive)+2.*v1/2.*ES_AvArrive);
            ESk2 =Scenario1_2+Scenario2_2+Scenario3_2+2.*Scenario1.*(Scenario2+Scenario3)+2.*Scenario2.*Scenario3;
            varsk2 = ESk2- ESk.^2;   %ok
            
            CRN1.Es = ESk ./ Ps;
            
            
            CRN1.Es2   = ESk2./Ps  + (ESk.^2) .* (2*(1-Ps)./(Ps.^2));
            
            CRN1.var_s = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));
            
            if sim.control.debug_active
                %paper notations
                ES_2 =(1 ./ Ps) .* (s1  + (u.^2 ./(2.*(v+u))));  %another notation
                ES2_2 = (1./Ps).*((v./(v+u)) .* s1.^2 + (u./(v+u)) .* (s1.^2 + 2*u.^2/3  + u.* s1)) +  ...
                    (2*(1-Ps)./(Ps.^2)) .* (s1 + (u.^2 ./(2*(u+v)))).^2;
                varS_2 = CRN1.Es2 - CRN1.Es.^2;   %for extr check
                %closed form
                varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
                    (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
                %different notation
                varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                    (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
                    ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));
                
                %coefficient of variation
                CS = sqrt(CRN1.var_s) ./ CRN1.Es;
                CS_2 = sqrt( 1 - 2.*Ps + Ps .*(s1.^2 + (u./(v+u)) .* (2 .* u.^2 /3 + u.*s1))/((s1 + u.^2/(2.*(u+v))).^2));
                
                
                %check various notations
                if mean(abs(CRN1.Es-ES_2)./abs(CRN1.Es+ES_2)) > 0.01
                    warning('Potential error in calculating E(CRN:S) Err:%f<>%f\n', mean(CRN1.Es), mean(ES_2));
                end
                
                if mean(abs(CRN1.Es2-ES2_2)./abs(CRN1.Es2+ES2_2)) > 0.01
                    warning('Potential error in calculating E(CRN:S2) Err:%f<>%f\n', mean(CRN1.Es2), mean(ES2_2));
                end
                
                if mean(abs(CRN1.var_s-varS_2)./abs(CRN1.var_s+varS_2)) > 0.01 || ...
                        mean(abs(CRN1.var_s-varS_3)./abs(CRN1.var_s+var_3)) > 0.01 || ...
                        mean(abs(CRN1.var_s-varS_4)./abs(CRN1.var_s+varS_4)) > 0.01
                    warning('Potential error in calculating var(CRN:S) Err:%f<>%f<>%f<>%f \n', mean(CRN1.var_s), mean(CRN1.var_s_2), mean(CRN1.var_s3), mean(CRN1.var_s4));
                end
                
            end
            
            CRN1.ET = SN.ET; CRN1.ET2 = SN.ET2; CRN1.varT = SN.varT; CRN1.CT2=SN.CT2;  CRN1.Ef=SN.Ef; %remaind unchanged
            CRN1.ru = CRN1.Es ./ CRN1.ET
            CRN1.Cs2 = CRN1.var_s./(CRN1.Es.^2);
             CRN1.Ew = 0.5 * (CRN1.ru ./ (1-CRN1.ru)) .* (CRN1.CT2 + CRN1.Cs2) .* CRN1.Es;  % % % % %     if adjust_kingman, Ew  = Ew - ET2 /(2*ET);  end
            CRN1.Ew = CRN1.ru.*0;
            
            CRN1.Ed =  CRN1.Es +CRN1.Ew+ CRN1.Ef;
            
            %save tp3; disp('tp3 saved in fanalyzesimNB'); pause;
        end
        
    else
        fprintf('Error Invalid CRN approximation mode:%d \n', sim.CRN.Approx);
    end
end

sim.Res.A_SN = SN;
if sim.cognitive && sim.CRN.Approx == 1
    sim.Res.A_CRN1 = CRN1;
end
if sim.cognitive && sim.CRN.Approx == 0
    sim.Res.A_CRN1 = CRN1;
end
if sim.control.debug_active
    if sim.CodedSystem   %Modifications
        tit=fprintf('Coded Heq:%d  Neq:%d  alpha_eq:%1.4f  Rch:%f Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f] ', H,N,a, sim.ch.Rch, sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD);
    else
        tit = fprintf('Uncoded H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f', H,N,a, B, sim.ch.Rch);
    end
end
% pause;
if 0 && sim.control.runtime_plot_active
    if sim.CodedSystem   %Modifications
        tit=sprintf('Coded Heq:%d  Neq:%d  alpha_eq:%1.4f  Rch:%f Code Rate[%f %f] Tail Bits[%d %d]  BER imp[%1.3f %1.3f] ', H,N,a, sim.ch.Rch, sim.coding.CRH, sim.coding.CRD, sim.coding.nTbitsH,sim.coding.nTbitsD,   sim.coding.impBERH,sim.coding.impBERD);
    else
        tit = sprintf('Uncoded H:%d  N:%d  alpha:%1.4f   BER:%f  Rch:%f', H,N,a, B, sim.ch.Rch);
    end
    
    for nmode=0:sim.cognitive
        if nmode==0   %single network
            Res = sim.Res.A_SN;
            tit =[' SN ', tit];
        else    %cognitive network
            Res = sim.Res.A_CRN1;
            tit =[' CRN ', tit];
        end
        
        
        figure; lgnd={'BER:0', ['BER:',num2str(sim.ch.BER)]};
        subplot(231); semilogy(Kv, Res.Es*lambda, 'r-x');  title('E[s]'); legend(lgnd);  %was h=stem(), hb = get(h,'Baseline');set(hb,'Visible','off');set(gca,'yscal','log')
        subplot(232); semilogy(Kv, Res.var_s*lambda, 'r-x');  title('Var[s]');  legend(lgnd);
        subplot(233); semilogy(Kv, Res.var_s ./ (Res.Es.^2), 'r-x');  title('Var[s]/E[s]^2'); legend(lgnd);
        
        subplot(234); stem(Kv, Res.ET*lambda, 'r-x');  title('E[T]'); legend(lgnd);
        subplot(235); stem(Kv, Res.varT*lambda, 'r-x');  title('Var[T]'); legend(lgnd);
        subplot(236); stem(Kv, Res.varT./ (Res.ET.^2), 'r-x');  title('Var[T]/E[T]^2'); legend(lgnd);
        pause;
        
        
        figure;
        semilogy(Kv, Res.Ef*lambda,'b-x'); hold on;
        semilogy(Kv, Res.Es*lambda,'g-x');
        semilogy(Kv, Res.Ew*lambda,'c-x');
        semilogy(Kv, Res.Ed*lambda,'r-x');
        legend({'E[f]', 'E[s]', 'E[w]', 'E[d]'});
        title(['Error Ch.: ', tit]); pause;
        
        title(tit); pause;
        %ylim([0, 100000]);
        % hold on;
    end
end
%             %%%%debug remove after debug%$%%%
%            load('sim1.mat')
%              ESk2 = (p.*s1).^2+((1-p).^2).*(ES_AvArrive+s1+v1./2).^2+2.*p.*(1-p).*(s1).*(ES_AvArrive+s1+v1./2)
%
%                 CRN1.Es = ESk;
%
%                 CRN1.Es2   = ESk2;
%
%                 CRN1.var_s = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));
%
%          sim.Res.S_CRN= sim1.Res.S_CRN;
%          plot_data_fileNB(sim.datapath,sim.fname, 0,'');
%
%             %%%%%
return