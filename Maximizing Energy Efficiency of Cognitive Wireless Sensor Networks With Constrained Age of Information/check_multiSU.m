clear all
P_s = 10;
P_p = 60;
v=1;
u=2;
SUn=[1 5 10 15:10:55 100 500];
for i=1:length(SUn)
    m = SUn(i);
    trafficmode = 2;
 N = 10;  H = 80;  BER=5e-3; 
Kv = 1:30;
Plenv = Kv * N + H;
    if(trafficmode == 1) %heavy
        Rch = 800;
       lambda = 5;
    else
        Rch = 1000;
        lambda = 1;
    end
s1=Plenv./Rch;
P_t = (lambda*s1)./Kv;
if(m>1)
    P_c = 1- m.*P_t.*(1-P_t).^(m-1);
    P_p = P_p+(P_p*.2*m);
else
    P_c=1;
end
    BER = P_c * BER;
PERv = 1 - (1 - BER) .^ Plenv;




%sim.KvA, sim.KvA, sim.PERv,sim.Plenv,disp('del1'); pause;

H0=H; N0=N;
B0 = BER;
a0 = 1 - B0;

a = a0; B=B0;
Lenv = Kv.*N+H;
PERv = 1-(1 - B).^ Lenv;

SN=[];
SN.Lenv_anal = Lenv;
SN.PERv_anal=PERv;
SN.Ef= 0.5*(Kv-1)./ lambda; %packet formation time for number based policy

%moments of inter-arrival time for Gamma(shape=k, rate=lambda, scale=1/lambda)
SN.ET = Kv / lambda;
SN.ET2 = (Kv .*(Kv+1)) ./ (lambda.^2);
SN.varT=Kv / lambda^2;

%moments of service time :L/Rch * gemoetric(Ps = 1-PER = (1-b)^L= (1-b)^(KN+H)
Ps= 1-PERv;

SN.Es = (s1) ./ Ps;  %Geometric distributed

additive_delay = SN.Es*((m-1)/2);
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
Er = ((1 - BER) .^( -Plenv));
EE = ((Er.*(Kv.*N+H).*(P_s))+(P_p))./Kv.*N;

p              = exp(-s1./v);
v_e             = ((((-v.*s1-1).*exp(-v.*s1))+1)./v)./(1-exp(-s1.*v));
E_uv_2         = 2.*u.^2+((-(s1.^2+2.*s1.*v+2.*v.^2).*exp(-s1./v)+2.*v.^2)./(1-exp(-s1./v)))+2.*u.*v_e;
a              = 1;
ES_BusyArrive  = ((1./p)-a).*((u)+(v_e));
ES_AvArrive    = ((1./p).*(u))+((((1./p)-a)).*(v_e));
E_bu2          = E_uv_2 .*((1./p)-a)+(((2-p)./p.^2)-((1+2.*a)./p)+a.^2+a).*(u+v_e).^2;
ES_BusyArrive2 = 2.*u.^2./3 + E_bu2+s1.^2 +2.*s1.*(u./2+ES_BusyArrive)+ 2.*ES_BusyArrive.*u./2;
ES_AvArrive2   = ES_BusyArrive2+2.*u.^2+2.*u.*ES_BusyArrive;
E_sav2         = E_bu2+(2.*u.^2)+2.*u.*((1./p)-a).*(u+v_e);
var_bu         = E_bu2-ES_BusyArrive.^2;
var_sav        = E_sav2-ES_BusyArrive.^2; % it should be (E_sav2- ES_AvArrive^2);
Scenario1 = (p) .* s1; %arrive at available time slot and it has enough time to do service
Scenario2 =.3*(u./(v+u)).* (s1 + u./2+ ES_BusyArrive);%arrive at busy time slot and it has enough time to do service
Scenario3 = ((1-p)).* (ES_AvArrive+s1+v_e./2);%arrive at available time slot but it doesnt have  enough time to do service

Scenario1_2 = ((p).^2 .* s1.^2);
Scenario2_2 = .3*((u./(v+u)).^2.* (s1.^2 + 2.*u.^2/3 + ES_BusyArrive2 + 2.*s1.*(u/2 +ES_BusyArrive)+2.*u/2.*ES_BusyArrive));
Scenario3_2 = ((1-p)).^2.* (s1.^2 + 2.*v_e.^2/3 + ES_AvArrive2 + 2.*s1.*(v_e/2+ES_AvArrive)+2.*v_e/2.*ES_AvArrive);
if(trafficmode == 1) %heavy
    ESk2 =Scenario1_2+Scenario3_2+2.*Scenario1.*(Scenario3);
    ESk = Scenario1+Scenario3;
else
    ESk = Scenario1+Scenario2+Scenario3;
    ESk2 =Scenario1_2+Scenario2_2+Scenario3_2+2.*Scenario1.*(Scenario2+Scenario3)+2.*Scenario2.*Scenario3;
end

varsk2 = ESk2- ESk.^2;   %ok
Es = ESk ./ Ps;
Es2   = ESk2./Ps  + (ESk.^2) .* (2*(1-Ps)./(Ps.^2));
var_s = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));

ET = SN.ET; CRN1.ET2 = SN.ET2; CRN1.varT = SN.varT; CT2=SN.CT2;  Ef=SN.Ef; %remaind unchanged
ru = Es ./ ET ;
Cs2 = var_s./(Es.^2);
Ew = 0.5 * (ru ./ (1-ru)) .* (CT2 + Cs2) .*  Es;  % % % % %     if 1, Ew  = Ew - ET2 /(2*ET);  end
Ew(find( ru >= 1)) = Ew(find( ru <1,1,'first'))+(10./Kv(find( ru >= 1))).*Ew(find( ru <1,1,'first'));
Ed =  Ew +  Es +  Ef+additive_delay;
plot(Ed);
% hold on
figure;
plot(EE);
[Ed_A,Kd_A] = min(Ed);
[EE,Ke] = min(EE);
fprintf ('m: %f Kd:%f Ed:%f Ke:%f EE: %f\n',m,Kd_A,Ed_A,Ke,EE);
end