function [u_rnd,v_rnd,BER_rnd,sim] = expert_update(u_rnd,v_rnd,BER_rnd,Best_prediction,sim,EE_diff_best,diff)
sim.alpha = (.5+(rand/2));%chi;% coeffient of velocity
sim.beta = sim.chi*rand*sim.phi1; %coeffient of Memory Best %method 1 has 0.1 coefficient for all these three
sim.gamma = sim.chi*rand*sim.phi2; %coeffient of current best

temp_v= v_rnd;
temp_u = u_rnd;
temp_BER = BER_rnd ;
if(isempty(sim.v_past))
    sim.v_past = temp_v;
    sim.u_past = temp_u;
    sim.BER_past = temp_BER ;
    sim.Best_EE_diff = ones(1,length(u_rnd)).*Inf;
    sim.u_Mbest = u_rnd;
    sim.v_Mbest = v_rnd;
    sim.BER_Mbest = BER_rnd;
    sim.u_veloc = zeros(size(u_rnd));
    sim.v_veloc = zeros(size(v_rnd)) ;
    sim.BER_veloc = zeros(size(BER_rnd));
    
end
% Global best
sim.u_best = u_rnd(Best_prediction);
sim.v_best = v_rnd(Best_prediction);
sim.BER_best = BER_rnd(Best_prediction);
%memory Best

Best_EE_diff = sim.Best_EE_diff;
u_Mbest = sim.u_Mbest;
v_Mbest = sim.v_Mbest;
BER_Mbest = sim.BER_Mbest;
for ll=1:length(u_rnd)
    if(diff(ll)<Best_EE_diff(ll))
        Best_EE_diff(ll) = diff(ll);
        u_Mbest(ll) = u_rnd(ll);
        v_Mbest(ll) = v_rnd(ll);
        BER_Mbest(ll) = BER_rnd(ll);
    else
        ali=1;
    end
end
sim.u_Mbest = u_Mbest;
sim.v_Mbest = v_Mbest;
sim.BER_Mbest = BER_Mbest;
sim.Best_EE_diff = Best_EE_diff;
%sim.alpha = .99*sim.alpha;
uBEST = u_rnd(Best_prediction);
vBEST = v_rnd(Best_prediction);
BBest = BER_rnd(Best_prediction);
sim.u_veloc = sim.alpha.*sim.u_veloc+sim.beta.*((ones(size(u_rnd)).*sim.u_Mbest)-u_rnd)+sim.gamma.*(u_rnd(Best_prediction)-u_rnd);
sim.v_veloc = sim.alpha.*sim.v_veloc+sim.beta.*((ones(size(v_rnd)).*sim.v_Mbest)-v_rnd)+sim.gamma.*(v_rnd(Best_prediction)-v_rnd);
sim.BER_veloc = sim.alpha.*sim.BER_veloc+sim.beta.*((ones(size(BER_rnd)).*sim.BER_Mbest)-BER_rnd)+sim.gamma.*(BER_rnd(Best_prediction)-BER_rnd);

u_rnd = u_rnd+sim.u_veloc;
v_rnd = v_rnd+sim.v_veloc;
BER_rnd = BER_rnd+sim.BER_veloc;
u_rnd(Best_prediction) = uBEST;
v_rnd(Best_prediction) = vBEST;
BER_rnd(Best_prediction) = BBest;
if((length(find(BER_rnd<0))>0)||(length(find(v_rnd<0))>0)||(length(find(u_rnd<0))>0))
    ali=1;
    u_rnd(find(u_rnd<sim.umin))= sim.umin;
    v_rnd(find(v_rnd<sim.vmin))=sim.vmin;
    u_rnd(find(u_rnd>sim.umax))= sim.umax;
    v_rnd(find(v_rnd>sim.vmax))=sim.vmax;
    %    BER_rnd(find(BER_rnd<0))=temp_BER(find(BER_rnd<0));
end

% v_rnd = v_rnd+(1-sim.alpha).*(v_rnd(Best_prediction)-v_rnd);
% BER_rnd = BER_rnd+(1-sim.alpha).*(BER_rnd(Best_prediction)-BER_rnd);

sim.v_past = temp_v;
sim.u_past = temp_u;
sim.BER_past = temp_BER ;
%     sim.alpha = .5+(rand/2);%chi;% coeffient of velocity
%     sim.beta =  sim.chi*rand*sim.phi1; %coeffient of Memory Best
%     sim.gamma = sim.chi*rand*sim.phi2; %coeffient of current best
% temp_v= v_rnd;
% temp_u = u_rnd;
% temp_BER = BER_rnd ;
% if(isempty(sim.v_past))
% sim.v_past = temp_v;
% sim.u_past = temp_u;
% sim.BER_past = temp_BER ;
% sim.Best_EE_diff = Inf;
% end
%
% sim.u_best = u_rnd(Best_prediction);
% sim.v_best = v_rnd(Best_prediction);
% sim.BER_best = BER_rnd(Best_prediction);
% if(EE_diff_best<sim.Best_EE_diff)
% sim.Best_EE_diff = EE_diff_best;
% [sim.u_Mbest] = u_rnd;
% [sim.v_Mbest] = v_rnd;
% [sim.BER_Mbest] = BER_rnd;
% end
%
% u_rnd = u_rnd+sim.alpha.*(u_rnd - sim.u_past)+sim.beta.*((ones(size(u_rnd)).*sim.u_Mbest)-u_rnd)+sim.gamma.*(u_rnd(Best_prediction)-u_rnd);
% v_rnd = v_rnd+sim.alpha.*(v_rnd - sim.v_past)+sim.beta.*((ones(size(v_rnd)).*sim.v_Mbest)-v_rnd)+sim.gamma.*(v_rnd(Best_prediction)-v_rnd);
% BER_rnd = BER_rnd+sim.alpha.*(BER_rnd - sim.BER_past)+sim.beta.*((ones(size(BER_rnd)).*sim.BER_Mbest)-BER_rnd)+sim.gamma.*(BER_rnd(Best_prediction)-BER_rnd);
% u_rnd(u_rnd<=0)=1;
% v_rnd(v_rnd<=0)=1;
% % v_rnd = v_rnd+(1-sim.alpha).*(v_rnd(Best_prediction)-v_rnd);
% % BER_rnd = BER_rnd+(1-sim.alpha).*(BER_rnd(Best_prediction)-BER_rnd);
%
% sim.v_past = temp_v;
% sim.u_past = temp_u;
% sim.BER_past = temp_BER ;
end