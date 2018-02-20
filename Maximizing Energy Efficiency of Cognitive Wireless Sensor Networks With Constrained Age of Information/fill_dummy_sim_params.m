function [sim] = fill_dummy_sim_params(sim)
%fill with dummy data
%sim.T_TB=[]; %for backward compatibility

z=zeros(1,length(sim.Tv));
sim.Rbv= z;
sim.Plenv= z;
sim.PERv= z;

sim.Rbv2= z;
sim.PERv2= z;
sim.Plenv2= z;

sim.ET_approx= z;
sim.ru_approx= z;
sim.Es_approx= z;
sim.Es2_approx= z;
sim.Ew_approx= z;
sim.Ed_approx= z;
sim.ET= z;
sim.ru= z;
sim.Es= z;
sim.Es2= z;
sim.Ew= z;
sim.Ed= z;
sim.ET_0= z;
sim.ru_0= z;
sim.Es_0= z;
sim.Es2_0= z;
sim.Ew_0= z;
sim.Ed_0= z;
end

