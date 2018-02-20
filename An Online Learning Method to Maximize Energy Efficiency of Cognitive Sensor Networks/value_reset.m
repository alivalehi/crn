function [ sim ] = value_reset( sim )
%VALUE_RESET Summary of this function goes here
%   Detailed explanation goes here
    sim.v_past = [];
    sim.u_past = [];
    sim.BER_past = [] ;
    sim.Best_EE_diff = [];
    sim.u_Mbest = [];
    sim.v_Mbest = [];
    sim.BER_Mbest = [];
    sim.u_veloc = [];
    sim.v_veloc = [];
    sim.BER_veloc = [];

end

