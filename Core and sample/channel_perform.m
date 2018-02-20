function [sim] = channel_perform(sim,nch)
 v = sim.CRN.Ch_AvailMeanTime(1); u=sim.CRN.Ch_NotAvailMeanTime(1);
BERch =[];
chf = [0];
sim.utrack=[];
sim.vtrack=[];
sim.slottrack  = [];

avt  = exprnd(v , [1, nch]);     %duration with available channel
avnt = exprnd(u, [1, nch]);  %duration with no available channel
avnvt = reshape([avt; avnt], [1, 2*nch]);
chf = reshape([ones(1,nch); zeros(1,nch)], [1, 2*nch]);
cht = zeros(size(avnvt));
for i = 1: length(cht)-1
    cht(i+1)=cht(i)+avnvt(i);
end
display('Channel perform done!\n');
sim.chf = chf;
sim.cht = cht;
sim.chs = avnvt;

