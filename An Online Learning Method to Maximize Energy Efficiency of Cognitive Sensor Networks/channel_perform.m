function [sim] = channel_perform(sim)
u=sim.u ;
v=   sim.v;
BER =  sim.ch.BER;
nchlist = sim.nch;
kanalarray = [4 6 7 8 9 8 7 6 4];%[11 10 10 9 9 8 8 7 7 6 5 4 4 4 5 6 7 7 8 8 9 9 10 10 11 ];
%for Rch=10 N16 H80[23 23 22 21 20 19 18 17 16 14 13 11 9 7 4 7 9 11 13 14 16 17 18 19 20 21 22 23 23];

BERch =[];
chf = [0];
chs=[];
cht=[];
npackets=1;
rv =[];
sv1= [];
sim.utrack=[];
sim.vtrack=[];
kanal =[];
sim.slottrack  = [];
for l =1: sim.m
    %             nch = 100* ceil(max(tp) / (a+b));
    %             nch = min(1e8, max(1e5, nch));
    %generate channel availability durations
  if(length(sim.nch)>1)
   nch = nchlist(l)*1e4;
  else
      nch = nchlist*1e7;
  end 
      a=v(l);
    b=u(l);
    if(length(BER)>1)
    c = BER(l);
    else
    c = BER;    
    end
    
    avt  = exprnd(a , [1, nch]);     %duration with available channel
    avnt = exprnd(b, [1, nch]);  %duration with no available channel
    vtracktemp = ones(1, nch).*a; 
    utracktemp = ones(1, nch).*b;
    avnvt = reshape([avt; avnt], [1, 2*nch]);
    BERch_t = ones(1,2*nch)*c;
    %channel availability flag
    chf_t = reshape([ones(1,nch); zeros(1,nch)], [1, 2*nch]);
    
    %channel availability process
    cht_t = nan(size(avnvt));
    if(isempty(cht))
        cht_t(1) = 0;
    else
        cht_t(1) = cht(end);
    end
    for i = 1: length(cht_t)-1
        cht_t(i+1)=cht_t(i)+avnvt(i);
    end
    kanal = [kanal, ones(1, 2*nch).*kanalarray(l)];
    chf = [chf, chf_t];
    chcont = [1:length(chf)];
    chs = [chs,avnvt];
    cht = [cht,cht_t];
    BERch = [BERch,BERch_t];
    sim.utrack=[sim.utrack,utracktemp];
    sim.vtrack=[sim.vtrack,vtracktemp ];
    fprintf('Channel perform set %d done!\n',l);
end
sim.slotuv = reshape([sim.vtrack; sim.utrack], [1, length(cht)]);

sim.kanal = kanal;
chf(end)= [];
sim.chf = chf;
sim.cht = cht;
sim.BERch = BERch;
%chcont = [1:length(chf)];
%chs = [chs,avnvt];
%sv_t = (Plen .* rv)/sim.ch.Rch;           %service time
%     rv = [rv,ones(1,sim.m)*];

