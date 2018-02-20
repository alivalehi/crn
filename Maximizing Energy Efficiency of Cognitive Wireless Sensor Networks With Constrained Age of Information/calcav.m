function [yav] = calcav(t1,y, dt)
% calculate time average of process y
% it inserts a pulse function with at each epoch
% whose length is y(i) and duration is interarrival time t1(i+1)-t1(i)
% t1 should be sorted, smaller dt gives more accurate average

t = t1(1) : dt :t1(end);
if (length(t) / length(t1)) < 10, warning('Possibly small dt !!!'); end
for i= 1:length(t)
    for j = 1:length(t1)-1
        if(t(i) >= t1(j)) &&(t(i) < t1(j+1))
            yt(i) = y(j);
        elseif  t(i) == t1(j+1)
            yt(i) = y(j+1);
        end
    end
end
yav = sum(yt)/length(yt);
