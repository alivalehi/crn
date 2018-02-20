function [output_args] = splot(x,y,attr)
% step  plot: time points[x] values[y]
if ~issorted(x), warning('x is not sorted, possibly graph misrepresentation!!!\n'); end

xprev = 0;
yprev = 0;
for i = 1: length(x)
    xh = [xprev x(i)];
    yh = [yprev yprev];
    plot(xh,yh, attr);
    hold on;
    
    xv = [x(i) x(i)];
    yv = [yprev y(i)];
    plot(xv,yv, attr);
    
    xprev=x(i);
    yprev=y(i);
end
