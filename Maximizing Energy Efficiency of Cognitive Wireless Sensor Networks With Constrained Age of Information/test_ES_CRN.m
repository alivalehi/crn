
N=1e6; u=5; v=9;  a=1; s1=3; ps =  exp(-s1/v); %ps = probability of sucess[available time being greater than s1 p[TV>s1]


fprintf('This program verifies E[S] for differet situations in a cognitiveradio\n');
fprintf('Params:  TU~expo(u): Busy Time Slot Length\n         TV~expo(v):Available Time Slot Length\n         s1: Single Packet transmissiontime\n         SV: Service time for Available Sate\n         SW: Servicetime for Available State (in the beginning of TS)\n         SU: Servicetime for Busy Sate\n');
fprintf('Parameters are set [v:%1.2f] [u:%1.2f] [s1:%1.2f]Pr(TV>s1)[ps:%1.2f]\n', u,v,s1, ps);
pause;


%%
k = 1 + geornd(ps, [1, N]);
fprintf('E[k~geo(ps)]: %1.3f   %1.3f   \n', mean(k), 1/ps);
fprintf('E[k2~geo(ps)]: %1.3f   %1.3f  \n', mean(k.^2), (2-ps)/(ps^2));
fprintf('var[k~geo(ps)]: %1.3f   %1.3f \n\n', var(k), (1-ps)/(ps^2));


TU = exprnd(u, [max(k), N]); TV = exprnd(v, [max(k), N]);
TU1 = TU(:); TV1=TV(:);
fprintf('E[x~expo(u)]: %1.3f   %1.3f   \n', mean(TU1), u);
fprintf('E[x2~expo(u)]: %1.3f   %1.3f   \n', mean(TU1.^2), 2*u^2);
fprintf('Var[x2~expo(u)]: %1.3f   %1.3f   \n\n', var(TU1), u^2);

fprintf('E[x~expo(v)]: %1.3f   %1.3f   \n', mean(TV1), v);
fprintf('E[x2~expo(v)]: %1.3f   %1.3f   \n', mean(TV1.^2), 2*v^2);
fprintf('Var[x2~expo(v)]: %1.3f   %1.3f   \n\n', var(TV1), v^2);

T = TU + TV;
%
% fprintf('E[T]: %1.3f  = %1.3f   \n', mean(T), u+v);
% fprintf('E[T]: %1.3f  = %1.3f    \n', mean(T.^2), 2*(u^2+v^2+u*v));
% fprintf('var[T]: %1.3f  = %1.3f    \n\n', var(T),2*(u^2+v^2+u*v)-(u+v)^2);pause;

k=k-a;
yu=zeros(1,N);
for i=1:N
    yu(i) = sum(TU([1:k(i)], i));
end


yv=zeros(1,N);
for i=1:N
    yv(i) = sum(TV([1:k(i)], i));
end






fprintf('Verify moments of X=sum(Xi)[i=1:...:k-a], where Xi is expo(u=%f),a=%d and k is Geo(Psucess=%1.2f)\n', u, a, ps); pause;

fprintf('E[yu]: %1.3f  = %1.3f   \n', mean(yu), u*(1/ps-a));
fprintf('E[yu2]: %1.3f  = %1.3f    \n', mean(yu.^2), u^2*(2/ps^2 - 2*a/ps +a^2-a));
fprintf('var[yu2]: %1.3f  = %1.3f    \n\n', var(yu), u^2*(1/ps^2 - a));

fprintf('E[yv]: %1.3f  = %1.3f   \n', mean(yv), v*(1/ps-a));
fprintf('E[yv2]: %1.3f  = %1.3f    \n', mean(yv.^2), v^2*(2/ps^2 - 2*a/ps +a^2-a));
fprintf('var[yv2]: %1.3f  = %1.3f    \n\n', var(yv), v^2*(1/ps^2 - a));




yt = yu + yv;

fprintf('Verify moments of yt = Xi + Xj [i,j=1:...:k-a], where Xi is sum of k number of expo(u=%f) and Xj is  sum of k number of expo(v=%f), a=%d and kis Geo(Psucess=%1.2f)\n', u,v, a, ps); pause;

fprintf('E[yt]: %1.3f  = %1.3f   \n', mean(yt), v*(1/ps-a)+u*(1/ps-a));
fprintf('E[yt2]: %1.3f  = %1.3f    \n', mean(yt.^2),2*(u^2+v^2+u*v)*((1/ps)-a)+(((2-ps)/ps^2)-((1+2*a)/ps)+a^2+a)*(u+v)^2);
fprintf('var[yt2]: %1.3f  = %1.3f    \n\n', var(yt),(2*(u^2+v^2+u*v)*((1/ps)-a)+(((2-ps)/ps^2)-((1+2*a)/ps)+a^2+a)*(u+v)^2)-(v*(1/ps-a)+u*(1/ps-a))^2);

[k(1:10); y(1:10); x([1:min(size(x,1),5)],[1:10])]


%
fprintf('w=uniform[0, t], t~expo(u=%f)\n', u);
alpha = rand(1,N);
w=alpha.*x1;
fprintf('E[w]: %1.3f   %1.3f   \n', mean(w), u/2);
fprintf('E[w2]: %1.3f   %1.3f   \n', mean(w.^2), 2*u^2/3);
fprintf('Var[w]: %1.3f   %1.3f   \n\n', var(w), 5*u^2/12);



%