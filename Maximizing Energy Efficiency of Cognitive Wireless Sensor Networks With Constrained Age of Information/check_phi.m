clear all
v=34;
%vv = [1:1000];
u=1;
n=1e6;
%s1=4
%s=1000;
s= [1:10:600, 600:1000];
for n = 1:length(s)
%for n = 1:length(vv)

       s1 = s(n);
  % v= vv(n);
phiprime=rand([1,n]).*exprnd(v,[1,n]);
phi = phiprime(find(phiprime<(s1)));
alpha(n) = ((1-exp(-s1/v)+(s1/v)*expint(s1/v)));
phianalytic(n) = (v/2-(1/2)*(s1+v)*exp(-s1/v)+((s1^2)/(2*v))*expint(s1/v))/alpha(n);
phianalytic2(n) = (2*(v.^2)/3-(1/3)*((s1^2)+(2*v.^2)+(2*v*s1))*exp(-s1/v)+((s1^3)/(3*v))*expint(s1/v))/alpha(n);
end
plot(phianalytic);
hold on
plot(ones(1,length(phianalytic)).*(v/2)+.001)
figure
plot(phianalytic2);
hold on
plot(ones(1,length(phianalytic2)).*(2*v^2/3)+.001)
set(gca,'fontsize',22)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
fprintf('mean phi'':%f analytical 1:%f \n',mean(phiprime),v/2)
fprintf('mean2 phi'':%f analytical 1:%f\n',mean(phiprime.^2),(2*v^2)/3)
fprintf('mean phi:%f analytical 1:%f \n',mean(phi),phianalytic)
