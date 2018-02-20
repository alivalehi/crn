clear all
v=3;
u=1;
n=1e6;
s1= 2;
% for j =1: length(u)
% x=exprnd(v,[1,n]);
% y=exprnd(u(j),[1,n]);
% a=rand(1,n);
% w = a.*x;
% %w = w(find(w>s));
% %result = zeros(1,n);
% %result = x-w;
% for i=1:n
%     
%     q(i)=length(find(w>(s1+y(i))));
% end
% result(j) = mean(q./n);
% end
% %plot(result);
% log(mean(result))
% fprintf('%f %f',mean(result),v/2);
x=exprnd(1/v,[1,n]);
X = x(find(x<(s1)));
fprintf('mean X:%f analytical 1:%f analytical 2:%f\n',mean(X),(((((-v*s1-1)*exp(-v*s1))+1)/v)/(1-exp(-s1*v)))^2,v)
fprintf('mean2 X:%f analytical 1:%f analytical 2:%f\n',mean(X.^2),(((-(s1^2)-(2*s1/v)-(2/(v^2)))*exp(-v*s1))+2/(v^2))/(1-exp(-v*s1)),((((-v*s1-1)*exp(-v*s1))+1)/v)/(1-exp(-s1*v))^2)
fprintf('mean X:%f analytical 1:%f analytical 2:%f',var(X),(v^2)*exp(-2*s1*v)+exp(-s1*v)*(-s1-6*v^2-4*v*s1)+3*(v^2),v)
(v^2)*exp(-2*s1*v)+exp(-s1*v)*(-s1-6*v^2-4*v*s1)+3*(v^2)