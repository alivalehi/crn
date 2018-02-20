if 0
%check c+x distribution where x is expo(a) and compare with expo(a+c)
 
   % close all; clc;

    if 0
    clear;
    mu=1;
    xv = exprnd(mu, [1, 1e8]);
    [fx , x] = hist(xv, 30);

    xv2 = 1+xv;
    [fx2 , x2] = hist(xv2, 30);


    xv3 = exprnd(mu+1, [1, 1e8]);
    [fx3 , x3] = hist(xv3, 30);
    end

    [mean(xv) mean(xv2) mean(xv3)]
    [std(xv) std(xv2) std(xv3)]
    [mean(xv.^2) mean(xv2.^2) mean(xv3.^2)]
    [std(xv.^2) std(xv2.^2) std(xv3.^2)]
    [mean(xv.^3) mean(xv2.^3) mean(xv3.^3)]


    figure;
    % subplot(311);
    plot(x, fx, 'b'); hold on; pause;
    % subplot(312);  
    plot(x2-1, fx2, 'g'); pause;
    % subplot(313);
    plot(x3, fx3, 'r');
end


if 0  %a=mean channel avail time,   b=mean channel not avail time
    %s = s1 (with prob. a/a+b)   s=s1+x (with prob b/a+b) and x=expo(b)
    clear; clc;
    a=1; b=1;  n=1e5;
    s1 = 3;
 
    sv = s1 * ones(1, floor(n * v/(v+u)));
    su = s1 + exprnd(b, [1, floor(n * u/(v+u))]);
    sk = [sv, su];
 
    ESk = (v/(v+u)) * s1 + (u/(v+u)) * (s1 + b);
    ESk2 = (v/(v+u)) * s1^2 + (u/(v+u)) * (s1^2 + 2*b^2 + 2*s1*b);
    varsk = b^3*(2*a+b) / (a+b)^2;
    varsk2 = ESk2-ESk^2;
    if abs(varsk-varsk2)/abs(varsk+varsk2) > 0.01, warning('Potential error in calculating var(CRN:S) Err:%f<>%f \n', varsk1, varsk2); end
 
    fprintf('a:%f b:%f   s1:%f   ES:[%f  %f],  ES2::[%f  %f],  var(s):[%f  %f  %f] \n', a,b,s1, ESk, mean(sk), ESk2, mean(sk.^2), varsk, varsk2, var(sk));
 

    %Approx 1
    %SK = s1+s2+...+sk, where k is geometric(Ps)
    Ps = 1;
    K = 1 + geornd(Ps,[1,n]);
    maxK = max(K);
 
    S1 = s1 * ones(maxK, n);
    S2 = s1 + exprnd(b, [maxK, n]);
    u = rand(maxK, n);
    mask1 = (u <= v/(v+u)); mask2=1-mask1;
    S12 = mask1.* S1 + mask2 .* S2;
    SK=zeros(1,n);
    for i=1:n
       SK(i) =  sum(S12([1:K(i)], [i]));
    end
 
    CRN1.ES = ESk / Ps;
    CRN1.ES2 = ESk2/Ps + ESk^2 *(2*(1-Ps)/Ps^2);
    CRN1.varS = ESk2/Ps + ESk^2 *((1-2*Ps)/Ps^2);
    CRN1.varS_2 = CRN1.ES2 - CRN1.ES^2;
    if abs(CRN1.varS-CRN1.varS_2)/abs(CRN1.varS+CRN1.varS_2) > 0.01, warning('Potential error in calculating var(CRN:S) Err:%f<>%f \n', CRN1.varS, CRN1.varS_2); end
 
    fprintf('a:%f b:%f PS:%f  s1:%f   ESK:[%f  %f],  ESK2::[%f  %f],  var(SK):[%f  %f  %f] \n', a,b,Ps, s1, CRN1.ES, mean(SK), CRN1.ES2, mean(SK.^2), CRN1.varS, CRN1.varS_2, var(SK));

end



if 0  %a=mean channel avail time,   b=mean channel not avail time
    %s = s1 (with prob. a/a+b)   s=s1+x * ux (with prob b/a+b) and x=expo(b), ux=uniform(0,x)
    clear; clc;
    v=1;  %mean ch. available interval time
    u=4;  %mean ch. busy interval time
    n=1e6;
    s1 = 3;
 
    nv=floor(n * v/(v+u)); nu=floor(n * u/(v+u));
    sv = s1 * ones(1, nv);          %avail interval
    su = s1 + exprnd(u, [1, nu]);   %busy interval
 
 
    %split with a uniform splitting
    ux = rand([1, nu]);
    w = exprnd(u, [1, nu]).*ux; suu=s1 + w;
 
 
    Ew = u/2; Ew2 = 2*u^2/3;
    fprintf('average waiting time in busy interval w:  Ew:[%f  %f],  Ew2::[%f  %f]\n', Ew, mean(w),Ew2, mean(w.^2)); pause;
 

    %average delay
    sk = [sv, suu];
 
    ESk = (v./(v+u)) .* s1   + (u./(v+u)) .* (s1 + u/2);  %ok
    ESk2 = (v./(v+u)) .* s1^2 + (u./(v+u)) .* (s1^2 + 2*u.^2/3 + s1.*u);  %ok
    varsk =(u.^3) .* (8*v+5*u)./(12 *((v+u).^2));   %ok
    varsk2 = ESk2-ESk.^2;
    if abs(varsk-varsk2)/abs(varsk+varsk2) > 0.01, warning('Potential error in calculating var(CRN:S) Err:%f<>%f \n', varsk, varsk2); end
 
    fprintf('Si analysis: v:%f u:%f   s1:%f   ES:[%f  %f],  ES2::[%f  %f],  var(s):[%f  %f  %f] \n', v,u,s1, ESk, mean(sk), ESk2, mean(sk.^2), varsk, varsk2, var(sk));
 
    %Approx 1
    %SK = s1+s2+...+sk, where k is geometric(Ps)
    Ps = 0.2;
    K = 1 + geornd(Ps,[1,n]);
    maxK = max(K);
 
    S1 = s1 * ones(maxK, n);
    S2 = s1 + exprnd(u, [maxK, n]).* rand([maxK, n]);
    uu = rand(maxK, n);
    mask1 = (uu <= v/(v+u)); mask2=1-mask1;
    S12 = mask1.* S1 + mask2 .* S2;
    SK=zeros(1,n);
    for i=1:n
       SK(i) =  sum(S12([1:K(i)], [i]));
    end
 

    CRN1.ES = ESk ./ Ps;
ES_2 =(1 ./ Ps) * (s1  + (u.^2 ./(2*(v+u))));  %another notation
 
 
 
    CRN1.ES2   = ESk2./Ps  + (ESk.^2) .* (2*(1-Ps)./(Ps.^2));
ES2_2 = (1./Ps).*((v./(v+u)) .* s1.^2 + (u./(v+u)) .* (s1.^2 + 2*u.^2/3  + u.* s1)) +  ... %another notation
    (2*(1-Ps)./(Ps.^2)) .* (s1 + (u.^2 ./(2*(u+v)))).^2;
 
CRN1.varS = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));
varS_2 = CRN1.ES2 - CRN1.ES.^2;   %for extr check
%closed form
varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
    (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
%different notation
varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                  (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
               ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));
varS_2 = CRN1.ES2 - CRN1.ES.^2;   %for extr check
%closed form
varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
    (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
%different notation
varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                  (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
               ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));

%check var(sum(s)) = E[sum(var)]   =E[r]var ??? [not verified by
%simulations]
varS_5 = varsk ./ Ps;
varS_6 = (u.^3*(8*v+5*u)./(12 *(v+u).^2)) ./ Ps;


   CS = sqrt(CRN1.varS) ./ CRN1.ES;
   CS_2 = sqrt( 1 - 2*Ps + Ps *(s1.^2 + (u./(v+u)) .* (2 * u.^2 /3 + u.*s1))/((s1 + u.^2/(2*(u+v))).^2));
 
 
   %check various notations
   if mean(abs(CRN1.ES-ES_2)./abs(CRN1.ES+ES_2)) > 0.01
      warning('Potential error in calculating E(CRN:S) Err:%f<>%f\n', mean(CRN1.ES), mean(ES_2));
   end
 
   if mean(abs(CRN1.ES2-ES2_2)./abs(CRN1.ES2+ES2_2)) > 0.01
      warning('Potential error in calculating E(CRN:S2) Err:%f<>%f\n', mean(CRN1.ES2), mean(ES2_2));
   end
 
   if mean(abs(CS-CS_2)./abs(CS+CS_2)) > 0.01
      warning('Potential error in calculating C(CRN:S) Err:%f<>%f\n', mean(CS), mean(CS_2));
   end
 
   if mean(abs(CRN1.varS-varS_2)./abs(CRN1.varS+varS_2)) > 0.01 || ...
      mean(abs(CRN1.varS-varS_3)./abs(CRN1.varS+varS_3)) > 0.01 || ...                    
      mean(abs(CRN1.varS-varS_4)./abs(CRN1.varS+varS_4)) > 0.01 
           warning('Potential error in calculating var(CRN:S) Err:%f<>%f<>%f<>%f \n', mean(CRN1.varS), mean(CRN1.varS_2), mean(CRN1.varS3), mean(CRN1.varS4));
   end
 
    fprintf('S=sum[Si](1:r) analysis: v:%f u:%f PS:%f  s1:%f \n  ESK:[%f  %f  s:%f], \n  ESK2:[%f  %f    s:%f], \n  var(SK):[%f  %f   %f  %f  s:%f]  \n  C(SK):[%f  %f   s:%f] \n', ...
        v,u,Ps, s1, CRN1.ES, ES_2, mean(SK), CRN1.ES2, ES2_2, mean(SK.^2), CRN1.varS, varS_2, varS_3, varS_4,  var(SK), CS, CS_2, std(SK)./mean(SK));
end
if 1
    clear all;
    k=1000;
    u=4;v=1;
    p=0.1;
    N = geornd(p,[k,1])+1;
    X = zeros(k,100);
    Y = zeros(k,1);
    Z = zeros(1,k);
   
     if 0 %check exponential 

          for l=1:k
            X  =   exprnd(u, [1,10000]);   %arrive at busy
          end
       histogram(X);
        hold;
        Z = exppdf(1:10000,u);
        plot(Z);

     end
     if 0 %check theorem for exponential R.V
        for l=1:k
            X(l,1:1000)  =   exprnd(u, [1,1000]);   %arrive at busy
        end

        for l=1:k
            S=N(l,1);
            S=1/p;
            if(S>1000)
               fprintf('N exceeds max'); pause;
            end  
            Y(l)  =   sum(X(l,1:S));   %arrive at busy
        end

         Z  =   gamrnd(1/p,u,[1,1000]);   %arrive at busy
        close all;
        figure;
        histogram(Y,'Normalization','pdf');
        hold all
        histogram(Z,'Normalization','pdf');
        legend('Simulation', 'Analytical')
        mean(Y)
        mean(Z)
        var(Y)
        var(Z)
     end
    
    if 1 %check theorem for gamma R.V using Welch–Satterthwaite equation
       kl=1;
        if 1 
        for L=1:kl %first set of sum of exp R.V 
             N = geornd(p,[k,1])+1;
              for l=1:k
                  X(l,1:k)  =   exprnd(u, [1,k]);   %arrive at busy
              end

              for l=1:k
                  S=N(l,1);
                  if(S>k)
                    S = k;
                  end  
                  Y(l,L) = sum(X(l,1:S));   %arrive at busy
              end
          end
         for L=1:kl  %second set of sum of exp R.V
             N = geornd(p,[k,1])+1;
              for l=1:k
                  X2(l,1:k)  =   exprnd(v, [1,k]);   %arrive at busy
              end

              for l=1:k
                  S=N(l,1);
                  if(S>k)
                    S = k;
                  end  
                  Y2(l,L) = sum(X2(l,1:S));   %arrive at busy
              end
          end
        Yt= Y+Y2;
        else 
            load('matlab.mat');
        end
    
        %%% in this block we consider some of two exponential random
        %%% variable with different lambda
%             a=1:k;
%             Z= (u*v/u-v)*(exp(-v*a)- exp(-u*a));
%             mean(Z)
%             var(Z)
        %%
        %%using Welch–Satterthwaite equation
%             ksum=(u+v)^2/(p*(u^2+v^2))
%             fisum=(u^2+v^2)/(u+v)
%             Z  =   gamrnd(ksum,fisum,[1,1000]);
%             mean(Z)
%             var(Z)
        %%
        %%Expected value and variance using charchteristic functions
            K =1/p^2; 
            Egamma=(K)*((u)+(v))
            Egamma2= ((K)*((K)+1))*((u)^2+(v)^2)+(2*((K)^2)*((u)*(v)));
            vargamma=Egamma2-Egamma^2
        %%
%         close all;
%         figure;
%         histogram(Yt,'Normalization','pdf');
%         hold all
%         histogram(Z,'Normalization','pdf');
        mean(Yt)
        var(Yt)
     end
 end
if 0  %a=mean channel avail time,   b=mean channel not avail time
    %s = s1 (with prob. a/a+b)   s=s1+x * ux (with prob b/a+b) and x=expo(b), ux=uniform(0,x)
    clear; clc;
    v=1;  %mean ch. available interval time
    u=4;  %mean ch. busy interval time
    n=1e6;
    s1 = 3;
    p  = exp(-s1/u)
    nv=floor(n *p* v/(v+u));nsv=floor(n *(1-p)* v/(v+u)); nu=floor(n * u/(v+u));
    sv = s1 * ones(1, nv);          %avail interval
    p  = exp(-s1/u)%*(u/(u+v));
    N = geornd(p)+1;
    ssv = zeros(1,nsv);
    su = zeros(1,nu);
   
    for k=1:nu
        su(k)  =   sum(exprnd(u, [1,N])) + sum(exprnd(v, [1, N]));   %arrive at busy
    end
    for k=1:nsv
        ssv(k) =   sum(exprnd(u, [1,N])) + sum(exprnd(v, [1, N+1]));   %arrive at small available
    end
    %split with a uniform splitting
    ux = rand([1, nu]);
    wu = exprnd(u, [1, nu]).*ux; suu=s1 + wu + su;
    vx = rand([1, nsv]);
    wv = exprnd(u, [1, nsv]).*vx; svv=s1 + wv + ssv;
    
 
    % Func =(j*a*b*(-(exp(b*x)*expint(t*(a - j*x))) + exp(a*x)*expint[x*(b - j*x)]))/ ((a - b)*exp((a + b)*x))
    %average delay
    sk = [sv, suu,svv];
    ES_Ep= (1/p)-1;
    %     ES_BusyArrive = ((1./(u.*p))+(1./(v.*(p))));
    %     ES_AvArrive = ((1./(u.*p))+(1./(v.*(p))));
    ES_BusyArrive = (1/p)*((u)+(v));
    ES_AvArrive = ((1/p)*(u))+(((1/p)-1)*(v));
    fprintf('S A :sigmau:%f %f sigmav:%f %f\n',mean(su),ES_BusyArrive,mean(sv),ES_AvArrive);
    Scenario1 = (p* v/(v+u)) .* s1
    Scenario2 =(u./(v+u)).* (s1 + u/2+ ES_BusyArrive);
    Scenario3 = ((v./(v+u)).*(1-p)).* (ES_AvArrive+s1+v/2);
    ESk = Scenario1 +Scenario2+Scenario3;
    lamda1p1 = 1/(p*u);
    lamda2p2 = 1/(p*v);
    ESk2 =(-2/(lamda1p1*lamda2p2))-(2/(lamda1p1)^2)- (2/(lamda2p2)^2);
    varsk2 = ESk2- ESk.^2;   %ok
 
   % if abs(varsk-varsk2)/abs(varsk+varsk2) > 0.01, warning('Potential error in calculating var(CRN:S) Err:%f<>%f \n', varsk, varsk2); end
 
    fprintf('Si analysis: v:%f u:%f   s1:%f   ES:[%f  %f],  ES2::[%f  %f],  var(s):[%f  %f] \n', v,u,s1, ESk, mean(sk), ESk2, mean(sk.^2), varsk2, var(sk)); pause;
 
    %Approx 1
    %SK = s1+s2+...+sk, where k is geometric(Ps)
    Ps = 0.2;
    K = 1 + geornd(Ps,[1,n]);
    maxK = max(K);
 
    S1 = s1 * ones(maxK, n);
    S2 = s1 + exprnd(u, [maxK, n]).* rand([maxK, n]);
    uu = rand(maxK, n);
    mask1 = (uu <= v/(v+u)); mask2=1-mask1;
    S12 = mask1.* S1 + mask2 .* S2;
    SK=zeros(1,n);
    for i=1:n
       SK(i) =  sum(S12([1:K(i)], [i]));
    end
 

    CRN1.ES = ESk ./ Ps;
ES_2 =(1 ./ Ps) * (s1  + (u.^2 ./(2*(v+u))));  %another notation
 
 
 
    CRN1.ES2   = ESk2./Ps  + (ESk.^2) .* (2*(1-Ps)./(Ps.^2));
ES2_2 = (1./Ps).*((v./(v+u)) .* s1.^2 + (u./(v+u)) .* (s1.^2 + 2*u.^2/3  + u.* s1)) +  ... %another notation
    (2*(1-Ps)./(Ps.^2)) .* (s1 + (u.^2 ./(2*(u+v)))).^2;
 
CRN1.varS = ESk2./Ps  + (ESk.^2) .* ((1-2*Ps)./(Ps.^2));
varS_2 = CRN1.ES2 - CRN1.ES.^2;   %for extr check
%closed form
varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
    (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
%different notation
varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                  (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
               ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));
varS_2 = CRN1.ES2 - CRN1.ES.^2;   %for extr check
%closed form
varS_3 = (1./Ps).*(-s1.^2 - s1.*(u.^2)./(u+v) + (2*u.^3) ./(3*(u+v)) -(u.^4)./(2*(u+v).^2) ) + ...
    (1./(Ps.^2)) .* (s1.^2 + s1.*(u.^2) ./(u+v) + (u.^4) ./ (4*(u+v).^2));
%different notation
varS_4 =  ((1-Ps)./(Ps.^2))   .* (s1.^2 + s1.* (u.^2) ./(u+v)) + ...
                  (1./Ps)          .*  ((2*u.^3)./(3*(u+v))) + ...
               ((1-2*Ps)./(Ps.^2)) .* ((u.^4)./(4*(u+v).^2));

%check var(sum(s)) = E[sum(var)]   =E[r]var ??? [not verified by
%simulations]
varS_5 = varsk ./ Ps;
varS_6 = (u.^3*(8*v+5*u)./(12 *(v+u).^2)) ./ Ps;


   CS = sqrt(CRN1.varS) ./ CRN1.ES;
   CS_2 = sqrt( 1 - 2*Ps + Ps *(s1.^2 + (u./(v+u)) .* (2 * u.^2 /3 + u.*s1))/((s1 + u.^2/(2*(u+v))).^2));
 
 
   %check various notations
   if mean(abs(CRN1.ES-ES_2)./abs(CRN1.ES+ES_2)) > 0.01
      warning('Potential error in calculating E(CRN:S) Err:%f<>%f\n', mean(CRN1.ES), mean(ES_2));
   end
 
   if mean(abs(CRN1.ES2-ES2_2)./abs(CRN1.ES2+ES2_2)) > 0.01
      warning('Potential error in calculating E(CRN:S2) Err:%f<>%f\n', mean(CRN1.ES2), mean(ES2_2));
   end
 
   if mean(abs(CS-CS_2)./abs(CS+CS_2)) > 0.01
      warning('Potential error in calculating C(CRN:S) Err:%f<>%f\n', mean(CS), mean(CS_2));
   end
 
   if mean(abs(CRN1.varS-varS_2)./abs(CRN1.varS+varS_2)) > 0.01 || ...
      mean(abs(CRN1.varS-varS_3)./abs(CRN1.varS+varS_3)) > 0.01 || ...                    
      mean(abs(CRN1.varS-varS_4)./abs(CRN1.varS+varS_4)) > 0.01 
           warning('Potential error in calculating var(CRN:S) Err:%f<>%f<>%f<>%f \n', mean(CRN1.varS), mean(CRN1.varS_2), mean(CRN1.varS3), mean(CRN1.varS4));
   end
 
    fprintf('S=sum[Si](1:r) analysis: v:%f u:%f PS:%f  s1:%f \n  ESK:[%f  %f  s:%f], \n  ESK2:[%f  %f    s:%f], \n  var(SK):[%f  %f   %f  %f  s:%f]  \n  C(SK):[%f  %f   s:%f] \n', ...
        v,u,Ps, s1, CRN1.ES, ES_2, mean(SK), CRN1.ES2, ES2_2, mean(SK.^2), CRN1.varS, varS_2, varS_3, varS_4,  var(SK), CS, CS_2, std(SK)./mean(SK));
end