if 1
    
    %check that summation of n independent exponential variables are gamma distribution with scale parameter and shape parameter 
    %solution 1-Induction
    %solution 2-MGF


    %conflicting refernces:
    %link 1: http://math.stackexchange.com/questions/300816/exponential-random-variable-question?rq=1
    %[sum (X1,X2,...Xn) is also an exponetial with mean n/lambda !!!

    %link 2: http://math.stackexchange.com/questions/655302/gamma-distribution-out-of-sum-of-exponential-random-variables
    %[sum (X1,X2,...Xn) is Gamma(shape:n, rate:lambda or scale:1/lambda) !!!


    clear; close all; clc; N=3;
    lambda=2; beta=1/lambda; mu=beta; %parameters mu=mean[x]=1/lambda rate:lambda
    hx = 0.1:0.3:10; dx = hx(2)-hx(1);
    x = exprnd(mu,[N,1e6]);
    figure;
    for n=1: N
        y = sum(x(1:n,:),1);
        [pY] = hist(y, hx); pY=pY/(dx*sum(pY));
        subplot(N,1,n); %subplot(N,2,1+2*(n-1)); 
        plot(hx, pY, 'b-*'); title(['Simulation n=', num2str(n)]); hold on; 


        %subplot(N,2,2*n); 
        if n==1    %expo pdf
            pY2 = exppdf(hx,n*mu); 
            plot(hx, pY2, 'r-o'); title(['Analytic n=', num2str(n)]); hold on; 
        else      %gamma pdf
            shape=n; scale=1/lambda; rate=lambda;
            pY2 = gampdf(hx,shape, scale); 
            plot(hx, pY2, 'r-o'); title(['Analytic n=', num2str(n)]); hold on; 
        end
    end
end