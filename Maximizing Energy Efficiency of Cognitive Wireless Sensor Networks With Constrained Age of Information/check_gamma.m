if 1
    
    %check that combining two set of n independent gamma variables are gamma distribution with scale parameter and shape parameter
    
    
    
    clear; close all; clc; N=1;
    lambda=3; beta=1/lambda; mu=beta; %parameters mu=mean[x]=1/lambda rate:lambda
    hx = 0:50; dx = hx(2)-hx(1);
    k1 = 80;
    k2 = 10;
    x1 = gamrnd(k1,mu,[N,1e6]);
    x2 = gamrnd(k1,mu,[N,1e6]);
    x3 = gamrnd(k1,mu,[N,1e6]);
    x4 = gamrnd(k1,mu,[N,1e6]);
    x5 = gamrnd(k1,mu,[N,1e6]);
    x6 = gamrnd(k1,mu,[N,1e6]);
    x7 = gamrnd(k1,mu,[N,1e6]);
    x8 = gamrnd(k1,mu,[N,1e6]);
    xx = [x1,x2,x4,x5];
    cht1 = nan(size(x1));cht2 = nan(size(x2));cht3 = nan(size(x3));cht4 = nan(size(x4));
    cht1(1) = 0;cht2(1) = 0;cht3(1) = 0;cht4(1) = 0;cht5(1) = 0;cht6(1) = 0;cht7(1) = 0;cht8(1) = 0;
    for i = 1: length(cht1)
        cht1(i+1)=cht1(i)+x1(i);
        cht2(i+1)=cht2(i)+x2(i);
        cht3(i+1)=cht3(i)+x3(i);
        cht4(i+1)=cht4(i)+x4(i);
        cht5(i+1)=cht5(i)+x5(i);
        cht6(i+1)=cht6(i)+x6(i);
        cht7(i+1)=cht7(i)+x7(i);
        cht8(i+1)=cht8(i)+x8(i);
    
    end
    cht = [cht1,cht2];%x=,cht3];%,cht4];%,cht5];%,cht6];%,cht7,cht8];
    %cht = [cht1,cht2,cht3,cht4,cht5,cht6,cht7,cht8];
    
    cht = sort(cht);
    for p = 1: length(cht)-1
        inter_pack_time(p) = cht(p+1)-cht(p);
    end
    [pY1]=hist(x1);
    
    [pY2]=hist(x2);
    a=mean(x1)
    b=mean(x2)/8
    c=mean(inter_pack_time)
    d= var(x1)
     e = var(inter_pack_time)
    f=(mean(inter_pack_time.^2)/var(inter_pack_time))-1
    g = var(inter_pack_time)/mean(inter_pack_time)
    f*g
    [pY3]=hist(inter_pack_time);
    
    [pY4]=hist(x3);
    %     plot(pY1);
    %     hold on;
    %     plot(pY2);
 %   plot(pY3);
%     hold on
%     plot(pY4)
%     
end