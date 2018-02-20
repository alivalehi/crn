function [EE_prediction_vector,P_s,perr] = EE_Analisor(k,Best_prediction,sim,u_rnd,v_rnd,BER_rnd,kflag)
kv= k;
if kflag
%     BER_rnd= BER_rnd(Best_prediction);
%     u_rnd = u_rnd(Best_prediction);
%     v_rnd = v_rnd(Best_prediction);
 else
     kv = kv(1);
end
for j=1: length(kv)
    k= kv(j);
    Plen = k * sim.N + sim.H;
    s = ones(1,length(u_rnd)).*Plen/sim.ch.Rch;
    p1 = (1-((1 - BER_rnd) .^( Plen)));%retransmisson
    ps = exp(-s./(v_rnd));%success
    a11 = p1.*ps;
    a12 = 0;
    a13 = p1.*(1-ps);
    a21 = ps;
    a22 = 0;
    a23 = 1-ps;
    a31 = ps;
    a32 =0;
    a33 = 1-ps;
    vv = v_rnd./(v_rnd+u_rnd);
    uu = u_rnd./(v_rnd+u_rnd);
    pp1 = [vv.*(ps)];
    pp2 = [vv.*(1-ps)];
    pp3 = [uu];
    for l = 1:(length(ps))
       
        Q(:,:,l) =  [a11(l) a12 a13(l);a21(l) a22 a23(l);a31(l) a32 a33(l)];
        N(:,:,l) = inv(eye(length(Q(:,:,l)))-Q(:,:,l));
        x21 = pp1.*N(1,2,l);
        x22(l) = pp2(l).*N(2,2,l);
        x23(l) = pp3(l).*N(3,2,l);
        x11(l) = pp1(l).*N(1,1,l);
        x12(l) = pp2(l).*N(2,1,l);
        x13(l) = pp3(l).*N(3,1,l);
        x31(l) = pp1(l).*N(1,3,l);
        x32(l) = pp2(l).*N(2,3,l);
        x33(l) = pp3(l).*N(3,3,l);
        %%%test
%         x21 = pp1.*N(2,1,l);
%         x22(l) = pp2(l).*N(2,2,l);
%         x23(l) = pp3(l).*N(2,3,l);
        %P_s(l) = xx3(l)+xx2(l);
    end  
    
    P_s1 = x13+x23+x33;%+(p1.*(1-ps));


%%%%%temperory
%kv = kv(1);
%###############
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for j=1: length(kv)
                k= kv(j);
                Plen = k * sim.N + sim.H;
                s = ones(1,length(u_rnd)).*Plen/sim.ch.Rch;
                p1 = (1-((1 - BER_rnd) .^( Plen)));%retransmisson
                ps = exp(-s./(v_rnd));%success
                a11 = 0;
                a12 = 0;
                a13 = 0;
                a21 = ps;
                a22 = 1-ps;
                a23 = 0;
                a31 = ps;
                a32 =(1-ps);
                a33 = 0;
                vv = v_rnd./(v_rnd+u_rnd);
                uu = u_rnd./(v_rnd+u_rnd);
                pp1 = [vv.*(ps)];
                pp2 = [vv.*(1-ps)];
                pp3 = [uu];
                for l = 1:(length(ps))

                    M(:,:,l) =  [a11 a12 a13;a21(l) a22(l) a23;a31(l) a32(l) a33];
                    N(:,:,l) = inv(eye(length(M(:,:,l)))-M(:,:,l));
                    x21 = pp1.*N(1,2,l);
                    x22(l) = pp2(l).*N(2,2,l);
                    x23(l) = pp3(l).*N(3,2,l);
                    x11(l) = pp1(l).*N(1,1,l);
                    x12(l) = pp2(l).*N(2,1,l);
                    x13(l) = pp3(l).*N(3,1,l);
                    x31(l) = pp1(l).*N(1,3,l);
                    x32(l) = pp2(l).*N(2,3,l);
                    x33(l) = pp3(l).*N(3,3,l);
                    %%%test
            %         x21 = pp1.*N(2,1,l);
            %         x22(l) = pp2(l).*N(2,2,l);
            %         x23(l) = pp3(l).*N(2,3,l);
                    %P_s(l) = xx3(l)+xx2(l);
                end  

                P_s = x23+x22+x21;%+(p1.*(1-ps));
                
                
                
               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%       k= kv(j);
%     Plen = k * sim.N + sim.H;
%     s = ones(1,length(u_rnd)).*Plen/sim.ch.Rch;
%     p1 = (1-((1 - BER_rnd) .^( Plen)));%retransmisson
%     ps = exp(-s./(v_rnd));%success
%     a11 = 1-ps;
%     a12 = 0;
%     a13 = 0;
%     a21 = 1-ps;
%     a22 = 0;
%     a23 = 0;
%     a31 = ps;
%     a32 =(1-ps);
%     a33 = 0;
%     vv = v_rnd./(v_rnd+u_rnd);
%     uu = u_rnd./(v_rnd+u_rnd);
%     pp1 = [vv.*(ps)];
%     pp2 = [vv.*(1-ps)];
%     pp3 = [uu];
%     for l = 1:(length(a11))
%        
%         M(:,:,l) =  [a11(l) a12;a21(l) a22];
%         N(:,:,l) = inv(eye(length(M(:,:,l)))-M(:,:,l));
%         x21 = pp1.*N(1,2,l);
%         x22(l) = pp2(l).*N(2,2,l);
%         x23(l) = pp3(l).*N(3,2,l);
%         x11(l) = pp1(l).*N(1,1,l);
%         x12(l) = pp2(l).*N(2,1,l);
%         x13(l) = pp3(l).*N(3,1,l);
%         x31(l) = pp1(l).*N(1,3,l);
%         x32(l) = pp2(l).*N(2,3,l);
%         x33(l) = pp3(l).*N(3,3,l);
%         %%test
%        x21 = pp1.*N(2,1,l);
%        x22(l) = pp2(l).*N(2,2,l);
%        x23(l) = pp3(l).*N(2,3,l);
%         P_s(l) = xx3(l)+xx2(l);
%     end  
%     
%     P_s = x22+x21+(p1.*(1-ps));
  %  P_s2= x11+x12+x13;    
   % P_s3= x11+x12+x13;
    if ~kflag
        P_s;
    %    P_s2
    %    P_s3
   % fprintf('sAV:%f Av:%f B:%f\n',P_s(Best_prediction),P_s2(Best_prediction),P_s3(Best_prediction));
    end
    %  P_s = (pp1.*N(1,2)+pp2.*N(2,2)+pp3.*N(3,2));
    %pp1 = [vv*(1-ps) vv*(ps) uu];
    % pp2 = [vv*(ps) vv*(1-ps) uu];
    % pp1*N;
    % aaaa = pp2*N;
    % P_s  = aaaa(2)%(1./((v_rnd./(v_rnd+u_rnd)).*(1-exp(-s./(v_rnd)))))%number of arrive at small av
    %P_s = (exp(s./(v_rnd)));
    
    perr = ((1 - BER_rnd) .^( -Plen));
    if(length(kv)<=1)
       % EE_prediction_vector = ((((P_s./2)+ ((1 - BER_rnd) .^( -Plen))).*Plen.*(sim.power_bit))+(sim.power_sense))./(k.*sim.N);
        EE_prediction_vector = (((((1 - BER_rnd) .^( -Plen))).*(P_s./2+1).*Plen.*(sim.power_bit))+(sim.power_sense))./(k.*sim.N);
    else
        %   EE_prediction_vector(j) = ((((((1 - BER_rnd) .^( -Plen)))).*Plen.*(sim.power_bit))+(sim.power_sense))./k.*sim.N;
      %  EE_prediction_vector(j) = ((((P_s./2)+ ((1 - BER_rnd) .^( -Plen))).*Plen.*(sim.power_bit))+(sim.power_sense))./(k.*sim.N);
% try
%         EE_prediction_vector(j) = (((((1 - BER_rnd) .^( -Plen))).*(P_s./2+1).*Plen.*(sim.power_bit))+(sim.power_sense))./(k.*sim.N);
% catch
 EE_prediction_vector(j,:) = (((((1 - BER_rnd) .^( -Plen))).*(P_s./2+1).*Plen.*(sim.power_bit))+(sim.power_sense))./(k.*sim.N);

    end
end
coeff1 = (((1 - BER_rnd) .^( -Plen))).*(P_s1./2+1).*Plen;
P_s1(end);
P_s(end);


coeff1(end);
% P_s = 2*(pp(1)*N(1,2)+pp(2)*N(2,2)+pp(3)*N(3,2))
% pp(1)*N(1,2)
% pp(2)*N(2,2)
% pp(3)*N(3,2)
% pp(2)*N(2,2)+pp(3)*N(3,2)
ali=1;
end