            Av_mean = 10; 
            Bu_mean = 40;
            k = 20; % max number of symbol per package
            n=1e6;
            lambda = 30;
            for l = 5:k
                inter_arrival = poissrnd(lambda,1,n);
                number_of_packet = floor(n/l);
                tp = zeros(1,number_of_packet);
                tp_s(1) = 0;
                for j = 1:number_of_packet-1
                    tp_s(j+1) = sum(inter_arrival(1:l*j)) ;
                    tp_e(j) = tp_s(j) + sum(inter_arrival(((j*l)+1):((j*l)+(l-1))));
                end
                pause;
                av_time = exprnd(Av_mean,1,n);
                bu_time = exprnd(Bu_mean,1,n);
                channel = reshape([av_time; bu_time], [1, 2*n]);
                chf = reshape([ones(1,n); zeros(1,n)], [1, 2*n]);
                slot(1) = 0;
                for i=1:2*n
                    slot(i+1) = slot(i) + channel(i);
                end
                CRNwnv = zeros(1,npackets); CRNsv = zeros(1,npackets); 
               CRNwnv(1) = 0; 
                for i=1:length(tp_e)
                 sst0 = tp_s(i) + CRNwnv(i);   %service start time for all retransmissions
                 sst = sst0;               %service start time
                 CRNsv(i) = sst - sst0;
                 inter_pack_time = tp(p+1)-tp(p);
                 CRNwnv(i+1) = max(0, CRNwnv(i) + CRNsv(i) - inter_pack_time);  %inter-arrival time
                 current_pointer = sst;
                 next_start =  min(find(slot>current_pointer));
                 for j=next_start:length(slot)
                    previous_end = max(find(slot<=current_pointer));
                    next_start =  min(find(slot>current_pointer));
               
                if(chf(next_start)==1)% busy
                    service_time(i) =service_time(i)+ (next_start-current_pointer);
                    current_pointer = next_start
                else % available
                    if((next_start-current_pointer)<s1)%not enough length
                        service_time(i) = service_time(i)+ (next_start-current_pointer);
                        current_pointer = next_start
                        current_length = next_start- slot(i)
                    else
                        service_time(i) = service_time(i)+ s1;
                        break;
                    end
               end
                end   
            end
            end   