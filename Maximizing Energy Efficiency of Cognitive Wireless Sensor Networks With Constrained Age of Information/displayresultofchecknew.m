 L=sim2(sim2~=0);
                L2=sim3(sim3~=0);
                opt1 = find(status_flag==2);
                fprintf('Number of arrive at busy time slot:%f  Analytical:%f \n',mean(P(opt1)),(1/p)-a);
                opt2 = find(status_flag==3);
                fprintf('Number of arrive at small available time slot:%f  Analytical:%f\n',mean(P(opt2)),(1/p)-a);    
                opt3 = find(status_flag==1);
                fprintf('mean of service time of available time slots:%f  Analytical:%f\n', mean(sim1(opt3)),s1);   
                opt4 = find(status_flag==2);
                fprintf('mean of service time of busy time slots:%f  Analytical:%f\n', mean(sim1(opt4)),ES_BusyArrive+s1 + u/2); 
                opt5 = find(status_flag==3);
                fprintf('mean of service time of small available time slots:%f  Analytical:%f\n', mean(sim1(opt5)),ES_AvArrive+s1+mean(L2) );
                fprintf('total service time A:%f S:%f \n', ESk ,mean(sim1));
                fprintf('scenario %f %f %f\n',Scenario1,Scenario2,Scenario3);
                last_slot2 = last_slot;
                ww = find(first_slot~=0);
                last_slot2 = last_slot2(ww);
                qq = find(first_slot2~=0);
                last_slot = last_slot(qq);
                first_slot = first_slot(find(first_slot~=0));
                first_slot2 = first_slot2(find(first_slot2~=0));
               fprintf('mean of service time of small available time includes p=1  :%f  Analytical:%f\n', mean(slot(last_slot2)-slot(first_slot)),ES_AvArrive)
                fprintf('mean of service time of small available time slots from end of end of first busy to last moment:%f  Analytical:%f\n', mean(slot(last_slot)-slot(first_slot2)),ES_AvArrive-u )
               % fprintf('mean of service time of small available time slots from end of slot that we arrive to end of next slot(first busy):%f  Analytical:%f\n',mean(slot(first_slot2)-slot(first_slot)),u )
               % mean(slot(first_slot2)-slot(first_slot))+s1+mean(P(opt2))*6                                     
                
                 fprintf('result %f %f\n',mean(L),mean(L2));

                 fprintf('Si analysis: v:%f u:%f   s1:%f   ES:[%f  %f],  ES2::[%f  %f],var(s):[%f  %f] \n', v,u,s1, ESk, mean(sim1), ESk2, mean(sim1.^2), varsk2,var(sim1));