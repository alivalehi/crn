function [tio, ni]= combine_arriv_depart(ti, to)

%This function combines arrival and departure times to 
% and the result is time ephocs and number of packets at the time
    
N1 = length(ti);
N2= length(to);

ti1 = 2;
to1 = 1;
tio1= 1; %instants of arrivals and departures

ni(tio1) = 1; % arrival at origin
tio(tio1) = ti(1);

not_finished = 1;
while ((ti1 <= N1) && (to1 <= N2))&& not_finished
    fprintf('.');
    if (ti1 == N1) && (to1 == N2)
        disp('lastround');
        not_finished= 0
       ti(ti1) 
       to(to1)
    end
    
    if ti(ti1) == to(to1)  %simultaneous arrival and departure
        tio1 =tio1+1;
        ni(tio1) = ni(tio1-1); 
        tio(tio1) = to(to1);
        ti1 = ti1+1;
        to1 = to1+1;
    elseif ti(ti1) < to(to1)  %arrival
        tio1 =tio1+1;
        ni(tio1) = ni(tio1-1) + 1; 
        tio(tio1) = ti(ti1);
        ti1= ti1+1;

        
    else  %departure
        tio1 =tio1+1;
        ni(tio1) = ni(tio1-1) - 1; 
        tio(tio1) = to(to1);
        to1= to1+1;
    end
    if ti1 > N1 %only departuters left
        n=length(to)-to1+1;
        tio = [tio, to(to1:end)]; 
        ni = [ni , ni(end) - 1: -1 : ni(end) - n];
        ti1 = N1+1;
    end
    if to1 > N2 %only departuters left
        n=length(ti)-ti1+1;
        tio = [tio, ti(ti1:end)]; 
        ni = [ni , ni(end) + 1: 1 : ni(end) + n];
        to1 = N2+1;
    end
    
    if any(ni<0), fprintf('Potentially wrong arrival or departure processes: departure event occurs without any waiting item!!!\n'); end
end
