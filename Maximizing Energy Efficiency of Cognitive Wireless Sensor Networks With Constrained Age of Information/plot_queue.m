function [] = plot_queue(queue, Kind) %(tp, kp, rv, sv, wnv, sv2, wnv2, cht, chf) 
    %nplot = min(Nt,100);
    kplot=min(length(queue.kp), 10);
    
    figure;
%     subplot(221);
%     stem(t(1:nplot),ones(1,nplot));
%     title('original arrival symbols');
% 
%     subplot(222);
%     splot(t(1:nplot),[1:nplot], 'b'); 
%     title('Number of accumulated arrival symbols');
%     
%     subplot(223);
    stem(queue.tp(1:kplot),queue.kp(1:kplot));
    title('Packets ');

    subplot(224);
    splot(queue.tp(1:kplot),[1:kplot], 'b'); 
    title('Packets');
    


    plt(Kind) = stem(queue.wnv);
    plt(Kind) = plot(queue.wnv);
    lgn{Kind} = sprintf('T=%1.3f / L', T*a );
    hold on;
    pause;

    xv = (0: 0.1*T : 5*T);
    hist (queue.wnv, xv);



    plt(Kind) = stem(queue.wnv);
    plt(Kind) = plot(queue.wnv);
    lgn{Kind} = sprintf('K=%d', K );
    hold on; %pause;

end