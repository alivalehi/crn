function [] = plotdef(square_window)
%close all



start_x = 200;%500;
start_y = 100;%300;


if square_window == 1
    width = 750*.8;
    height = 650*.8;
    defpos = [start_x, start_y, width, height];  %start x    start y    width height
    set(gcf,'position',defpos);  
elseif square_window ==  2
    scrsz = get(0,'ScreenSize');
    defpos = scrsz + [50 50 -100 -150];
    set(gcf,'position',defpos);  
elseif square_window == 3 %full screen
    
    scrsz = get(0,'ScreenSize');
    set(gcf,'Position',[scrsz(1) scrsz(2) scrsz(3) scrsz(4)])    
else
    
    width = 750;
    height = 550;
    defpos = [start_x, start_y, width, height];  %start x    start y    width height
    set(gcf,'position',defpos);  
end

set(gcf,'defaultaxesfontsize',14);
set(gcf,'defaultaxesfontsize',14);
set(gcf,'defaultaxeslinewidth',1);
set(gcf,'defaultlinelinewidth',2.5);


%AX=legend('One');
%LEG = findobj(AX,'type','text');
%set(LEG,'FontSize',12)
%set(gca,'FontSize',16,'FontName','Times') % changes the size and font of the axis numbers
%set(gca,'XLim',[0 5],'YLim',[0 1]); % changes the x and y plot range
%set(gca,'XTick',[0:2:10],'XGrid','on'); % sets new ticks on the x-axis and adds grid lines


end
