function [] = ali_plotdef(square_window)
set(gcf, 'WindowStyle', 'undocked') 
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
set(gca,'fontsize',22)
set(findall(gca, 'Type', 'Line'),'LineWidth',3);
end
