%%%%% making a movie of the given snowfall %%%%%


%     titleOption = 1;
    
    % defaults for plotting 
    fontSize=14; lineWidth=2; markerSize=7; 
    set(0,'DefaultLineMarkerSize',markerSize);
    set(0,'DefaultLineLineWidth',lineWidth);
    set(0,'DefaultAxesFontSize',fontSize);
    set(0,'DefaultLegendFontSize',fontSize);


%     run(filename);

    v = VideoWriter('snowMovie.avi');
    open(v);
    
    for n=1:length(t)
%         plot(sfLocation(:,1,n),sfLocation(:,2,n),'ko','MarkerSize',7,'MarkerFaceColor','w');
        plot(sfLocation(:,1,n),sfLocation(:,2,n),...
            'w*','MarkerSize',8,'LineWidth',2);
%         title(sprintf('snowflakes: t=%2.2f',t(1,n)));
        xlim([.1,.9]); ylim([.1,.9]);
%         xlabel('x');ylabel('y or z');
        
        axis tight; 
        hold on;
        I = imread('xmastree.jpeg'); 
        h = image([.1,.9],[.1,.9],I); 
        uistack(h,'bottom');
        xlim([.1,.9]); ylim([.1,.9]);
        hold off;
        drawnow;
        
        frame = getframe(gcf);
        writeVideo(v,frame);
    end



    

    close(v);

