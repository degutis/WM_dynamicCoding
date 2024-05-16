function plotSubspaceTriangleMultipleROI(noDist,orDist,noiseDist,angles,timePointIndices,namesROI,alpha,resultsSave,nameSave)

    numRois  = size(noDist,3);
    indexMatrix = nan(size(noDist,1),size(noDist,2));
    
    for t=1:length(timePointIndices)
        indexMatrix(timePointIndices{t}(1),timePointIndices{t}(2))=1;
    end
    
    indexMatrix = repmat(indexMatrix,1,1,numRois);
    
    noDist = indexMatrix.*noDist;
    orDist = indexMatrix.*orDist;
    noiseDist = indexMatrix.*noiseDist;
    
    idxNoDist = noDist<alpha;
    idxOrDist = orDist<alpha;
    idxNoiseDist = noiseDist<alpha;
    
    titleCombo = {'No Distractor','Orientation Distractor','Noise Distractor'};
    
    colorsMap = cbrewer('seq', 'Oranges', 50);
    colorsMap = colorsMap(10:50,:);
    
    colors = cbrewer('div','PRGn',11);
    colors = [colors(3,:);colors(8,:);colors(10,:)];
    
    
    fig=figure(1)
    set(fig, 'PaperUnits', 'inches');
    x_width=14 ;y_width=7;
    set(fig, 'PaperPosition', [0 0 x_width y_width]); %
    for r = 1:numRois
        subplot(numRois,3,(r-1)*3+1)
        h=imagesc(angles(:,:,r,1))  
        set(h, 'AlphaData', ~isnan(noDist(:,:,r)))
        set(gca, 'ydir', 'normal');
        axis square;
        hold on
        for a=1:size(noDist,1)
            for b=1:size(noDist,2)
                if idxNoDist(a,b,r)==1
                    plot(b,a,'*c')
                    hold on
                end
            end
        end
    
        title([namesROI{r},': ',titleCombo{1}],'Fontsize',8)
        colormap(colorsMap)
        set(gca,'TickDir','out','XTick',1:6,'XTickLabel',{'','7.2s','','12s','','16.8s'},'Fontsize',5);
        set(gca,'TickDir','out','YTick',1:6,'YTickLabel',{'4.8s','7.2s','9.6s','12s','14.4s','16.8s'},'Fontsize',5);
        handles = colorbar;
        caxis([-0.3 8.5])
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'Principal Angle';
        handles.Label.FontSize = 5;
        drawnow;
    
        subplot(numRois,3,(r-1)*3+3)
        h=imagesc(angles(:,:,r,2))  
        set(h, 'AlphaData', ~isnan(orDist(:,:,r)))
        set(gca, 'ydir', 'normal');
        axis square;
        hold on
        for a=1:size(noDist,1)
            for b=1:size(noDist,2)
                if idxOrDist(a,b,r)==1
                    plot(b,a,'*c')
                    hold on
                end
            end
        end
    
        title([namesROI{r},': ',titleCombo{2}],'Fontsize',8)
        colormap(colorsMap)
        set(gca,'TickDir','out','XTick',1:6,'XTickLabel',{'','7.2s','','12s','','16.8s'},'Fontsize',5);
        set(gca,'TickDir','out','YTick',1:6,'YTickLabel',{'4.8s','7.2s','9.6s','12s','14.4s','16.8s'},'Fontsize',5);
        handles = colorbar;
        caxis([-0.3 8.5])
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'Principal Angle';
        handles.Label.FontSize = 5;
        drawnow;
        
        subplot(numRois,3,(r-1)*3+2)
        h=imagesc(angles(:,:,r,3))  
        set(h, 'AlphaData', ~isnan(noiseDist(:,:,r)))
        set(gca, 'ydir', 'normal');
        axis square;
        hold on
        for a=1:size(noDist,1)
            for b=1:size(noDist,2)
                if idxNoiseDist(a,b,r)==1
                    plot(b,a,'*c')
                    hold on
                end
            end
        end
    
        title([namesROI{r},': ',titleCombo{3}],'Fontsize',8)
        colormap(colorsMap)
        set(gca,'TickDir','out','XTick',1:6,'XTickLabel',{'4.8s','','9.6s','','14.4s',''},'Fontsize',5);
        set(gca,'TickDir','out','YTick',1:6,'YTickLabel',{'4.8s','7.2s','9.6s','12s','14.4s','16.8s'},'Fontsize',5);
        handles = colorbar;
        caxis([-0.3 8.5])
        handles.TickDirection = 'out';
        handles.Box = 'off';
        handles.Label.String = 'Principal Angle';
        handles.Label.FontSize = 5;
        drawnow;
    
    end
    
    saveas(fig,fullfile(resultsSave,nameSave))
    close(fig)

end