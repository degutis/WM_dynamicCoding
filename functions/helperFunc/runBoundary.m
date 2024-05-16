function runBoundary(cluster_above_chance,colorPlot,lineWidth)

    if nargin==2
        colorPlotVal = colorPlot;
        lineWidthVal = 1;
    elseif nargin==3
        colorPlotVal = colorPlot;
        lineWidthVal =lineWidth;
    else
        colorPlotVal='k';
        lineWidthVal =1;
    end
    
    boundaries = bwboundaries(imresize(logical(cluster_above_chance'),32));
    numBound = size(boundaries,1);
    for k=1:numBound
        curBound = boundaries{k};
        curBound = round(curBound/32);
        curBound(:,1) = curBound(:,1)+0.5;
        curBound(:,2) = curBound(:,2)+0.5;
        hold on;
        p=plot(curBound(:,2), curBound(:,1),'Color',colorPlotVal,'LineWidth',lineWidthVal);
        p.LineStyle = '-';
    end

end