function [pval_WMspace,WMspace]=projectDistractor(numRois,my_rois_plus,resultsSave,currentWM_orDist_cell,currentWM_distractor_cell)

%% Projection distractor

numSubs = size(currentWM_orDist_cell,1);
numIter = size(currentWM_distractor_cell,4);
numBins = size(currentWM_distractor_cell{1},1);

orWM_space_or = zeros(numBins,2,numRois);
orWM_space_distractor = zeros(numBins,2,numRois);
distractor_space_distractor = zeros(numBins,2,numRois);
distractor_space_or = zeros(numBins,2,numRois);

for rr = 1:numRois   
    currentROI_WM_train = cat(2,currentWM_orDist_cell{:,rr,1,1});
    currentROI_WM_test = cat(2,currentWM_orDist_cell{:,rr,2,1});

    currentROI_distractor_train = cat(2,currentWM_distractor_cell{:,rr,1,1});
    currentROI_distractor_test = cat(2,currentWM_distractor_cell{:,rr,2,1});

    currentROI_WM_train = mean(currentROI_WM_train,3);
    currentROI_WM_train = currentROI_WM_train - mean(currentROI_WM_train,1);
    
    currentROI_distractor_train = mean(currentROI_distractor_train,3);
    currentROI_distractor_train = currentROI_distractor_train - mean(currentROI_distractor_train,1);

    [~,W]=pca(currentROI_WM_train','NumComponents',2);
    [~,W_distractor]=pca(currentROI_distractor_train','NumComponents',2);
    
    currentROI_distractor_test = mean(currentROI_distractor_test,3);
    currentROI_distractor_test = currentROI_distractor_test - mean(currentROI_distractor_test,1);           

    currentROI_WM_test = mean(currentROI_WM_test,3);
    currentROI_WM_test = currentROI_WM_test - mean(currentROI_WM_test,1);           
        
    orWM_space_distractor(:,:,rr) = currentROI_distractor_test*W;
    orWM_space_or(:,:,rr) = currentROI_WM_test*W;

    distractor_space_distractor(:,:,rr) = currentROI_distractor_test*W_distractor;
    distractor_space_or(:,:,rr) = currentROI_WM_test*W_distractor;
end 


for roi=1:numRois       

    current_orWM_space_distractor = orWM_space_distractor(:,:,roi);
    current_orWM_space_or = orWM_space_or(:,:,roi);
    
    current_combo  = cat(3,current_orWM_space_distractor,current_orWM_space_or);
    maxTwoDim_dist = max(max(abs(current_combo),[],3));

    for dim=1:2
        current_orWM_space_distractor(:,dim) = current_orWM_space_distractor(:,dim)/maxTwoDim_dist(dim);
        current_orWM_space_or(:,dim) = current_orWM_space_or(:,dim)/maxTwoDim_dist(dim);
    end
    current_orWM_space_distractor = current_orWM_space_distractor*10;        
    current_orWM_space_or = current_orWM_space_or*10;       

    orWM_space_distractor(:,:,roi)=current_orWM_space_distractor;
    orWM_space_or(:,:,roi)=current_orWM_space_or;


    current_distractor_space_distractor = distractor_space_distractor(:,:,roi);   
    current_distractor_space_or = distractor_space_or(:,:,roi);
    
    current_combo_2  = cat(3,current_distractor_space_distractor,current_distractor_space_or);
    maxTwoDim_dist_2 = max(max(abs(current_combo_2),[],3));

    for dim=1:2
        current_distractor_space_distractor(:,dim) = current_distractor_space_distractor(:,dim)/maxTwoDim_dist_2(dim);
        current_distractor_space_or(:,dim) = current_distractor_space_or(:,dim)/maxTwoDim_dist_2(dim);
    end
    
    current_distractor_space_distractor = current_distractor_space_distractor*10;
    distractor_space_distractor(:,:,roi)=current_distractor_space_distractor;

    current_distractor_space_or = current_distractor_space_or*10;
    distractor_space_or(:,:,roi)=current_distractor_space_or;

end

colors = cbrewer('qual', 'Set1', 6);

fig2 = figure(2);    
for a= 1:numRois
    set(fig2, 'PaperUnits', 'inches');
    x_width=10 ;y_width=10;
    set(fig2, 'PaperPosition', [0 0 x_width y_width]); %
    
    subplot(2,numRois,a)
    
    % Distractor in distractor
    
    currentPlot = polyshape(squeeze(distractor_space_distractor(:,:,a)));
    pg2 = plot(currentPlot);
    pg2.FaceAlpha = 0.75;
    pg2.FaceColor = 'r';
    
    hold on
    
    p1=scatter(squeeze(distractor_space_distractor(1,1,a)),squeeze(distractor_space_distractor(1,2,a)),5,'filled','HandleVisibility','off');
    p1.MarkerFaceColor = colors(1,:);
    p1.MarkerEdgeColor = colors(1,:);
    p1.MarkerFaceAlpha = 1;
    p1.LineWidth = 2;

    hold on 

    p2=scatter(squeeze(distractor_space_distractor(2,1,a)),squeeze(distractor_space_distractor(2,2,a)),5,'filled','HandleVisibility','off');
    p2.MarkerFaceColor = colors(2,:);
    p2.MarkerEdgeColor = colors(2,:);
    p2.MarkerFaceAlpha = 1;
    p2.LineWidth = 2;

    hold on 
    
    p3=scatter(squeeze(distractor_space_distractor(3,1,a)),squeeze(distractor_space_distractor(3,2,a)),5,'filled','HandleVisibility','off');
    p3.MarkerFaceColor = colors(3,:);
    p3.MarkerEdgeColor = colors(3,:);
    p3.MarkerFaceAlpha = 1;
    p3.LineWidth = 2;

    hold on 
    
    p4=scatter(squeeze(distractor_space_distractor(4,1,a)),squeeze(distractor_space_distractor(4,2,a)),5,'filled','HandleVisibility','off');
    p4.MarkerFaceColor = colors(4,:);
    p4.MarkerEdgeColor = colors(4,:);
    p4.MarkerFaceAlpha = 1;
    p4.LineWidth = 2;
    
    % Wm in distractor
    
    currentPlot = polyshape(squeeze(distractor_space_or(:,:,a)));
    pg2 = plot(currentPlot);
    pg2.FaceAlpha = 0.75;
    pg2.FaceColor = 'b';

    hold on
        
    p1=scatter(squeeze(distractor_space_or(1,1,a)),squeeze(distractor_space_or(1,2,a)),5,'filled','HandleVisibility','off');
    p1.MarkerFaceColor = colors(1,:);
    p1.MarkerEdgeColor = colors(1,:);
    p1.MarkerFaceAlpha = 1;
    p1.LineWidth = 2;

    hold on 

    p2=scatter(squeeze(distractor_space_or(2,1,a)),squeeze(distractor_space_or(2,2,a)),5,'filled','HandleVisibility','off');
    p2.MarkerFaceColor = colors(2,:);
    p2.MarkerEdgeColor = colors(2,:);
    p2.MarkerFaceAlpha = 1;
    p2.LineWidth = 2;

    hold on 
    
    p3=scatter(squeeze(distractor_space_or(3,1,a)),squeeze(distractor_space_or(3,2,a)),5,'filled','HandleVisibility','off');
    p3.MarkerFaceColor = colors(3,:);
    p3.MarkerEdgeColor = colors(3,:);
    p3.MarkerFaceAlpha = 1;
    p3.LineWidth = 2;

    hold on 
    
    p4=scatter(squeeze(distractor_space_or(4,1,a)),squeeze(distractor_space_or(4,2,a)),5,'filled','HandleVisibility','off');
    p4.MarkerFaceColor = colors(4,:);
    p4.MarkerEdgeColor = colors(4,:);
    p4.MarkerFaceAlpha = 1;
    p4.LineWidth = 2;
    
    axis square
    title(my_rois_plus{a},' distractor')

    xlim([-10, 10])
    ylim([-10, 10])

    
    subplot(2,numRois,numRois+a)

    % WM in WM
    
    currentPlot = polyshape(squeeze(orWM_space_or(:,:,a)));
    pg2 = plot(currentPlot,"LineWidth",1,'HandleVisibility','off');
    pg2.FaceAlpha = 0.75;
    pg2.FaceColor = 'b';
    axis square
    hold on 

    hold on

    p1=scatter(squeeze(orWM_space_or(1,1,a)),squeeze(orWM_space_or(1,2,a)),5,'filled','HandleVisibility','off');
    p1.MarkerFaceColor = colors(1,:);
    p1.MarkerEdgeColor = colors(1,:);
    p1.MarkerFaceAlpha = 1;
    p1.LineWidth = 2;

    hold on 

    p2=scatter(squeeze(orWM_space_or(2,1,a)),squeeze(orWM_space_or(2,2,a)),5,'filled','HandleVisibility','off');
    p2.MarkerFaceColor = colors(2,:);
    p2.MarkerEdgeColor = colors(2,:);
    p2.MarkerFaceAlpha = 1;
    p2.LineWidth = 2;

    hold on 
    
    p3=scatter(squeeze(orWM_space_or(3,1,a)),squeeze(orWM_space_or(3,2,a)),5,'filled','HandleVisibility','off');
    p3.MarkerFaceColor = colors(3,:);
    p3.MarkerEdgeColor = colors(3,:);
    p3.MarkerFaceAlpha = 1;
    p3.LineWidth = 2;

    hold on 
    
    p4=scatter(squeeze(orWM_space_or(4,1,a)),squeeze(orWM_space_or(4,2,a)),5,'filled','HandleVisibility','off');
    p4.MarkerFaceColor = colors(4,:);
    p4.MarkerEdgeColor = colors(4,:);
    p4.MarkerFaceAlpha = 1;
    p4.LineWidth = 2;

    hold on

    % Distractor in WM

    currentPlot = polyshape(squeeze(orWM_space_distractor(:,:,a)));
    pg = plot(currentPlot,"LineWidth",1,'HandleVisibility','off');
    pg.FaceAlpha = 0.75;
    pg.FaceColor = 'r';

    p11=scatter(squeeze(orWM_space_distractor(1,1,a)),squeeze(orWM_space_distractor(1,2,a)),5,'filled','HandleVisibility','off');
    p11.MarkerFaceColor = colors(1,:);
    p11.MarkerEdgeColor = colors(1,:);
    p11.MarkerFaceAlpha = 1;
    p11.LineWidth = 2;

    hold on 

    p22=scatter(squeeze(orWM_space_distractor(2,1,a)),squeeze(orWM_space_distractor(2,2,a)),5,'filled','HandleVisibility','off');
    p22.MarkerFaceColor = colors(2,:);
    p22.MarkerEdgeColor = colors(2,:);
    p22.MarkerFaceAlpha = 1;
    p22.LineWidth = 2;

    hold on 
    
    p33=scatter(squeeze(orWM_space_distractor(3,1,a)),squeeze(orWM_space_distractor(3,2,a)),5,'filled','HandleVisibility','off');
    p33.MarkerFaceColor = colors(3,:);
    p33.MarkerEdgeColor = colors(3,:);
    p33.MarkerFaceAlpha = 1;
    p33.LineWidth = 2;

    hold on 
    
    p44=scatter(squeeze(orWM_space_distractor(4,1,a)),squeeze(orWM_space_distractor(4,2,a)),5,'filled','HandleVisibility','off');
    p44.MarkerFaceColor = colors(4,:);
    p44.MarkerEdgeColor = colors(4,:);
    p44.MarkerFaceAlpha = 1;
    p44.LineWidth = 2;

    xlim([-10, 10])
    ylim([-10, 10])
    title(my_rois_plus{a},' WM space')
end

saveas(fig2,fullfile(resultsSave,'distractorSpace.svg'))
close(fig2)


%% per subject

Theta = zeros(size(currentWM_orDist_cell,1),numRois,numIter);

for b=1:numIter
    for sub=1:size(currentWM_orDist_cell,1)
        for rr = 1:numRois   
            Theta(sub,rr,b) = runPCAangle_OrientationDistractor_splitHalf(currentWM_orDist_cell,currentWM_distractor_cell,sub,rr,b);
        end
    end
end

WMspace = mean(Theta,3);

pval_WMspace = runStatsPermutationAgainstNull(WMspace,[],[],1,4,14011);

fig1 = figure(1)    
set(fig1, 'PaperUnits', 'inches');
x_width=3 ;y_width=3;
set(fig1, 'PaperPosition', [0 0 x_width y_width]); %

for rr=1:numRois
    h1 = ploterr(rr-0.15, mean(WMspace(:,rr)), [], std(WMspace(:,rr))/sqrt(numSubs), 'k.', 'abshhxy', 0);
    set(h1,'Color','b')
    hold on
    mysigstar(gca,rr-0.15,15,pval_WMspace(rr))
end
box off
ylim([-5 15])
xlim([0 numRois+0.3])
xticks([1:8])
xticklabels(my_rois_plus)
hold on

saveas(fig1,fullfile(resultsSave,'anglesDistractor_OrWM.svg'))
close(fig1)
