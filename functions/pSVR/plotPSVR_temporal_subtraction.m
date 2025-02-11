function plotPSVR_temporal_subtraction(experiment,expandIPS)

global path2data

selectTime =1:21;
timeDim = length(selectTime);

nameAdd = 'ReliableVox';

resultsSave = fullfile(path2data,'pSVR_analysis',experiment,'temporal',nameAdd);
dirReliable = dir(fullfile(resultsSave,'temporal_runwise_PSVM_*mat'));

mkdir(resultsSave)

if expandIPS==1
    my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS','LO1','LO2'};
else
    my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
end

if strcmp('E1',experiment)
    N=6;
    dirReliable2 = dir(fullfile(resultsSave,'temporal_Extra_distractor_runwise_*mat'));
elseif strcmp('E2',experiment)
    N=7;
    dirReliable2 = dir(fullfile(resultsSave,'temporal_Extra_*mat'));
end

Roi=numel(my_rois_plus);

wmDist_reliable = cell(N,Roi);
wmDist_orOnDist_reliable = cell(N,Roi);
wmDist_distOnOr_reliable = cell(N,Roi);


wmNoDist_reliable = cell(N,Roi);

wmOrDist_reliable = cell(N,Roi);
wmOrDistonMem_reliable = cell(N,Roi);

wmNoiseDist_reliable = cell(N,Roi);
wmNoiseDistonMem_reliable = cell(N,Roi);


wmDist_reliable_MEAN = cell(N,1);
wmDist_orOnDist_reliable_MEAN = cell(N,1);
wmDist_distOnOr_reliable_MEAN = cell(N,1);

wmNoDist_reliable_MEAN = cell(N,1);
wmDist_reliable_MEAN = cell(N,1);

wmOrDist_reliable_MEAN = cell(N,1);
wmOrDistonMem_reliable_MEAN = cell(N,1);

wmNoiseDist_reliable_MEAN = cell(N,1);
wmNoiseDistonMem_reliable_MEAN = cell(N,1);


for f = 1:numel(dirReliable)
    subNow = load(fullfile(dirReliable(f).folder,dirReliable(f).name));
    subNow_extra = load(fullfile(dirReliable2(f).folder,dirReliable2(f).name));
    
    wmNoDist_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow.wmnoDist_fcaMatrix(f,:), 'UniformOutput', false);
    wmOrDist_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow.wmOr_fcaMatrix(f,:), 'UniformOutput', false);
    wmOrDistonMem_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow.wmOronMem_fcaMatrix(f,:), 'UniformOutput', false);       
    wmNoiseDist_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow.wmNoise_fcaMatrix(f,:), 'UniformOutput', false);          
    wmNoiseDistonMem_reliable(f,:) = subNow.wmNoiseonMem_fcaMatrix(f,:);
    wmNoiseDistonMem_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow.wmNoiseonMem_fcaMatrix(f,:), 'UniformOutput', false);          

    wmDist_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow_extra.wmDist_fcaMatrix_distOnDist(f,:), 'UniformOutput', false);
    wmDist_orOnDist_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow.wmDist_fcaMatrix(f,:), 'UniformOutput', false);
    wmDist_distOnOr_reliable(f,:) = cellfun(@(x) mean(x(selectTime,selectTime,:),3), subNow_extra.wmDist_fcaMatrix_genDisttoOr(f,:), 'UniformOutput', false);
    
end

Cluster_noDist = zeros(timeDim,timeDim,Roi);

Cluster_Dist = zeros(timeDim,timeDim,Roi);
Cluster_OronDist = zeros(timeDim,timeDim,Roi);
Cluster_DistonOr = zeros(timeDim,timeDim,Roi);

Cluster_orDist = zeros(timeDim,timeDim,Roi);
Cluster_noiseDist = zeros(timeDim,timeDim,Roi);
Cluster_orDistonMem = zeros(timeDim,timeDim,Roi);
Cluster_noiseDistonMem = zeros(timeDim,timeDim,Roi);

if isfile(fullfile(resultsSave,['temporalClusterPerm_',nameAdd,'.mat']))
    load(fullfile(resultsSave,['temporalClusterPerm_',nameAdd,'.mat']))
else
    for r = 1:Roi
        disp(['Roi: ',num2str(r)])

        % P-values

        wmNoDist_Pval =reshape(cell2mat(wmNoDist_reliable(:,r)'),[timeDim timeDim N]);    
        
        wmDist_Pval =reshape(cell2mat(wmDist_reliable(:,r)'),[timeDim timeDim N]);
        wmOronDist_Pval =reshape(cell2mat(wmDist_orOnDist_reliable(:,r)'),[timeDim timeDim N]);
        wmDistonOr_Pval =reshape(cell2mat(wmDist_distOnOr_reliable(:,r)'),[timeDim timeDim N]);       
        
        wmOrDist_Pval =reshape(cell2mat(wmOrDist_reliable(:,r)'),[timeDim timeDim N]);
        wmNoiseDist_Pval =reshape(cell2mat(wmNoiseDist_reliable(:,r)'),[timeDim timeDim N]);
        wmOrDistonMem_Pval =reshape(cell2mat(wmOrDistonMem_reliable(:,r)'),[timeDim timeDim N]);
        wmNoiseDistonMem_Pval =reshape(cell2mat(wmNoiseDistonMem_reliable(:,r)'),[timeDim timeDim N]);


        Cluster_noDist(:,:,r)=permutationTest_cluster_psvr(wmNoDist_Pval-50);
        
        Cluster_Dist(:,:,r)=permutationTest_cluster_psvr(wmDist_Pval-50);
        Cluster_OronDist(:,:,r)=permutationTest_cluster_psvr(wmOronDist_Pval-50);
        Cluster_DistonOr(:,:,r)=permutationTest_cluster_psvr(wmDistonOr_Pval-50);

        Cluster_orDist(:,:,r)=permutationTest_cluster_psvr(wmOrDist_Pval-50);
        Cluster_noiseDist(:,:,r)=permutationTest_cluster_psvr(wmNoiseDist_Pval-50);
        Cluster_orDistonMem(:,:,r)=permutationTest_cluster_psvr(wmOrDistonMem_Pval-50);
        Cluster_noiseDistonMem(:,:,r)=permutationTest_cluster_psvr(wmNoiseDistonMem_Pval-50);
        
    end
    save(fullfile(resultsSave,['temporalClusterPerm_',nameAdd,'.mat']),'Cluster_noDist','Cluster_Dist','Cluster_OronDist','Cluster_DistonOr','Cluster_orDist','Cluster_noiseDist','Cluster_orDistonMem','Cluster_noiseDistonMem');    
end


%% Dynamic coding

Dynamic_noDist = zeros(timeDim,timeDim,Roi);

Dynamic_Dist = zeros(timeDim,timeDim,Roi);
Dynamic_OronDist = zeros(timeDim,timeDim,Roi);
Dynamic_DistonOr = zeros(timeDim,timeDim,Roi);

Dynamic_orDist = zeros(timeDim,timeDim,Roi);
Dynamic_noiseDist = zeros(timeDim,timeDim,Roi);
Dynamic_orDistonMem = zeros(timeDim,timeDim,Roi);
Dynamic_noiseDistonMem = zeros(timeDim,timeDim,Roi);

if isfile(fullfile(resultsSave,['temporalDynamicPerm_',nameAdd,'.mat']))
    load(fullfile(resultsSave,['temporalDynamicPerm_',nameAdd,'.mat']))
else
    for r = 1:Roi
        disp(['Dynamic... Roi: ',num2str(r)])

        % P-values

        wmNoDist_Pval =reshape(cell2mat(wmNoDist_reliable(:,r)'),[timeDim timeDim N]);    

        wmDist_Pval =reshape(cell2mat(wmDist_reliable(:,r)'),[timeDim timeDim N]);
        wmOronDist_Pval =reshape(cell2mat(wmDist_orOnDist_reliable(:,r)'),[timeDim timeDim N]);
        wmDistonOr_Pval =reshape(cell2mat(wmDist_distOnOr_reliable(:,r)'),[timeDim timeDim N]);       
        
        wmOrDist_Pval =reshape(cell2mat(wmOrDist_reliable(:,r)'),[timeDim timeDim N]);
        wmNoiseDist_Pval =reshape(cell2mat(wmNoiseDist_reliable(:,r)'),[timeDim timeDim N]);
        wmOrDistonMem_Pval =reshape(cell2mat(wmOrDistonMem_reliable(:,r)'),[timeDim timeDim N]);
        wmNoiseDistonMem_Pval =reshape(cell2mat(wmNoiseDistonMem_reliable(:,r)'),[timeDim timeDim N]);
        
        Dynamic_noDist(:,:,r)=dynamicCoding_signInversion(wmNoDist_Pval);
        Dynamic_Dist(:,:,r)=dynamicCoding_signInversion(wmDist_Pval);
        Dynamic_OronDist(:,:,r)=dynamicCoding_signInversion(wmOronDist_Pval);
        Dynamic_DistonOr(:,:,r)=dynamicCoding_signInversion(wmDistonOr_Pval);

        Dynamic_orDist(:,:,r)=dynamicCoding_signInversion(wmOrDist_Pval);
        Dynamic_noiseDist(:,:,r)=dynamicCoding_signInversion(wmNoiseDist_Pval);
        Dynamic_orDistonMem(:,:,r)=dynamicCoding_signInversion(wmOrDistonMem_Pval);
        Dynamic_noiseDistonMem(:,:,r)=dynamicCoding_signInversion(wmNoiseDistonMem_Pval);

    end
    save(fullfile(resultsSave,['temporalDynamicPerm_',nameAdd,'.mat']),'Dynamic_noDist','Dynamic_Dist','Dynamic_OronDist','Dynamic_DistonOr','Dynamic_orDist','Dynamic_noiseDist','Dynamic_orDistonMem','Dynamic_noiseDistonMem');    
end







%% Mean
for r = 1:Roi
    wmNoDist_reliable_MEAN(r) = {mean(reshape(cell2mat(wmNoDist_reliable(:,r)'),[timeDim timeDim N]),3)};
    wmDist_reliable_MEAN(r) = {mean(reshape(cell2mat(wmDist_reliable(:,r)'),[timeDim timeDim N]),3)};
    wmDist_orOnDist_reliable_MEAN(r) = {mean(reshape(cell2mat(wmDist_orOnDist_reliable(:,r)'),[timeDim timeDim N]),3)};
    wmDist_distOnOr_reliable_MEAN(r) = {mean(reshape(cell2mat(wmDist_distOnOr_reliable(:,r)'),[timeDim timeDim N]),3)};

    wmOrDist_reliable_MEAN(r) = {mean(reshape(cell2mat(wmOrDist_reliable(:,r)'),[timeDim timeDim N]),3)};
    wmOrDistonMem_reliable_MEAN(r) = {mean(reshape(cell2mat(wmOrDistonMem_reliable(:,r)'),[timeDim timeDim N]),3)};

    wmNoiseDist_reliable_MEAN(r) = {mean(reshape(cell2mat(wmNoiseDist_reliable(:,r)'),[timeDim timeDim N]),3)};
    wmNoiseDistonMem_reliable_MEAN(r) = {mean(reshape(cell2mat(wmNoiseDistonMem_reliable(:,r)'),[timeDim timeDim N]),3)};
end


%% Calculate dynamicism index (Spaak et al 2017)

dynamicism_orDist = zeros(length(selectTime),Roi);
dynamicism_noDist= zeros(length(selectTime),Roi);
dynamicism_noiseDist= zeros(length(selectTime),Roi);

for r=1:Roi
    dynamicism_orDist(:,r) = dynamicism(Dynamic_orDist(:,:,r).*dynamicsMultiply(Cluster_orDist(:,:,r)));
    dynamicism_noDist(:,r) = dynamicism(Dynamic_noDist(:,:,r).*dynamicsMultiply(Cluster_noDist(:,:,r)));    
    dynamicism_noiseDist(:,r) = dynamicism(Dynamic_noiseDist(:,:,r).*dynamicsMultiply(Cluster_noiseDist(:,:,r)));
end

%% Plotting
Roi_plot =numel(my_rois_plus); 

fig11 = figure(11)
colors = cbrewer('div','PRGn',11);
colors = [colors(3,:);colors(8,:);colors(10,:)];
set(fig11, 'PaperUnits', 'inches');
x_width=5 ;y_width=20;
set(fig11, 'PaperPosition', [0 0 x_width y_width]); %

for rr = 1:Roi_plot
    subplot(8,1,rr)
    
    plot(selectTime,dynamicism_noDist(:,rr),'Color',colors(1,:),'LineWidth',4);
    hold on
    plot(selectTime,dynamicism_orDist(:,rr),'Color',colors(2,:),'LineWidth',4);
    hold on
    plot(selectTime,dynamicism_noiseDist(:,rr),'Color',colors(3,:),'LineWidth',4);
    
    lgnd{1} = sprintf('No distractor')
    lgnd{2} = sprintf('Orientation distractor')
    lgnd{3} = sprintf('Noise distractor')    
    
    title(my_rois_plus{rr})
    ylim([0 0.4])
    xlim([0 timeDim])
    box off
    ylabel('Dynamicism index') 
end

suplabel('Time (s)','x')
lh=legend(lgnd);
legnames = {'No distractor','Orientation distractor','Noise distractor'};
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
lpos = lh.Position;
lh.Position = lpos;

saveas(fig11,[resultsSave,'/',nameAdd,'_VWM_dynamicismIndex.svg'])    
close(fig11)       



%% Orientation

colorsMap = cbrewer('div', 'RdYlBu', 64);
colorsMap = flipud(colorsMap); 
colorsMap(colorsMap<0) = 0;
colorsMap(colorsMap>1)=1;


fig = figure(11)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmOrDist_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_orDist(:,:,rr))
    runBoundary(Dynamic_orDist(:,:,rr).*dynamicsMultiply(Cluster_orDist(:,:,rr)),'b')
    
    dotplot = Dynamic_orDist(:,:,rr).*dynamicsMultiply(Cluster_orDist(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end
    
    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})
    
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_OrDistractor_onlyTrainSelf.svg'])    
close(fig)    


fig = figure(11)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmOrDistonMem_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_orDistonMem(:,:,rr))
    runBoundary(Dynamic_orDistonMem(:,:,rr).*dynamicsMultiply(Cluster_orDistonMem(:,:,rr)),'b')

    dotplot = Dynamic_orDistonMem(:,:,rr).*dynamicsMultiply(Cluster_orDistonMem(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})

    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_OrDistractor_onlyTrainNoDist.svg'])    
close(fig)    

%% Noise

fig = figure(22)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmNoiseDist_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_noiseDist(:,:,rr))
    runBoundary(Dynamic_noiseDist(:,:,rr).*dynamicsMultiply(Cluster_noiseDist(:,:,rr)),'b')
    dotplot = Dynamic_noiseDist(:,:,rr).*dynamicsMultiply(Cluster_noiseDist(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})   
    
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_NoiseDistractor_onlyTrainSelf.svg'])    
close(fig)    


fig = figure(22)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmNoiseDistonMem_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_noiseDistonMem(:,:,rr))
    runBoundary(Dynamic_noiseDistonMem(:,:,rr).*dynamicsMultiply(Cluster_noiseDistonMem(:,:,rr)),'b')
    dotplot = Dynamic_noiseDistonMem(:,:,rr).*dynamicsMultiply(Cluster_noiseDistonMem(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})   
    
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_NoiseDistractor_onlyTrainNoDist.svg'])    
close(fig)    

%% No distractor

fig = figure(33)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmNoDist_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_noDist(:,:,rr))
    runBoundary(Dynamic_noDist(:,:,rr).*dynamicsMultiply(Cluster_noDist(:,:,rr)),'b')
    dotplot = Dynamic_noDist(:,:,rr).*dynamicsMultiply(Cluster_noDist(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end
    
    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})
   
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_NoDistractor_onlyTrainSelf.svg'])    
close(fig)    

%% Sensory distractor

fig = figure(44)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmDist_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_Dist(:,:,rr))
    runBoundary(Dynamic_Dist(:,:,rr).*dynamicsMultiply(Cluster_Dist(:,:,rr)),'b')
    dotplot = Dynamic_Dist(:,:,rr).*dynamicsMultiply(Cluster_Dist(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})    
    
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_Distractor.svg'])    
close(fig)       

%% Train on distractor test on WM

fig = figure(44)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmDist_distOnOr_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_DistonOr(:,:,rr))
    runBoundary(Dynamic_DistonOr(:,:,rr).*dynamicsMultiply(Cluster_DistonOr(:,:,rr)),'b')
    dotplot = Dynamic_DistonOr(:,:,rr).*dynamicsMultiply(Cluster_DistonOr(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})    
    
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_Distractor_trainDisttestOr.svg'])    
close(fig)       


%% Train on VWM test on distractor

fig = figure(44)
set(fig, 'PaperUnits', 'inches');
x_width=20 ;y_width=5;
set(fig, 'PaperPosition', [0 0 x_width y_width]); %
for rr = 1:Roi_plot
    subplot(2,4,rr)
    imagesc(wmDist_orOnDist_reliable_MEAN{rr}')
    set(gca, 'ydir', 'normal');
    axis square;

    runBoundary(Cluster_OronDist(:,:,rr))
    runBoundary(Dynamic_OronDist(:,:,rr).*dynamicsMultiply(Cluster_OronDist(:,:,rr)),'b')
    dotplot = Dynamic_OronDist(:,:,rr).*dynamicsMultiply(Cluster_OronDist(:,:,rr));
    
    for a=1:size(dotplot,1)
        for b =1:size(dotplot,2)
            if dotplot(a,b)==1
                plot(a,b,'.','Color','b','MarkerSize',2)
                hold on
            end
        end
    end

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([5,10,15,20])
    yticklabels({'4','8','12','16'})    
    
    title(my_rois_plus{rr})
    handles = colorbar;
    caxis([45 75])
    handles.TickDirection = 'out';
    handles.Box = 'off';
    handles.Label.String = '% Decoding accuracy';
    handles.Label.FontSize = 5;
    drawnow;
    colormap(colorsMap)

end
suplabel('Test time (s)','y')
suplabel('Train time (s)','x')
saveas(fig,[resultsSave,'/',nameAdd,'_Distractor_trainOrtestDist.svg'])    
close(fig)       


%% not cross temporal 

pValuesNoise = zeros(timeDim,Roi);

colors = cbrewer('qual','Set1',8);
colors = [colors(1,:);colors(2,:);colors(7,:)]

fig5 = figure(5)
set(fig5, 'PaperUnits', 'inches');
x_width=10 ;y_width=5;
set(fig5, 'PaperPosition', [0 0 x_width y_width]); %
time=1:timeDim;
for rr = 1:Roi_plot
    subplot(2,4,rr)
          
    stdDiag1 = std(cell2mat(cellfun(@(x) diag(x), wmNoiseDist_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);
    stdDiag2 = std(cell2mat(cellfun(@(x) diag(x), wmNoiseDistonMem_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);

    Diag1 = reshape(cell2mat(cellfun(@(x) diag(x), wmNoiseDist_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
    Diag2 = reshape(cell2mat(cellfun(@(x) diag(x), wmNoiseDistonMem_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
   
    pValuesNoise(:,rr) = permutationTest_timecourse(Diag1-Diag2);
    
    pvaluesCh = pValuesNoise(:,rr);
    
    pvaluesCh(pvaluesCh==0)=NaN;
    
    b1 = boundedline(time,diag(wmNoiseDist_reliable_MEAN{rr})-50',stdDiag1', ...
        time,diag(wmNoiseDistonMem_reliable_MEAN{rr})-50', stdDiag2', ...
        'cmap',colors, ...
        'alpha');
    hold on 
    
    plot(linspace(0.75,max(time)+0.25,max(time)*2), repelem(pvaluesCh,2)*-3,'k')

    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([0,10,20])
    
    text(22,-3, 'p <0.05','FontSize',5)

    title(my_rois_plus{rr})
    ylim([-5 25])
    xlim([0 timeDim])
    yline(0,'--')
    
end
suplabel('Time (s)','x')
suplabel('FCA above chance','y')
sgtitle('Noise distractor: trained on itself and no distractor trials')
lh=legend(b1);
legnames = {'Trained on itself','Trained on no-distractor trials'};
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
lpos = lh.Position;
lpos(2) = lpos(2) -0.05;
lpos(1) = lpos(1) +0.4;
lh.Position = lpos;
box off
saveas(fig5,[resultsSave,'/',nameAdd,'_NoiseDistractor_Diag.svg'])    
close(fig5)       


pValuesOr = zeros(timeDim,Roi);

fig6 = figure(6)
set(fig6, 'PaperUnits', 'inches');
x_width=10 ;y_width=5;
set(fig6, 'PaperPosition', [0 0 x_width y_width]); %
time=1:timeDim;
for rr = 1:Roi_plot
    subplot(2,4,rr)
          
    stdDiag1 = std(cell2mat(cellfun(@(x) diag(x), wmOrDist_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);
    stdDiag2 = std(cell2mat(cellfun(@(x) diag(x), wmOrDistonMem_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);

    Diag1 = reshape(cell2mat(cellfun(@(x) diag(x), wmOrDist_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
    Diag2 = reshape(cell2mat(cellfun(@(x) diag(x), wmOrDistonMem_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
   
    pValuesOr(:,rr) = permutationTest_timecourse(Diag1-Diag2);
    
    pvaluesCh = pValuesOr(:,rr);
    
    pvaluesCh(pvaluesCh==0)=NaN;
    
    b1 = boundedline(time,diag(wmOrDist_reliable_MEAN{rr})-50',stdDiag1', ...
        time,diag(wmOrDistonMem_reliable_MEAN{rr})-50', stdDiag2', ...
        'cmap',colors, ...
        'alpha');

    hold on;
    
    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([0,10,20])
    
    plot(linspace(0.75,max(time)+0.25,max(time)*2), repelem(pvaluesCh,2)*-3,'k')
    
    text(22,-3, 'p <0.05','FontSize',5)

    title(my_rois_plus{rr})
    ylim([-5 25])
    xlim([0 timeDim])
    yline(0,'--')
    
end
suplabel('Time (s)','x')
suplabel('FCA above chance','y')
sgtitle('Orientation distractor: trained on itself and no distractor trials')
lh=legend(b1);
legnames = {'Trained on itself','Trained on no-distractor trials'};
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
lpos = lh.Position;
lpos(2) = lpos(2) -0.05;
lpos(1) = lpos(1) +0.4;
lh.Position = lpos;
box off
saveas(fig6,[resultsSave,'/',nameAdd,'_OrDistractor_Diag.svg'])    
close(fig6)       


colors = cbrewer('div','PRGn',11);
colors = [colors(9,:);colors(11,:);colors(2,:)];

pValuesBetween = zeros(timeDim,Roi);

fig7 = figure(7)
set(fig7, 'PaperUnits', 'inches');
x_width=10 ;y_width=5;
set(fig7, 'PaperPosition', [0 0 x_width y_width]); %
time=1:timeDim;
for rr = 1:Roi_plot
    subplot(2,4,rr)
          
    stdDiag1 = std(cell2mat(cellfun(@(x) diag(x), wmOrDist_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);
    stdDiag2 = std(cell2mat(cellfun(@(x) diag(x), wmNoiseDist_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);
    stdDiag3 = std(cell2mat(cellfun(@(x) diag(x), wmNoDist_reliable(:,rr), 'UniformOutput', false)'),[],2)/sqrt(6);


    Diag1 = reshape(cell2mat(cellfun(@(x) diag(x), wmOrDist_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
    Diag2 = reshape(cell2mat(cellfun(@(x) diag(x), wmNoiseDist_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
    Diag3 = reshape(cell2mat(cellfun(@(x) diag(x), wmNoDist_reliable(:,rr), 'UniformOutput', false)),timeDim,[])';
   
    pValuesBetween(:,rr) = permutationTest_timecourse(Diag2-Diag1);
    pValuesBetween_No_Noise(:,rr) = permutationTest_timecourse(Diag3-Diag2);
    pValuesBetween_No_Or(:,rr) = permutationTest_timecourse(Diag3-Diag1);
    
    pvaluesCh = pValuesBetween(:,rr);
    pvaluesCh(pvaluesCh==0)=NaN;

    pvaluesCh2 = pValuesBetween_No_Noise(:,rr);
    pvaluesCh2(pvaluesCh2==0)=NaN;
    
    pvaluesCh3 = pValuesBetween_No_Or(:,rr);
    pvaluesCh3(pvaluesCh3==0)=NaN;

    b1 = boundedline(time,diag(wmOrDist_reliable_MEAN{rr})-50',stdDiag1', ...
        time,diag(wmNoiseDist_reliable_MEAN{rr})-50', stdDiag2', ...
        time,diag(wmNoDist_reliable_MEAN{rr})-50', stdDiag3', ...
        'cmap',colors, ...
        'alpha');

    hold on;
    
    plot(linspace(0.75,max(time)+0.25,max(time)*2), repelem(pvaluesCh,2)*-3,'k')
    plot(linspace(0.75,max(time)+0.25,max(time)*2), repelem(pvaluesCh2,2)*-2.5,'r')
    plot(linspace(0.75,max(time)+0.25,max(time)*2), repelem(pvaluesCh3,2)*-2,'b')
    
    text(22,-3, 'Noise vs Orientation p <0.05','FontSize',5)
    text(22,-2.5, 'No-distractor vs Orientation p <0.05','FontSize',5)
    text(22,-2, 'No-distractor vs Noise p <0.05','FontSize',5)
    
    xticks([5,10,15,20])
    xticklabels({'4','8','12','16'})
    yticks([0,10,20])
 
    title(my_rois_plus{rr})
    ylim([-5 25])
    xlim([0 timeDim])
    yline(0,'--')
    
end
suplabel('Time (s)','x')
suplabel('FCA above chance','y')
sgtitle('Orientation and Noise distractor trials trained on themselves')
lh=legend(b1);

legnames = {'Orientation','Noise'};
for i = 1:length(legnames)
    str{i} = ['\' sprintf('color[rgb]{%f,%f,%f} %s', colors(i, 1), colors(i, 2), colors(i, 3), legnames{i})];
end
lh.String = str;
lh.Box = 'off';
lpos = lh.Position;
lpos(2) = lpos(2) -0.05;
lpos(1) = lpos(1) +0.4;
lh.Position = lpos;
box off

saveas(fig7,[resultsSave,'/',nameAdd,'_OrNoise_Diag.svg'])    
close(fig7)       

end



