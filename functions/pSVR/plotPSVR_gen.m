function plotPSVR_gen(experiment,expandIPS)

global path2data

nameAdd = 'ReliableVox';

resultsSave = fullfile(path2data,'pSVR_analysis',experiment,'temporal',nameAdd);
dirReliable = dir(fullfile(resultsSave,'gen_runwise*mat'));

mkdir(resultsSave)

if expandIPS==1
    my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS','LO1','LO2'};
else
    my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
end

if strcmp('E1',experiment)
    N=6;
elseif strcmp('E2',experiment)
    N=7;
end

Roi=numel(my_rois_plus);

for f = 1:numel(dirReliable)
    subNow = load(fullfile(dirReliable(f).folder,dirReliable(f).name));
    
    wmNoDist_reliable(f,:) = cellfun(@(x) mean(x(:,:),2), subNow.wmNoonNo_fcaMatrix(f,:), 'UniformOutput', false);
    wmOrDist_reliable(f,:) = cellfun(@(x) mean(x(:,:),2), subNow.wmOronOr_fcaMatrix(f,:), 'UniformOutput', false);
    wmNoiseDist_reliable(f,:) = cellfun(@(x) mean(x(:,:),2), subNow.wmNoiseonNoise_fcaMatrix(f,:), 'UniformOutput', false);          
    
    wmOrDistonMem_reliable(f,:) = cellfun(@(x) mean(x(:,:),2), subNow.wmOronMem_fcaMatrix(f,:), 'UniformOutput', false);          
    wmNoiseDistonMem_reliable(f,:) = cellfun(@(x) mean(x(:,:),2), subNow.wmNoiseonMem_fcaMatrix(f,:), 'UniformOutput', false);          
end

wmNoDist_reliable = cell2mat(wmNoDist_reliable);
wmOrDist_reliable = cell2mat(wmOrDist_reliable);
wmNoiseDist_reliable = cell2mat(wmNoiseDist_reliable);

wmOrDistonMem_reliable = cell2mat(wmOrDistonMem_reliable);
wmNoiseDistonMem_reliable = cell2mat(wmNoiseDistonMem_reliable);

wmOrDist_proportion = wmOrDist_reliable./wmNoDist_reliable;
wmNoiseDist_proportion = wmNoiseDist_reliable./wmNoDist_reliable;
wmOrDistonMem_proportion = wmOrDistonMem_reliable./wmNoDist_reliable;
wmNoiseDistonMem_proportion = wmNoiseDistonMem_reliable./wmNoDist_reliable;

% wmOrDist_proportion = wmOrDist_reliable./wmOrDist_reliable;
% wmNoiseDist_proportion = wmNoiseDist_reliable./wmNoiseDist_reliable;
% wmOrDistonMem_proportion = wmOrDistonMem_reliable./wmOrDist_reliable;
% wmNoiseDistonMem_proportion = wmNoiseDistonMem_reliable./wmNoiseDist_reliable;



stats_combo = [wmOrDist_proportion-wmOrDistonMem_proportion, wmNoiseDist_proportion-wmNoiseDistonMem_proportion];

pval = runStatsPermutationAgainstNull(stats_combo,[],[],1,4,1444);
pval_noise = pval(length(my_rois_plus)+1:length(my_rois_plus)*2);
pval_or = pval(1:length(my_rois_plus));


baseline_combo = [(wmNoiseDist_proportion-1)*-1, (wmNoiseDistonMem_proportion-1)*-1, (wmOrDist_proportion-1)*-1, (wmOrDistonMem_proportion-1)*-1];
pval_baseline = runStatsPermutationAgainstNull(baseline_combo,[],[],1,4,140);


pval_noise_noise_baseline = pval_baseline(1:8);
pval_noise_nodist_baseline = pval_baseline(9:16);

pval_or_or_baseline = pval_baseline(17:24);
pval_or_nodist_baseline = pval_baseline(25:32);


pvalTable_mod = [pval_noise_nodist_baseline,pval_noise_noise_baseline,pval_or_nodist_baseline,pval_or_or_baseline,pval_noise,pval_or];
pvalTable_mod = pvalTable_mod';
RowNames = ["Noise-between baseline", "Noise-within baseline", "Orientation-between baseline", "Orientation-within baseline", "Noise generalization", "Orientation generalization"];

tableFinal = table(RowNames',pvalTable_mod(:,1),pvalTable_mod(:,2),pvalTable_mod(:,3),pvalTable_mod(:,4),pvalTable_mod(:,5),pvalTable_mod(:,6),pvalTable_mod(:,7),pvalTable_mod(:,8));
tableFinal.Properties.VariableNames = ["Test",'V1','V2','V3','V3AB','V4','IPS','LO1','LO2'];
writetable(tableFinal,fullfile(resultsSave,"pvalTable_genCode.xlsx"),'Sheet',1)

fig1 = figure(1)
colors = cbrewer('div','PRGn',11);
colors = [colors(10,:);colors(10,:);colors(8,:);colors(8,:)];
set(fig1, 'PaperUnits', 'inches');
x_width=20 ;y_width=7;
set(fig1, 'PaperPosition', [0 0 x_width y_width]); %
for r=1:length(my_rois_plus)
    subplot(2,4,r)
%     violinplot([wmNoiseDist_proportion(:,r), wmNoiseDistonMem_proportion(:,r),wmOrDist_proportion(:,r),wmOrDistonMem_proportion(:,r)],[],'ViolinColor',colors)
    violinplot([wmNoiseDistonMem_proportion(:,r),wmNoiseDist_proportion(:,r),wmOrDistonMem_proportion(:,r),wmOrDist_proportion(:,r)],[],'ViolinColor',colors)
    for sub=1:N
%         plot([1,2],[wmNoiseDist_proportion(sub,r),wmNoiseDistonMem_proportion(sub,r)],"k")
%         plot([3,4],[wmOrDist_proportion(sub,r),wmOrDistonMem_proportion(sub,r)],"k")
        plot([1,2],[wmNoiseDistonMem_proportion(sub,r),wmNoiseDist_proportion(sub,r)],"k")
        plot([3,4],[wmOrDistonMem_proportion(sub,r),wmOrDist_proportion(sub,r)],"k")
    end
    
    xticklabels({'Noise-between','Noise-within','Orientation-between','Orientation-within'})
    xtickangle(45)
    yline(1,"--")
    title(my_rois_plus{r})
    mysigstar(gca,1.5,0.81,pval_noise(r))
    mysigstar(gca,3.5,0.81,pval_or(r))
    
    mysigstar(gca,1,1.1,pval_noise_noise_baseline(r))
    mysigstar(gca,2,1.1,pval_noise_nodist_baseline(r))
    mysigstar(gca,3,1.1,pval_or_or_baseline(r))
    mysigstar(gca,4,1.1,pval_or_nodist_baseline(r))   
    
    ylim([0.7,1.1])
%     yline(1)
end
saveas(fig1,fullfile(resultsSave,[nameAdd,'_generalizationViolin_prop.svg']))    
close(fig1)       

