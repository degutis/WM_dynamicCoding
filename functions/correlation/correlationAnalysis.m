function correlationAnalysis(experiment,fwhm,increaseIPS)

global path2data

results_new = fullfile(path2data,'correlation_analysis_GLM',experiment,['fwhm_',num2str(fwhm),'_shuffle_0']);
results_new_Sh = fullfile(path2data,'correlation_analysis_GLM',experiment,['fwhm_',num2str(fwhm),'_shuffle_1']);

results_Mat = dir(fullfile(results_new,['sub*Analysis.mat']));
results_Mat_shuffle= dir(fullfile(results_new_Sh,['sub*Analysis.mat']));

sub_MM_all = [];
sub_PP_all = [];  

for sub = 1:length(results_Mat)    
    load(fullfile(results_Mat(sub).folder, results_Mat(sub).name))
    
    sub_MM_part = [];
    sub_PP_part = [];
    
    for roi = 1:size(sub_MM_full,2)
        MM = sub_MM_full(:,roi)';
        MM = cell2mat(MM);
        MM_mean = atanh(MM(:,1));
        
        PP = sub_PP_full(:,roi)';
        PP = cell2mat(PP); 
        PP_mean = atanh(PP(:,1));
            
        sub_MM_part = [sub_MM_part,{MM_mean}];
        sub_PP_part = [sub_PP_part,{PP_mean}];
    end
    
    sub_MM_all = [sub_MM_all;sub_MM_part];
    sub_PP_all = [sub_PP_all;sub_PP_part];
    
end


sub_MM_all_Sh = [];
sub_PP_all_Sh = []; 

for part_sh = 1:length(results_Mat_shuffle)    
    load(fullfile(results_Mat_shuffle(part_sh).folder, results_Mat_shuffle(part_sh).name))
    
    sub_MM_part_sh = [];
    sub_PP_part_sh = [];

    
    for roi_sh = 1:size(sub_MM_full,2)
        MM_sh = cell2mat(sub_MM_full(:,roi_sh)');
        MM_sh = atanh(MM_sh(:));

        PP_sh = cell2mat(sub_PP_full(:,roi_sh)');
        PP_sh = atanh(PP_sh(:));
        
        sub_MM_part_sh = [sub_MM_part_sh,{MM_sh}];
        sub_PP_part_sh = [sub_PP_part_sh,{PP_sh}];
    end
    
    sub_MM_all_Sh = [sub_MM_all_Sh;sub_MM_part_sh];
    sub_PP_all_Sh = [sub_PP_all_Sh;sub_PP_part_sh];    
end

numSubs = size(sub_MM_all,1);


if increaseIPS==1
    my_rois = {'V1','V2','V3','V3AB','V4','IPS','LO1','LO2'};
    
    IPS_MM_all = cell(numSubs,1);
    IPS_PP_all = cell(numSubs,1);
    
    IPS_MM_all_Sh = cell(numSubs,1);
    IPS_PP_all_Sh = cell(numSubs,1);
    
    
    for subb = 1:size(sub_MM_all,1)
        IPS_MM_all(subb,1) = {cell2mat(sub_MM_all(subb,6:9)')};
        IPS_PP_all(subb,1) = {cell2mat(sub_PP_all(subb,6:9)')};
       
        IPS_MM_all_Sh(subb,1) = {cell2mat(sub_MM_all_Sh(subb,6:9)')};
        IPS_PP_all_Sh(subb,1) = {cell2mat(sub_PP_all_Sh(subb,6:9)')};
    end
    
    sub_MM_all = [sub_MM_all(:,1:5),IPS_MM_all,sub_MM_all(:,10:11)];
    sub_PP_all = [sub_PP_all(:,1:5),IPS_PP_all,sub_PP_all(:,10:11)];
    
    sub_MM_all_Sh = [sub_MM_all_Sh(:,1:5),IPS_MM_all_Sh,sub_MM_all_Sh(:,10:11)];
    sub_PP_all_Sh = [sub_PP_all_Sh(:,1:5),IPS_PP_all_Sh,sub_PP_all_Sh(:,10:11)];       
end

if increaseIPS==1
    results_new = fullfile([results_new,'_IPS']);
end

for ROI2 = 1:length(my_rois)
    
    results_ROI = fullfile(results_new,my_rois{ROI2}); 
    mkdir(results_ROI)

    numSubs = size(sub_MM_all,1);
           
    fig1=figure(1)
    for partPlot = 1:size(sub_PP_all,1)
        subplot(ceil(numSubs/2),ceil(size(sub_PP_all,1)/ceil(numSubs/2)),partPlot)
        histogram(tanh(sub_PP_all{partPlot,ROI2}),[-1:0.1:1],'Normalization','Probability')
        hold on
        h2=histogram(tanh(sub_PP_all_Sh{partPlot,ROI2}),[-1:0.1:1],'Normalization','Probability')
        set(h2,'facecolor','red')
        set(h2,'facealpha',0.5)
        box off
        xlabel('Pearson Correlation')
        ylabel('Proportion (voxels)')
        box off
        ylim([0 0.25])
        title(['Subject: ',num2str(partPlot)]);       
    end
    sgtitle(['Perception-perception correlation. ROI:',my_rois{ROI2}])
    saveas(fig1,fullfile(results_ROI,['Across_sub_PPcorr_',my_rois{ROI2},'.svg']));
    close(fig1)


    fig2=figure(2)
    for partPlot = 1:size(sub_MM_all,1)
        subplot(ceil(numSubs/2),ceil(size(sub_MM_all,1)/ceil(numSubs/2)),partPlot)
        histogram(tanh(sub_MM_all{partPlot,ROI2}),[-1:0.1:1],'Normalization','Probability')
        hold on
        h2=histogram(tanh(sub_MM_all_Sh{partPlot,ROI2}),[-1:0.1:1],'Normalization','Probability')
        set(h2,'facecolor','red')
        set(h2,'facealpha',0.5)
        box off
        xlabel('Pearson Correlation')
        ylabel('Proportion (voxels)')
        box off
        ylim([0 0.25])
        title(['Subject: ',num2str(partPlot)]);       
    end
    sgtitle(['Memory-memory correlation. ROI:',my_rois{ROI2}])
    saveas(fig2,fullfile(results_ROI,['Across_sub_MMcorr_',my_rois{ROI2},'.svg']));
    close(fig2)
    
end

%% Select all voxels above 75th percentile of the shuffled distribution 
sub_MM_all_percentile = cellfun(@(x) prctile(x,75), sub_MM_all_Sh, 'UniformOutput', false);
sub_PP_all_percentile = cellfun(@(x) prctile(x,75), sub_PP_all_Sh, 'UniformOutput', false);

sub_MM_all_voxSelect = cell(size(sub_MM_all,1),size(sub_MM_all,2));
sub_PP_all_voxSelect = cell(size(sub_PP_all,1),size(sub_PP_all,2));

for subj = 1:size(sub_MM_all,1)
    for ROI = 1:size(sub_MM_all,2)        
        sub_MM_all_voxSelect(subj,ROI) =  {sub_MM_all{subj,ROI}>sub_MM_all_percentile{subj,ROI}};
        sub_PP_all_voxSelect(subj,ROI) = {sub_PP_all{subj,ROI}>sub_PP_all_percentile{subj,ROI}};
    end
end

sub_voxSelect = cell(size(sub_MM_all_voxSelect,1),size(sub_MM_all_voxSelect,2));

for aa = 1:size(sub_MM_all_voxSelect,1)
    for bb=1:size(sub_MM_all_voxSelect,2)
        sub_voxSelect(aa,bb) = {sub_MM_all_voxSelect{aa,bb}==1 & sub_PP_all_voxSelect{aa,bb}==1};
    end
end
        
sub_MM_all_voxSelect_prop = cellfun(@(x) sum(x)/length(x), sub_MM_all_voxSelect, 'UniformOutput', false);
sub_PP_all_voxSelect_prop = cellfun(@(x) sum(x)/length(x), sub_PP_all_voxSelect, 'UniformOutput', false);

sub_voxSelect_prop = cellfun(@(x) sum(x)/length(x), sub_voxSelect, 'UniformOutput', false);


fig7=figure(7)
boxplot(cell2mat(sub_voxSelect_prop),'Labels',my_rois)
ylabel(['Proportion vox meeting 75th percentile corr cirteria in MM and PP'])
box off
sgtitle(['Proportion of voxels across ROIs'])
saveas(fig7,fullfile(results_new,['VoxROI_prop.svg']));
close(fig7)

fig8=figure(8)
try
    boxplot(cell2mat(sub_voxSelect_prop)','Labels',{'Sub1','Sub2','Sub3','Sub4','Sub5','Sub6'})
catch
    boxplot(cell2mat(sub_voxSelect_prop)','Labels',{'Sub1','Sub2','Sub3','Sub4','Sub5','Sub6','Sub7'})
end
ylabel(['Proportion vox meeting 75th percentile corr cirteria in MM and PP'])
box off
sgtitle(['Proportion of voxels across subjects'])
saveas(fig8,fullfile(results_new,['VoxSub_prop.svg']));
close(fig8)


save(fullfile(results_new,'FullResults_correlationAnalysis.mat'),'sub_MM_all','sub_MM_all_Sh','sub_MM_all_voxSelect','sub_MM_all_voxSelect_prop','sub_PP_all','sub_PP_all_Sh','sub_PP_all_voxSelect','sub_PP_all_voxSelect_prop','sub_voxSelect','sub_voxSelect_prop')


%% Plot selected voxels

% Select the voxels 

sub_voxSelect_MM = cell(size(sub_MM_all_voxSelect,1),size(sub_MM_all_voxSelect,2));
sub_voxSelect_PP = cell(size(sub_MM_all_voxSelect,1),size(sub_MM_all_voxSelect,2));

for cc=1:size(sub_voxSelect,1)
    for dd = 1:size(sub_voxSelect,2)
        sub_voxSelect_MM(cc,dd) = {sub_MM_all{cc,dd}(sub_voxSelect{cc,dd},:)};
        sub_voxSelect_PP(cc,dd) = {sub_PP_all{cc,dd}(sub_voxSelect{cc,dd},:)};
    end
end

my_x = linspace(1,length(my_rois), length(my_rois));
my_xs = [my_x-.2; my_x+.2]';
colors = cbrewer('qual', 'Set1', 3);
fig16=figure(16); 
set(fig16, 'PaperUnits', 'inches');
x_width=10 ;y_width=6;
set(fig16, 'PaperPosition', [0 0 x_width y_width]); %
subplot(2,2,[1 2])
for condition = 1:2
    if condition==1
        my_means = mean(cell2mat(cellfun(@(x) mean(tanh(x)), sub_voxSelect_MM, 'UniformOutput', false)))';
        my_ses = std(cell2mat(cellfun(@(x) mean(tanh(x)), sub_voxSelect_MM, 'UniformOutput', false)))'/sqrt(subj);
        h_1 = bar(my_xs(:,condition)',my_means,.4);
        hold on
        errorbar(my_xs(:,condition),my_means,my_ses,'k.');
        set(h_1, 'FaceColor',  colors(condition,:), 'EdgeColor', 'none');
        set(gca,'TickDir','out','XTick',my_x,'XTickLabel',...
        my_rois);
    elseif condition==2
        my_means = mean(cell2mat(cellfun(@(x) mean(tanh(x)), sub_voxSelect_PP, 'UniformOutput', false)))';
        my_ses = std(cell2mat(cellfun(@(x) mean(tanh(x)), sub_voxSelect_PP, 'UniformOutput', false)))'/sqrt(subj);
        h_2 = bar(my_xs(:,condition)',my_means,.4);      
        hold on
        errorbar(my_xs(:,condition),my_means,my_ses,'k.');
        set(h_2, 'FaceColor',  colors(condition,:), 'EdgeColor', 'none');
        set(gca,'TickDir','out','XTick',my_x,'XTickLabel',...
        my_rois);
    end    
        
end
legend([h_1,h_2],'Memory','Perception')
legend boxoff
ylabel('Average correlation','Fontsize',10);
xlabel('ROIs');
ylim([0.5 1])
box off

subplot(2,2,3)

my_means = mean(cell2mat(cellfun(@(x) length(x), sub_voxSelect_MM, 'UniformOutput', false)))';
my_ses = std(cell2mat(cellfun(@(x) length(x), sub_voxSelect_MM, 'UniformOutput', false)))'/sqrt(subj);
h_1 = bar(my_xs(:,condition)',my_means,.4);
hold on
errorbar(my_xs(:,condition),my_means,my_ses,'k.');
set(h_1, 'FaceColor',  colors(3,:), 'EdgeColor', 'none');
set(gca,'TickDir','out','XTick',my_x,'XTickLabel',...
my_rois);
ylim([0 800])
ylabel('Number of voxels','Fontsize',10);
xlabel('ROIs');
box off

subplot(2,2,4)

my_means = mean(cell2mat(sub_voxSelect_prop))';
my_ses = std(cell2mat(sub_voxSelect_prop))'/sqrt(subj);
h_1 = bar(my_xs(:,condition)',my_means,.4);
hold on
errorbar(my_xs(:,condition),my_means,my_ses,'k.');
set(h_1, 'FaceColor',  colors(3,:), 'EdgeColor', 'none');
set(gca,'TickDir','out','XTick',my_x,'XTickLabel',...
my_rois);
ylim([0 0.4])
ylabel('Proportion of total voxels','Fontsize',10);
xlabel('ROIs');
box off

saveas(fig16,fullfile(results_new,['MeanMM_PP_reliable_voxels.svg']));
close(fig16)



end



