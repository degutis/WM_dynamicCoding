function GLMCorr(experiment, shuffleOr,numBootstrapIter,runSmoothing,fwhm)

global path2data

my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS0','IPS1','IPS2','IPS3','LO1','LO2'};
numROIs = length(my_rois_plus);
results_new = fullfile(path2data,'correlation_analysis_GLM',experiment,['fwhm_',num2str(fwhm),'_shuffle_',num2str(shuffleOr)]);

mkdir(results_new)
mkdir(fullfile(results_new,'individSub'))

sampleFiles = dir(fullfile('../raw_data',experiment,'SampleFile_*'));
timeCourseFiles = dir(fullfile('../raw_data',experiment,'TimeCourses_*'));

memory_1_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);
memory_2_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);
perception_1_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);
perception_2_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);

memory_1_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';
memory_2_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';
perception_1_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';
perception_2_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';

memory_1_SNR_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);
memory_2_SNR_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);
perception_1_SNR_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);
perception_2_SNR_full = cell2table(cell(length(sampleFiles)*numBootstrapIter,numROIs), 'VariableNames', my_rois_plus);

memory_1_SNR_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';
memory_2_SNR_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';
perception_1_SNR_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';
perception_2_SNR_full.Iteration = repelem([1:numBootstrapIter],length(sampleFiles))';

rng(13011991+shuffleOr)


for sub = 1:length(sampleFiles)    
    load(fullfile(sampleFiles(sub).folder,sampleFiles(sub).name));
    load(fullfile(timeCourseFiles(sub).folder,timeCourseFiles(sub).name));

    %% Working memory
    
    wmTrialIndex = InfoMain(:,3)==0;
    wmTrials = InfoMain(InfoMain(:,3)==0,2);
    
    wmResponsesROI = cell2table(AvgResponsesMain);
    wmResponsesROI.Properties.VariableNames = my_rois;
    wmResponsesNoDist = varfun(@(x) x{:}(wmTrialIndex,:), wmResponsesROI);
    wmResponsesNoDist.Properties.VariableNames = my_rois;
   
    %% Localizer 
    
    if experiment=='E1'
        locTrialIndex_full = InfoLocalizer;
        locResponsesROI = cell2table(AvgResponsesLocalizer);

    elseif experiment == 'E2'
        locTrialIndex_full = InfoLocalizer2;
        locResponsesROI = cell2table(AvgResponsesLocalizer2);

    end
    
    locIndexDonut = locTrialIndex_full(:,2);
    locIndexDonut = locIndexDonut(locTrialIndex_full(:,1)==1);
    
    locResponsesROI.Properties.VariableNames = my_rois;
    locResponsesROI = varfun(@(x) x{:}(:,:), locResponsesROI);
    locResponsesROI.Properties.VariableNames = my_rois;
    
    %% Compute bins for actual estimation of tuning curve bins in degrees
    
    edges = [0:30:180]; % 6 orientation bins

    wmTrials_binned = wmTrials;
    for bin = 1:length(edges)-1
        wmTrials_binned(wmTrials_binned>=edges(bin)+1 & wmTrials_binned<=edges(bin+1)) = (edges(bin)+edges(bin+1))/2;
    end

    locTrials_binned = locIndexDonut;
    for bin = 1:length(edges)-1
        locTrials_binned(locTrials_binned>=edges(bin)+1 & locTrials_binned<=edges(bin+1)) = (edges(bin)+edges(bin+1))/2;
    end

    %% Run correlation analysis for n iterations and on all ROIs
    
    sub_MM_full = {};
    sub_PP_full = {};

    for iterB = 1:numBootstrapIter
 
        if shuffleOr == 1
            locTrials_binned = locTrials_binned(randperm(length(locTrials_binned)));
            wmTrials_binned = wmTrials_binned(randperm(length(wmTrials_binned)));
        end
        
        disp(['Sub: ',num2str(sub),' Iteration: ',num2str(iterB)])        
        
        binnedOr = unique(wmTrials_binned); 

        selectionMemory = zeros(length(wmTrials_binned),2,numBootstrapIter);   
        
        %% Working memory
           
        wm_c1 = cvpartition(wmTrials_binned(1:length(wmTrials_binned)/3),'KFold',2,'Stratify',true);
        ind_wm_train1 = training(wm_c1,1);
        ind_wm_test1 = test(wm_c1,1); 

        wm_c2 = cvpartition(wmTrials_binned(length(wmTrials_binned)/3+1:length(wmTrials_binned)/3*2),'KFold',2,'Stratify',true);
        ind_wm_train2 = training(wm_c2,1);
        ind_wm_test2 = test(wm_c2,1); 

        wm_c3 = cvpartition(wmTrials_binned(length(wmTrials_binned)/3*2+1:length(wmTrials_binned)),'KFold',2,'Stratify',true);
        ind_wm_train3 = training(wm_c3,1);
        ind_wm_test3 = test(wm_c3,1); 

        ind_wm_train = [ind_wm_train1;ind_wm_train2;ind_wm_train3];
        ind_wm_test = [ind_wm_test1;ind_wm_test2;ind_wm_test3];
           
        selectionMemory(:,:,iterB) = [ind_wm_train, ind_wm_test];
        
        if runSmoothing ==1
            [labels,ixlabels] = sort([wmTrials(ind_wm_train),wmTrials(ind_wm_test)],'ascend');
            
            wmResponsesNoDist_smoothed_train = table();
            wmResponsesNoDist_smoothed_test = table();
            
            for column = 1:length(my_rois)
                yt = wmResponsesNoDist.(column);
                            
                yt_1 = yt(ind_wm_train,:);    
                yt_1 = yt_1(ixlabels(:,1),:);
                
                yt_2 = yt(ind_wm_test,:);
                yt_2 = yt_2(ixlabels(:,2),:);

                wmResponsesNoDist_smoothed_train(:,column) = table(feature_space_smoothing(yt_1,labels(:,1),fwhm)); 
                wmResponsesNoDist_smoothed_test(:,column) = table(feature_space_smoothing(yt_2,labels(:,2),fwhm)); 
            end
            
            trialInput1 = wmTrials_binned(ind_wm_train);
            trialInput2 = wmTrials_binned(ind_wm_test);
            
            trialInput1 = trialInput1(ixlabels(:,1));
            trialInput2 = trialInput2(ixlabels(:,2));
                        
            [memory_1, memory_1_SNR] = estimateTuningCurveGLM(wmResponsesNoDist_smoothed_train,trialInput1,binnedOr,1);
            [memory_2, memory_2_SNR] = estimateTuningCurveGLM(wmResponsesNoDist_smoothed_test,trialInput2,binnedOr,1);
        
        else    
            
            [memory_1, memory_1_SNR] = estimateTuningCurveGLM(wmResponsesNoDist(ind_wm_train,:),wmTrials_binned(ind_wm_train),binnedOr,1);
            [memory_2, memory_2_SNR] = estimateTuningCurveGLM(wmResponsesNoDist(ind_wm_test,:),wmTrials_binned(ind_wm_test),binnedOr,1);
        end
        
    
        %% Perception
        
        locTrials_binned_select = locTrials_binned;
        
        if length(locTrials_binned_select)==300 || length(locTrials_binned_select)==216 || length(locTrials_binned_select)==342 || length(locTrials_binned_select)==270 || length(locTrials_binned_select)==180 || length(locTrials_binned_select)==234
            loc_c21 = cvpartition(locTrials_binned_select(1:length(locTrials_binned_select)/3),'KFold',2,'Stratify',true);
            ind_loc_train1 = training(loc_c21,1);
            ind_loc_test1 = test(loc_c21,1); 
            
            loc_c22 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/3+1:length(locTrials_binned_select)/3*2),'KFold',2,'Stratify',true);
            ind_loc_train2 = training(loc_c22,1);
            ind_loc_test2 = test(loc_c22,1); 
    
            loc_c23 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/3*2+1:length(locTrials_binned_select)),'KFold',2,'Stratify',true);
            ind_loc_train3 = training(loc_c23,1);
            ind_loc_test3 = test(loc_c23,1); 
            
            ind_loc_train = [ind_loc_train1;ind_loc_train2;ind_loc_train3];
            ind_loc_test = [ind_loc_test1;ind_loc_test2;ind_loc_test3];
        
        elseif length(locTrials_binned_select)==320 
            loc_c21 = cvpartition(locTrials_binned_select(1:length(locTrials_binned_select)/4),'KFold',2,'Stratify',true);
            ind_loc_train1 = training(loc_c21,1);
            ind_loc_test1 = test(loc_c21,1); 
            
            loc_c22 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/4+1:length(locTrials_binned_select)/2),'KFold',2,'Stratify',true);
            ind_loc_train2 = training(loc_c22,1);
            ind_loc_test2 = test(loc_c22,1); 
    
            loc_c23 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/2+1:length(locTrials_binned_select)/4*3),'KFold',2,'Stratify',true);
            ind_loc_train3 = training(loc_c23,1);
            ind_loc_test3 = test(loc_c23,1); 
            
            loc_c24 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/4*3+1:length(locTrials_binned_select)),'KFold',2,'Stratify',true);
            ind_loc_train4 = training(loc_c24,1);
            ind_loc_test4 = test(loc_c24,1); 
            
            ind_loc_train = [ind_loc_train1;ind_loc_train2;ind_loc_train3;ind_loc_train4];
            ind_loc_test = [ind_loc_test1;ind_loc_test2;ind_loc_test3;ind_loc_test4];
        
        elseif length(locTrials_binned_select)==340 || length(locTrials_binned_select)==360

            loc_c21 = cvpartition(locTrials_binned_select(1:length(locTrials_binned_select)/5),'KFold',2,'Stratify',true);
            ind_loc_train1 = training(loc_c21,1);
            ind_loc_test1 = test(loc_c21,1); 
            
            loc_c22 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/5+1:length(locTrials_binned_select)/5*2),'KFold',2,'Stratify',true);
            ind_loc_train2 = training(loc_c22,1);
            ind_loc_test2 = test(loc_c22,1); 
    
            loc_c23 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/5*2+1:length(locTrials_binned_select)/5*3),'KFold',2,'Stratify',true);
            ind_loc_train3 = training(loc_c23,1);
            ind_loc_test3 = test(loc_c23,1); 
            
            loc_c24 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/5*3+1:length(locTrials_binned_select)/5*4),'KFold',2,'Stratify',true);
            ind_loc_train4 = training(loc_c24,1);
            ind_loc_test4 = test(loc_c24,1); 

            loc_c25 = cvpartition(locTrials_binned_select(length(locTrials_binned_select)/5*4+1:length(locTrials_binned_select)),'KFold',2,'Stratify',true);
            ind_loc_train5 = training(loc_c25,1);
            ind_loc_test5 = test(loc_c25,1); 
            
            ind_loc_train = [ind_loc_train1;ind_loc_train2;ind_loc_train3;ind_loc_train4;ind_loc_train5];
            ind_loc_test = [ind_loc_test1;ind_loc_test2;ind_loc_test3;ind_loc_test4;ind_loc_test5];
        end           
       
        locResponsesROIs_select = locResponsesROI;
        
        
        if runSmoothing ==1
            [labels_loc,ixlabels_loc] = sort([locIndexDonut(ind_loc_train),locIndexDonut(ind_loc_test)],'ascend');
            
            locResponses_smoothed_train = table();
            locResponses_smoothed_test = table();
            
            for column = 1:length(my_rois)
                yt_loc = locResponsesROIs_select.(column);
            
                yt_loc_1 = yt_loc(ind_loc_train,:);
                yt_loc_1 = yt_loc_1(ixlabels_loc(:,1),:);
                
                yt_loc_2 = yt_loc(ind_loc_test,:);
                yt_loc_2 = yt_loc_2(ixlabels_loc(:,2),:);

                locResponses_smoothed_train(:,column) = table(feature_space_smoothing(yt_loc_1,labels_loc(:,1),fwhm)); 
                locResponses_smoothed_test(:,column) = table(feature_space_smoothing(yt_loc_2,labels_loc(:,2),fwhm)); 
            end
            
            trialInput1_loc = locTrials_binned(ind_loc_train);
            trialInput2_loc = locTrials_binned(ind_loc_test);
            
            trialInput1_loc = trialInput1_loc(ixlabels_loc(:,1));
            trialInput2_loc = trialInput2_loc(ixlabels_loc(:,2));
            
            [perception_1,perception_1_SNR] = estimateTuningCurveGLM(locResponses_smoothed_train,trialInput1_loc,binnedOr,1);
            [perception_2,perception_2_SNR] = estimateTuningCurveGLM(locResponses_smoothed_test,trialInput2_loc,binnedOr,1);
        
        else    
            [perception_1,perception_1_SNR] = estimateTuningCurveGLM(locResponsesROIs_select(ind_loc_train,:),locTrials_binned_select(ind_loc_train),binnedOr,1);
            [perception_2,perception_2_SNR] = estimateTuningCurveGLM(locResponsesROIs_select(ind_loc_test,:),locTrials_binned_select(ind_loc_test),binnedOr,1);
        end
 
        
        %% Correlation 
        
        memory_1_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = memory_1;
        memory_2_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = memory_2; 
        perception_1_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = perception_1;
        perception_2_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = perception_2;
        
        memory_1_SNR_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = memory_1_SNR;
        memory_2_SNR_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = memory_2_SNR; 
        perception_1_SNR_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = perception_1_SNR;
        perception_2_SNR_full(((iterB-1)*length(sampleFiles))+sub,1:end-1) = perception_2_SNR;
        
        
        memory_1 = cell2table(memory_1);
        memory_1.Properties.VariableNames = my_rois;
        
        memory_2 = cell2table(memory_2);
        memory_2.Properties.VariableNames = my_rois;

        perception_1 = cell2table(perception_1);
        perception_1.Properties.VariableNames = my_rois;

        perception_2 = cell2table(perception_2);
        perception_2.Properties.VariableNames = my_rois;
        
        sub_MM = {};
        sub_PP = {};
        sub_MP = {};
        
        for roiColumn = 1:size(locResponsesROI,2)
            roi_MM = [];
            roi_PP = [];

            for voxel = 1:size(memory_1.(roiColumn){:},1)
                corr_MM = corrcoef(memory_1.(roiColumn){:}(voxel,:),memory_2.(roiColumn){:}(voxel,:));
                corr_MM = corr_MM(2);
                
                corr_PP = corrcoef(perception_1.(roiColumn){:}(voxel,:),perception_2.(roiColumn){:}(voxel,:));
                corr_PP = corr_PP(2);
                            
                roi_MM = [roi_MM;corr_MM];
                roi_PP = [roi_PP;corr_PP];
            end
            sub_MM = [sub_MM,roi_MM];
            sub_PP = [sub_PP,roi_PP];
        end
 
        sub_MM_full(iterB,:) = sub_MM;
        sub_PP_full(iterB,:) = sub_PP;
    end

    save(fullfile(results_new,['sub',num2str(sub),'_correlationAnalysis.mat']),'wmTrials_binned','locTrials_binned','selectionMemory','sub_MM_full','sub_PP_full','my_rois','memory_1_SNR','memory_2_SNR','perception_1_SNR','perception_2_SNR')
    
    for roiColumn2 = 1:size(locResponsesROI,2)
        
        MM = sub_MM_full(:,roiColumn2)';
        MM = cell2mat(MM);
        MM_mean = mean(MM,2);
        MM_mean_acVox = mean(MM_mean);
        MM_std_acVox = std(MM_mean); 
        
        PP = sub_PP_full(:,roiColumn2)';
        PP = cell2mat(PP); 
        PP_mean = mean(PP,2);
        PP_mean_acVox = mean(PP_mean);
        PP_std_acVox = std(PP_mean);                
        
        fig1 = figure(1)
        scatter(MM_mean,PP_mean)
        xlim([-1 1])
        ylim([-1 1])
        xlabel('Memory-memory correlation')
        ylabel('Perception-perception correlation')
        box off
        axis equal
        title(['Subject:',num2str(sub),' MM PP voxels ',my_rois{roiColumn2}]);       
        saveas(fig1,fullfile(results_new,'individSub',['sub',num2str(sub),'CorrMatrix_MemPerc_',my_rois{roiColumn2},'.png']));
        close(fig1)
        
        fig2 = figure(2)
        bar([MM_mean_acVox, PP_mean_acVox])
        hold on
        errorbar([MM_mean_acVox, PP_mean_acVox],[MM_std_acVox, PP_std_acVox],'LineStyle','none')
        set(gca,'XTickLabel',[{'MM'},{'PP'}])
        xlabel('Correlations')
        ylabel('Correlation coefficient')
        box off
        title(['Subject:',num2str(sub),' Correlation coef across voxels ',my_rois{roiColumn2}]);       
        saveas(fig2,fullfile(results_new,'individSub',['sub',num2str(sub),'CorrAcVox',my_rois{roiColumn2},'.png']));
        close(fig2)
                
    end
   
end    

end
