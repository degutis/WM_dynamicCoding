function subspace_grandPCA_projection_validated(experiment,EVCcombine,numRois_dist,saveLargeFiles)

rng(111)

voxSelect=1;
IPSconcat = 1;
numBins=4;  
numIter=1000;

nameHoldOut = '50_50';
holdOut = 1/2;

global path2data

if voxSelect==0
    nameAdd = 'AllVoxels';
elseif voxSelect==1
    nameAdd = 'ReliableVox';
end

sampleFiles = dir(fullfile('../raw_data',experiment,'SampleFile_*'));
timeCourseFiles = dir(fullfile('../raw_data',experiment,'TimeCourses_*'));

if EVCcombine==1
    nameAdd = [nameAdd,'_bins',num2str(numBins),'_Projection_VALIDATED',nameHoldOut,'_EVCcombo']; %Added new name
else
    nameAdd = [nameAdd,'_bins',num2str(numBins),'_Projection_VALIDATED',nameHoldOut]; %Added new name
end

my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS','LO1','LO2'};
indexVox = load(fullfile(path2data,'correlation_analysis_GLM',experiment,'fwhm_30_shuffle_0_IPS','FullResults_correlationAnalysis.mat'),'sub_voxSelect');
nameAdd = [nameAdd,'_acrossAllVoxels'];

if EVCcombine==1
    my_rois_plus = {'V1-V3AB'};
end

resultsSave = fullfile(path2data,'PCA_SubspaceAnalysis',experiment,nameAdd);
mkdir(resultsSave)


trs_to_recon = [-1 1:22];  % 1 tr before stim onset, 22 tr's after
num_trs_to_recon = length(trs_to_recon);
num_trials = 108*3;

currentWM_noDist_cell = cell(numBins,length(my_rois_plus),2,numIter); 
currentWM_orDist_cell = cell(numBins,length(my_rois_plus),2,numIter); 
currentWM_distractor_cell = cell(numBins,length(my_rois_plus),2,numIter); 
currentWM_noiseDist_cell = cell(numBins,length(my_rois_plus),2,numIter); 

delayPeriodSelected = 6:23;
lengthDelay = length(delayPeriodSelected);

if isfile(fullfile(resultsSave,'splitHalf_distractor.mat'))==0

    for sub = 1:length(sampleFiles)
    
        % load data and such
        load(fullfile(sampleFiles(sub).folder,sampleFiles(sub).name));
        load(fullfile(timeCourseFiles(sub).folder,timeCourseFiles(sub).name));   
    
        if strcmp('E2',experiment)
            TimeCoursesLocalizer = TimeCoursesLocalizer2;
            AvgResponsesLocalizer = AvgResponsesLocalizer2;
        end
        
        if IPSconcat==1
            TimeCoursesLocalizer = [TimeCoursesLocalizer(1:5),{cat(2,TimeCoursesLocalizer{6:9})},TimeCoursesLocalizer(10:11)];
            TimeCoursesMain = [TimeCoursesMain(1:5),{cat(2,TimeCoursesMain{6:9})},TimeCoursesMain(10:11)];
            
            AvgResponsesLocalizer = [AvgResponsesLocalizer(1:5),{cat(2,AvgResponsesLocalizer{6:9})},AvgResponsesLocalizer(10:11)];
            AvgResponsesMain = [AvgResponsesMain(1:5),{cat(2,AvgResponsesMain{6:9})},AvgResponsesMain(10:11)];
        end
        
        wmTrialIndex = InfoMain(:,3)==0;
        wmTrials = InfoMain(InfoMain(:,3)==0,2);
        
        % E1 trials        
        if strcmp('E1',experiment)
            wmTrialIndex_noise = isnan(InfoMain(:,3));
            wmTrialIndex_orDist = InfoMain(:,3)>0;
    
            wmTrials_noise = InfoMain(isnan(InfoMain(:,3)),2);
            wmTrials_orDist = InfoMain(InfoMain(:,3)>0,2);
        
        % E2 trials
        elseif strcmp('E2',experiment)
            wmTrialIndex_noise = InfoMain(:,3)<1 & InfoMain(:,3)>0;
            wmTrialIndex_orDist = InfoMain(:,3)>=1;
            
            wmTrials_noise = InfoMain(InfoMain(:,3)<1 & InfoMain(:,3)>0,2);
            wmTrials_orDist = InfoMain(InfoMain(:,3)>=1,2);
        end
    
        distractor_WM = InfoMain(wmTrialIndex_orDist,3);
        
        if voxSelect==1
            for roiInit = 1:numel(TimeCoursesMain)            
                TimeCoursesMain(roiInit)={TimeCoursesMain{roiInit}(:,indexVox.sub_voxSelect{sub,roiInit})};
            end
        end
        
        edges = [0:45:180]; % 4 orientation bins

        wmTrials_binned = wmTrials;
        wmTrials_binned_orientation = wmTrials_orDist;
        wmTrials_binned_noise = wmTrials_noise;
        distractor_binned_orientation = distractor_WM;        

        for bin = 1:length(edges)-1
            wmTrials_binned(wmTrials_binned>=edges(bin)+1 & wmTrials_binned<=edges(bin+1)) = (edges(bin)+edges(bin+1))/2;
            wmTrials_binned_orientation(wmTrials_binned_orientation>=edges(bin)+1 & wmTrials_binned_orientation<=edges(bin+1)) = (edges(bin)+edges(bin+1))/2;
            wmTrials_binned_noise(wmTrials_binned_noise>=edges(bin)+1 & wmTrials_binned_noise<=edges(bin+1)) = (edges(bin)+edges(bin+1))/2;
            distractor_binned_orientation(distractor_binned_orientation>=edges(bin)+1 & distractor_binned_orientation<=edges(bin+1)) = (edges(bin)+edges(bin+1))/2;            
        end
        
        binsVal = unique(wmTrials_binned); 
       
        for roi = 1:numel(my_rois_plus)
            
            disp(['Subject: ',num2str(sub),' Roi: ',num2str(roi)]);
                
            if EVCcombine==1
                if roi==1
                    my_TC = cat(2,TimeCoursesMain{1:4});
                end
            else
                my_TC = TimeCoursesMain{roi};
            end
    
            % Make a matrix of num_trs_to_recon by num_voxels by num_trials
            rsp_wm = NaN(num_trs_to_recon,size(my_TC,2),num_trials);
            for trial = 1:num_trials
                rsp_wm(:,:,trial) = my_TC(TrialStartInd_main(trial)+trs_to_recon(1):TrialStartInd_main(trial)+trs_to_recon(end)-1,:);
            end
                
            for b=1:numIter
                seed = ceil((130191+roi+sub)*b);
                currentWM_noDist_cell(sub,roi,:,b) = extractTrialsMean_splitHalf(delayPeriodSelected,wmTrialIndex,rsp_wm,wmTrials_binned,binsVal,'HoldOut',holdOut,seed);
                currentWM_orDist_cell(sub,roi,:,b) = extractTrialsMean_splitHalf(delayPeriodSelected,wmTrialIndex_orDist,rsp_wm,wmTrials_binned_orientation,binsVal,'HoldOut',holdOut,seed);
                currentWM_distractor_cell(sub,roi,:,b) = extractTrialsMean_splitHalf(delayPeriodSelected,wmTrialIndex_orDist,rsp_wm,distractor_binned_orientation,binsVal,'HoldOut',holdOut,seed);
                currentWM_noiseDist_cell(sub,roi,:,b) = extractTrialsMean_splitHalf(delayPeriodSelected,wmTrialIndex_noise,rsp_wm,wmTrials_binned_noise,binsVal,'HoldOut',holdOut,seed);        
            end
        end
    end
    
    if saveLargeFiles ==1
        save(fullfile(resultsSave,'splitHalf_distractor.mat'),'currentWM_distractor_cell','-v7.3')
        save(fullfile(resultsSave,'splitHalf_no.mat'),'currentWM_noDist_cell','-v7.3')
        save(fullfile(resultsSave,'splitHalf_or.mat'),'currentWM_orDist_cell','-v7.3')
        save(fullfile(resultsSave,'splitHalf_noise.mat'),'currentWM_noiseDist_cell','-v7.3')
    end
else
    binsVal = [22.5;67.5;112.5;157.5]; 
    load(fullfile(resultsSave,'splitHalf_distractor.mat'))
    load(fullfile(resultsSave,'splitHalf_no.mat'))
    load(fullfile(resultsSave,'splitHalf_or.mat'))
    load(fullfile(resultsSave,'splitHalf_noise.mat'))
end


%%
if isfile(fullfile(resultsSave,'angles.mat'))==0
    
    timePointCell = {[1:3],[4:6]; ...
        [1:3],[7:9]; ...
        [1:3],[10:12];...
        [1:3],[13:15];...
        [1:3],[16:18];...
        [4:6],[7:9];...
        [4:6],[10:12];...
        [4:6],[13:15];...
        [4:6],[16:18];...
        [7:9],[10:12];...
        [7:9],[13:15];...
        [7:9],[16:18];...
        [10:12],[13:15];...
        [10:12],[16:18];...
        [13:15],[16:18]};
    
    timePointIndex = {[1,2],[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6],[3,4],[3,5],[3,6],[4,5],[4,6],[5,6]};
    
    %% Angle calculation across time
    
    NoDistAngles = zeros(size(currentWM_noDist_cell,1),size(currentWM_noDist_cell,2),length(timePointCell),numIter);
    OrDistAngles = zeros(size(currentWM_noDist_cell,1),size(currentWM_noDist_cell,2),length(timePointCell),numIter);
    NoiseDistAngles = zeros(size(currentWM_noDist_cell,1),size(currentWM_noDist_cell,2),length(timePointCell),numIter);
    
    numSubs=size(currentWM_noDist_cell,1);
    
    for sub=1:size(currentWM_noDist_cell,1)
        for roi=1:size(currentWM_noDist_cell,2)
            disp(['sub: ',num2str(sub),' roi: ',num2str(roi)])
            for b=1:numIter
                [NoDistAngles(sub,roi,:,b)]=runPCAangleSplit(currentWM_noDist_cell(:,:,:,b),sub,roi,timePointCell,timePointIndex);
                [OrDistAngles(sub,roi,:,b)]=runPCAangleSplit(currentWM_orDist_cell(:,:,:,b),sub,roi,timePointCell,timePointIndex);
                [NoiseDistAngles(sub,roi,:,b)]=runPCAangleSplit(currentWM_noiseDist_cell(:,:,:,b),sub,roi,timePointCell,timePointIndex);
            end
        end
    end
    
    NoDistAngles = mean(NoDistAngles,4);
    OrDistAngles = mean(OrDistAngles,4);
    NoiseDistAngles = mean(NoiseDistAngles,4);
    
    comboAngles = cat(3,NoDistAngles,OrDistAngles,NoiseDistAngles);
    pval_combo = runStatsPermutationAgainstNull(comboAngles,[],timePointIndex,1,3,14011);
       
    pval = pval_combo(:,:,:,1);
    pval_Or = pval_combo(:,:,:,2);
    pval_Noise = pval_combo(:,:,:,3);

    % create and save p-val table

    pvalTable = zeros(length(timePointIndex),3)
    for type = 1:3
        for a = 1:length(timePointIndex)
            pvalTable(a,type) = pval_combo(timePointIndex{a}(1),timePointIndex{a}(2),1,type)
        end
    end

    pvalTable_mod = [pvalTable(:,1),pvalTable(:,3),pvalTable(:,2)];
    pvalTable = array2table(pvalTable_mod);
    pvalTable.Properties.VariableNames = ["No distractor","Noise distractor","Orientation distractor"];
    RowNames = ["4.8 - 7.2s","4.8 - 9.6s","4.8 - 12s","4.8 - 14.4s", "4.8 - 16.8s",...
        "7.2 - 9.6s","7.2 - 12s","7.2 - 14.4s","7.2 - 16.8s","9.6 - 12s","9.6 - 14.4s","9.6 - 16.8s","12 - 14.4",...
        "12 - 16.8s" , "14.4s - 16.8"];

    tableFinal = table(RowNames',pvalTable_mod(:,1),pvalTable_mod(:,2),pvalTable_mod(:,3));
    tableFinal.Properties.VariableNames = ["Time points","No distractor","Noise distractor","Orientation distractor"];

    writetable(tableFinal,fullfile(resultsSave,"pvalTable.xlsx"),'Sheet',1)

    meanNoDistAngles = mean(NoDistAngles);
    meanOrDistAngles = mean(OrDistAngles);
    meanNoiseDistAngles = mean(NoiseDistAngles);

    NoDistAngle = nan(size(pval));
    OrDistAngle = nan(size(pval));
    NoiseDistAngle = nan(size(pval));

    for t3 = 1:length(timePointIndex)  
        NoDistAngle(timePointIndex{t3}(1),timePointIndex{t3}(2),roi) = meanNoDistAngles(:,roi,t3);
        OrDistAngle(timePointIndex{t3}(1),timePointIndex{t3}(2),roi) = meanOrDistAngles(:,roi,t3);
        NoiseDistAngle(timePointIndex{t3}(1),timePointIndex{t3}(2),roi) = meanNoiseDistAngles(:,roi,t3);
    end

    save(fullfile(resultsSave,'angles.mat'),'NoDistAngles','OrDistAngles','NoiseDistAngles','timePointCell','timePointIndex','pval','pval_Or','pval_Noise')
 
    colors = cbrewer('div','PRGn',11);
    colors = [colors(3,:);colors(8,:);colors(10,:)];

    for rr=1:1
        fig1 = figure(1)    
        set(fig1, 'PaperUnits', 'inches');
        x_width=4 ;y_width=4;
        set(fig1, 'PaperPosition', [0 0 x_width y_width]); %
    
        for t=1:5
            h1 = ploterr(t-0.1, squeeze(mean(NoDistAngles(:,rr,t),1)), [], squeeze(std(NoDistAngles(:,rr,t),[],1))/sqrt(numSubs), 'k.', 'abshhxy', 0);
            set(h1,'Color',colors(1,:))
            hold on 
            mysigstar(gca,t-0.1,13,pval(timePointIndex{t}(1),timePointIndex{t}(2),rr))
            hold on
            h2 = ploterr(t, squeeze(mean(OrDistAngles(:,rr,t),1)), [], squeeze(std(OrDistAngles(:,rr,t),[],1))/sqrt(numSubs), 'k.', 'abshhxy', 0);
            set(h2,'Color',colors(2,:))
            hold on     
            mysigstar(gca,t,14,pval_Or(timePointIndex{t}(1),timePointIndex{t}(2),rr))
            hold on
            h3 = ploterr(t+0.1, squeeze(mean(NoiseDistAngles(:,rr,t),1)), [], squeeze(std(NoiseDistAngles(:,rr,t),[],1))/sqrt(numSubs), 'k.', 'abshhxy', 0);
            set(h3,'Color',colors(3,:)) 
            hold on     
            mysigstar(gca,t,15,pval_Noise(timePointIndex{t}(1),timePointIndex{t}(2),rr))
        end
        box off
        ylim([-5 30])
        xlim([0.8 6.3])
        xticklabels({'4.8s','7.2s','9.6s','12s','14.4s','16.8s'})
        suptitle(my_rois_plus{rr})
        saveas(fig1,fullfile(resultsSave,['ROIsubspaceAngle_plot_timepoint_',my_rois_plus{rr},'.svg']))
        close(fig1)
    end
        
    plotSubspaceTriangleMultipleROI(pval,pval_Or,pval_Noise,cat(4,NoDistAngle,OrDistAngle,NoiseDistAngle),timePointIndex,my_rois_plus,0.05,resultsSave,'ttest_Dynamic.svg');
end

%% Projection distractor

if isfile(fullfile(resultsSave,'pvalStableDistractor.mat'))==0 

    if isempty(currentWM_orDist_cell{1})
        load(fullfile(resultsSave,'splitHalf_distractor.mat'))
        load(fullfile(resultsSave,'splitHalf_or.mat'))
    end
    
    [pval_WMspace,OrDistValues,perc_explained,perc_explained_distractor]=projectDistractor(numRois_dist,my_rois_plus,resultsSave,currentWM_orDist_cell,currentWM_distractor_cell);

    save(fullfile(resultsSave,'pvalStableDistractor.mat'),'pval_WMspace','OrDistValues','perc_explained','perc_explained_distractor')
end


%% Projection of datapoints onto subsapces across time

timePointCell_averageProject = {[1:3],[4:6],[7:9],[10:12],[13:14],[16:18]};
timePointCell = {[4:6],[10:12],[16:18]};   
timePointName = {'EarlyDelay','MiddleDelay','LateDelay'};

perc_explained_no = zeros(numBins,length(timePointCell),length(my_rois_plus));
perc_explained_or = zeros(numBins,length(timePointCell),length(my_rois_plus));
perc_explained_noise = zeros(numBins,length(timePointCell),length(my_rois_plus));

for tPoint=1:length(timePointCell)

    roi_wm_projection_onto_Stable = zeros(numBins,2,length(timePointCell_averageProject),length(my_rois_plus));
    roi_wm_projection_onto_Stable_orientation = zeros(numBins,2,length(timePointCell_averageProject),length(my_rois_plus));
    roi_wm_projection_onto_Stable_noise = zeros(numBins,2,length(timePointCell_averageProject),length(my_rois_plus));
    
    for rr = 1:length(my_rois_plus)
        currentROI_train = cat(2,currentWM_noDist_cell{:,rr,1,1});
        currentROI_test = cat(2,currentWM_noDist_cell{:,rr,2,1});

        currentROI_stable_train = mean(currentROI_train(:,:,timePointCell{tPoint}),3);
        currentROI_stable_train = currentROI_stable_train - mean(currentROI_stable_train,1);

        [~,W,~,~,perc_explained_no(:,tPoint,rr)]=pca(currentROI_stable_train','NumComponents',2);

        for t=1:length(timePointCell_averageProject)
            timePointProject = currentROI_test(:,:,timePointCell_averageProject{t});
            timePointProject = mean(timePointProject,3);
            timePointProject = timePointProject - mean(timePointProject,1);           
            roi_wm_projection_onto_Stable(:,:,t,rr) = timePointProject*W;
        end

        % Orientation

        currentROI_or_train = cat(2,currentWM_orDist_cell{:,rr,1,1});
        currentROI_or_test = cat(2,currentWM_orDist_cell{:,rr,2,1});

        currentROI_stable_or_train = mean(currentROI_or_train(:,:,timePointCell{tPoint}),3);
        currentROI_stable_or_train = currentROI_stable_or_train - mean(currentROI_stable_or_train,1);

        [~,W_or,~,~,perc_explained_or(:,tPoint,rr)]=pca(currentROI_stable_or_train','NumComponents',2);

        for t=1:length(timePointCell_averageProject)
            timePointProject = currentROI_or_test(:,:,timePointCell_averageProject{t});
            timePointProject = mean(timePointProject,3);
            timePointProject = timePointProject - mean(timePointProject,1);           
            roi_wm_projection_onto_Stable_orientation(:,:,t,rr) = timePointProject*W_or;
        end    

        % Noise

        currentROI_noise_train = cat(2,currentWM_noiseDist_cell{:,rr,1,1});
        currentROI_noise_test = cat(2,currentWM_noiseDist_cell{:,rr,2,1});

        currentROI_stable_noise_train = mean(currentROI_noise_train(:,:,timePointCell{tPoint}),3);
        currentROI_stable_noise_train = currentROI_stable_noise_train - mean(currentROI_stable_noise_train,1);

        [~,W_noise,~,~,perc_explained_noise(:,tPoint,rr)]=pca(currentROI_stable_noise_train','NumComponents',2);

        for t=1:length(timePointCell_averageProject)
            timePointProject = currentROI_noise_test(:,:,timePointCell_averageProject{t});
            timePointProject = mean(timePointProject,3);
            timePointProject = timePointProject - mean(timePointProject,1);
            roi_wm_projection_onto_Stable_noise(:,:,t,rr) = timePointProject*W_noise;
        end           
        
    end

    %% Standardize for plotting
    
    for roi=1:length(my_rois_plus)
    
        currentROI_no = roi_wm_projection_onto_Stable(:,:,:,roi);
        maxTwoDim_no = max(max(abs(currentROI_no),[],3));
        
        currentROI_or = roi_wm_projection_onto_Stable_orientation(:,:,:,roi);
        maxTwoDim_or = max(max(abs(currentROI_or),[],3));

        currentROI_noise = roi_wm_projection_onto_Stable_noise(:,:,:,roi);
        maxTwoDim_noise = max(max(abs(currentROI_noise),[],3));
        
        for dim=1:2
            for t=1:length(timePointCell_averageProject)
                currentROI_no(:,dim,t) = currentROI_no(:,dim,t)/maxTwoDim_no(dim);
                currentROI_or(:,dim,t) = currentROI_or(:,dim,t)/maxTwoDim_or(dim);
                currentROI_noise(:,dim,t) = currentROI_noise(:,dim,t)/maxTwoDim_noise(dim);
            end
        end
        
        currentROI_no = currentROI_no*10;
        currentROI_or = currentROI_or*10;
        currentROI_noise = currentROI_noise*10;
        
        roi_wm_projection_onto_Stable(:,:,:,roi)=currentROI_no;
        roi_wm_projection_onto_Stable_orientation(:,:,:,roi)=currentROI_or;
        roi_wm_projection_onto_Stable_noise(:,:,:,roi)=currentROI_noise;
        
    end
        
    %%  PCA time

    colors = cbrewer('qual', 'Set1', 6);

    if EVCcombine==0
        l=4;
    else
        l=1;
    end
    fig2 = figure(2)    
    for a= 1:l
    
        set(fig2, 'PaperUnits', 'inches');
        x_width=10 ;y_width=10;
        set(fig2, 'PaperPosition', [0 0 x_width y_width]); %

        subplot(4,3,3*(a-1)+1)

        p12=plot(squeeze(roi_wm_projection_onto_Stable(1,1,:,a)),squeeze(roi_wm_projection_onto_Stable(1,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p1=scatter(squeeze(roi_wm_projection_onto_Stable(1,1,1:size(roi_wm_projection_onto_Stable,3)/2,a)),squeeze(roi_wm_projection_onto_Stable(1,2,1:size(roi_wm_projection_onto_Stable,3)/2,a)),5,'filled','HandleVisibility','off')
        p1.MarkerFaceColor = colors(1,:);
        p1.MarkerEdgeColor = colors(1,:);
        p1.MarkerFaceAlpha = 0.5;
        p1.LineWidth = 1;
        hold on   
        p11=scatter(squeeze(roi_wm_projection_onto_Stable(1,1,size(roi_wm_projection_onto_Stable,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable(1,2,size(roi_wm_projection_onto_Stable,3)/2:end,a)),5,'filled')
        p11.MarkerFaceColor = colors(1,:);
        p11.MarkerEdgeColor = colors(1,:);
        p11.LineWidth = 1;
        
        hold on   

        p22=plot(squeeze(roi_wm_projection_onto_Stable(2,1,:,a)),squeeze(roi_wm_projection_onto_Stable(2,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p2=scatter(squeeze(roi_wm_projection_onto_Stable(2,1,1:size(roi_wm_projection_onto_Stable,3)/2,a)),squeeze(roi_wm_projection_onto_Stable(2,2,1:size(roi_wm_projection_onto_Stable,3)/2,a)),5,'filled','HandleVisibility','off')
        p2.MarkerFaceColor = colors(2,:);
        p2.MarkerEdgeColor = colors(2,:);
        p2.MarkerFaceAlpha = 0.5;
        p2.LineWidth = 1;
        
        hold on   
        p21=scatter(squeeze(roi_wm_projection_onto_Stable(2,1,size(roi_wm_projection_onto_Stable,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable(2,2,size(roi_wm_projection_onto_Stable,3)/2:end,a)),5,'filled')
        p21.MarkerFaceColor = colors(2,:);
        p21.MarkerEdgeColor = colors(2,:);
        p21.LineWidth = 1;
        
        hold on   

        p32=plot(squeeze(roi_wm_projection_onto_Stable(3,1,:,a)),squeeze(roi_wm_projection_onto_Stable(3,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p3=scatter(squeeze(roi_wm_projection_onto_Stable(3,1,1:size(roi_wm_projection_onto_Stable,3)/2,a)),squeeze(roi_wm_projection_onto_Stable(3,2,1:size(roi_wm_projection_onto_Stable,3)/2,a)),5,'filled','HandleVisibility','off')
        p3.MarkerFaceColor = colors(3,:);
        p3.MarkerEdgeColor = colors(3,:);
        p3.MarkerFaceAlpha = 0.5;
        p3.LineWidth = 1;
        
        hold on   
        p31=scatter(squeeze(roi_wm_projection_onto_Stable(3,1,size(roi_wm_projection_onto_Stable,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable(3,2,size(roi_wm_projection_onto_Stable,3)/2:end,a)),5,'filled')
        p31.MarkerFaceColor = colors(3,:);
        p31.MarkerEdgeColor = colors(3,:);
        p31.LineWidth = 1;
        
        hold on   

        p42=plot(squeeze(roi_wm_projection_onto_Stable(4,1,:,a)),squeeze(roi_wm_projection_onto_Stable(4,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p4=scatter(squeeze(roi_wm_projection_onto_Stable(4,1,1:size(roi_wm_projection_onto_Stable,3)/2,a)),squeeze(roi_wm_projection_onto_Stable(4,2,1:size(roi_wm_projection_onto_Stable,3)/2,a)),5,'filled','HandleVisibility','off')
        p4.MarkerFaceColor = colors(4,:);
        p4.MarkerEdgeColor = colors(4,:);
        p4.MarkerFaceAlpha = 0.5;
        p4.LineWidth = 1;
        hold on   
        p41=scatter(squeeze(roi_wm_projection_onto_Stable(4,1,size(roi_wm_projection_onto_Stable,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable(4,2,size(roi_wm_projection_onto_Stable,3)/2:end,a)),5,'filled')
        p41.MarkerFaceColor = colors(4,:);
        p41.MarkerEdgeColor = colors(4,:);
        p41.LineWidth = 1;
        

        if numBins==6
            hold on   

            p52=plot(squeeze(roi_wm_projection_onto_Stable(5,1,:,a)),squeeze(roi_wm_projection_onto_Stable(5,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
            hold on
            p5=scatter(squeeze(roi_wm_projection_onto_Stable(5,1,1:size(roi_wm_projection_onto_Stable,3)/2,a)),squeeze(roi_wm_projection_onto_Stable(5,2,1:size(roi_wm_projection_onto_Stable,3)/2,a)),5,'filled','HandleVisibility','off')
            p5.MarkerFaceColor = colors(5,:);
            p5.MarkerEdgeColor = colors(5,:);
            p5.MarkerFaceAlpha = 0.5;
            p5.LineWidth = 1;
            
            hold on   
            p51=scatter(squeeze(roi_wm_projection_onto_Stable(5,1,size(roi_wm_projection_onto_Stable,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable(5,2,size(roi_wm_projection_onto_Stable,3)/2:end,a)),5,'filled')
            p51.MarkerFaceColor = colors(5,:);
            p51.MarkerEdgeColor = colors(5,:);
            p51.LineWidth = 1;
            
            hold on   

            p62=plot(squeeze(roi_wm_projection_onto_Stable(6,1,:,a)),squeeze(roi_wm_projection_onto_Stable(6,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
            hold on
            p6=scatter(squeeze(roi_wm_projection_onto_Stable(6,1,1:size(roi_wm_projection_onto_Stable,3)/2,a)),squeeze(roi_wm_projection_onto_Stable(6,2,1:size(roi_wm_projection_onto_Stable,3)/2,a)),5,'filled','HandleVisibility','off')
            p6.MarkerFaceColor = colors(6,:);
            p6.MarkerEdgeColor = colors(6,:);
            p6.MarkerFaceAlpha = 0.5;
            p6.LineWidth = 1;
            
            hold on   
            p61=scatter(squeeze(roi_wm_projection_onto_Stable(6,1,size(roi_wm_projection_onto_Stable,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable(6,2,size(roi_wm_projection_onto_Stable,3)/2:end,a)),5,'filled')
            p61.MarkerFaceColor = colors(6,:);
            p61.MarkerEdgeColor = colors(6,:);
            p61.LineWidth = 1;
        end

        axis square
        hold on   
        box off
        xlim([-10,10])
        ylim([-10,10])
        xticks([-10 0 10])
        yticks([-10 0 10])
        xlabel(strjoin([string(round(perc_explained_no(1,tPoint,rr))),"%"],""))
        ylabel(strjoin([string(round(perc_explained_no(2,tPoint,rr))),"%"],""))

        
        
        subplot(4,3,3*(a-1)+2)

        p12=plot(squeeze(roi_wm_projection_onto_Stable_orientation(1,1,:,a)),squeeze(roi_wm_projection_onto_Stable_orientation(1,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p1=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(1,1,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_orientation(1,2,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),5,'filled','HandleVisibility','off')
        p1.MarkerFaceColor = colors(1,:);
        p1.MarkerEdgeColor = colors(1,:);
        p1.MarkerFaceAlpha = 0.5;
        p1.LineWidth = 1;
        
        hold on   
        p11=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(1,1,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_orientation(1,2,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),5,'filled')
        p11.MarkerFaceColor = colors(1,:);
        p11.MarkerEdgeColor = colors(1,:);
        p11.LineWidth = 1;
        
        hold on   

        p22=plot(squeeze(roi_wm_projection_onto_Stable_orientation(2,1,:,a)),squeeze(roi_wm_projection_onto_Stable_orientation(2,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p2=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(2,1,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_orientation(2,2,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),5,'filled','HandleVisibility','off')
        p2.MarkerFaceColor = colors(2,:);
        p2.MarkerEdgeColor = colors(2,:);
        p2.MarkerFaceAlpha = 0.5;
        p2.LineWidth = 1;
        
        hold on   
        p21=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(2,1,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_orientation(2,2,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),5,'filled')
        p21.MarkerFaceColor = colors(2,:);
        p21.MarkerEdgeColor = colors(2,:);
        p21.LineWidth = 1;
        
        hold on   

        p32=plot(squeeze(roi_wm_projection_onto_Stable_orientation(3,1,:,a)),squeeze(roi_wm_projection_onto_Stable_orientation(3,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p3=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(3,1,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_orientation(3,2,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),5,'filled','HandleVisibility','off')
        p3.MarkerFaceColor = colors(3,:);
        p3.MarkerEdgeColor = colors(3,:);
        p3.MarkerFaceAlpha = 0.5;
        p3.LineWidth = 1;
        
        hold on   
        p31=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(3,1,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_orientation(3,2,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),5,'filled')
        p31.MarkerFaceColor = colors(3,:);
        p31.MarkerEdgeColor = colors(3,:);
        p31.LineWidth = 1;
        
        hold on   

        p42=plot(squeeze(roi_wm_projection_onto_Stable_orientation(4,1,:,a)),squeeze(roi_wm_projection_onto_Stable_orientation(4,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p4=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(4,1,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_orientation(4,2,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),5,'filled','HandleVisibility','off')
        p4.MarkerFaceColor = colors(4,:);
        p4.MarkerEdgeColor = colors(4,:);
        p4.MarkerFaceAlpha = 0.5;
        p4.LineWidth = 1;
        
        hold on   
        p41=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(4,1,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_orientation(4,2,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),5,'filled')
        p41.MarkerFaceColor = colors(4,:);
        p41.MarkerEdgeColor = colors(4,:);
        p41.LineWidth = 1;
        

        if numBins==6
            hold on   

            p52=plot(squeeze(roi_wm_projection_onto_Stable_orientation(5,1,:,a)),squeeze(roi_wm_projection_onto_Stable_orientation(5,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
            hold on
            p5=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(5,1,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_orientation(5,2,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),5,'filled','HandleVisibility','off')
            p5.MarkerFaceColor = colors(5,:);
            p5.MarkerEdgeColor = colors(5,:);
            p5.MarkerFaceAlpha = 0.5;
            p5.LineWidth = 1;
            
            hold on   
            p51=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(5,1,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_orientation(5,2,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),5,'filled')
            p51.MarkerFaceColor = colors(5,:);
            p51.MarkerEdgeColor = colors(5,:);
            p51.LineWidth = 1;
            
            hold on   

            p62=plot(squeeze(roi_wm_projection_onto_Stable_orientation(6,1,:,a)),squeeze(roi_wm_projection_onto_Stable_orientation(6,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
            hold on
            p6=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(6,1,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_orientation(6,2,1:size(roi_wm_projection_onto_Stable_orientation,3)/2,a)),5,'filled','HandleVisibility','off')
            p6.MarkerFaceColor = colors(6,:);
            p6.MarkerEdgeColor = colors(6,:);
            p6.MarkerFaceAlpha = 0.5;
            hold on   
            p61=scatter(squeeze(roi_wm_projection_onto_Stable_orientation(6,1,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_orientation(6,2,size(roi_wm_projection_onto_Stable_orientation,3)/2:end,a)),5,'filled')
            p61.MarkerFaceColor = colors(6,:);
            p61.MarkerEdgeColor = colors(6,:);
            p61.LineWidth = 1;
        end

        axis square
        hold on
        box off
        xlim([-10,10])
        ylim([-10,10])  
        xticks([-10 0 10])
        yticks([-10 0 10])
        xlabel(strjoin([string(round(perc_explained_or(1,tPoint,rr))),"%"],""))
        ylabel(strjoin([string(round(perc_explained_or(2,tPoint,rr))),"%"],""))
        

        subplot(4,3,3*(a-1)+3)

        p12=plot(squeeze(roi_wm_projection_onto_Stable_noise(1,1,:,a)),squeeze(roi_wm_projection_onto_Stable_noise(1,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p1=scatter(squeeze(roi_wm_projection_onto_Stable_noise(1,1,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_noise(1,2,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),5,'filled','HandleVisibility','off')
        p1.MarkerFaceColor = colors(1,:);
        p1.MarkerEdgeColor = colors(1,:);
        p1.MarkerFaceAlpha = 0.5;
        hold on   
        p11=scatter(squeeze(roi_wm_projection_onto_Stable_noise(1,1,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_noise(1,2,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),5,'filled')
        p11.MarkerFaceColor = colors(1,:);
        p11.MarkerEdgeColor = colors(1,:);
        p11.LineWidth = 1;
        
        hold on   

        p22=plot(squeeze(roi_wm_projection_onto_Stable_noise(2,1,:,a)),squeeze(roi_wm_projection_onto_Stable_noise(2,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p2=scatter(squeeze(roi_wm_projection_onto_Stable_noise(2,1,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_noise(2,2,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),5,'filled','HandleVisibility','off')
        p2.MarkerFaceColor = colors(2,:);
        p2.MarkerEdgeColor = colors(2,:);
        p2.MarkerFaceAlpha = 0.5;
        p2.LineWidth = 1;
        
        hold on   
        p21=scatter(squeeze(roi_wm_projection_onto_Stable_noise(2,1,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_noise(2,2,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),5,'filled')
        p21.MarkerFaceColor = colors(2,:);
        p21.MarkerEdgeColor = colors(2,:);
        p21.LineWidth = 1;
        
        hold on   

        p32=plot(squeeze(roi_wm_projection_onto_Stable_noise(3,1,:,a)),squeeze(roi_wm_projection_onto_Stable_noise(3,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p3=scatter(squeeze(roi_wm_projection_onto_Stable_noise(3,1,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_noise(3,2,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),5,'filled','HandleVisibility','off')
        p3.MarkerFaceColor = colors(3,:);
        p3.MarkerEdgeColor = colors(3,:);
        p3.MarkerFaceAlpha = 0.5;
        p3.LineWidth = 1;
        
        hold on   
        p31=scatter(squeeze(roi_wm_projection_onto_Stable_noise(3,1,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_noise(3,2,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),5,'filled')
        p31.MarkerFaceColor = colors(3,:);
        p31.MarkerEdgeColor = colors(3,:);
        p31.LineWidth = 1;
        
        hold on   

        p42=plot(squeeze(roi_wm_projection_onto_Stable_noise(4,1,:,a)),squeeze(roi_wm_projection_onto_Stable_noise(4,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
        hold on
        p4=scatter(squeeze(roi_wm_projection_onto_Stable_noise(4,1,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_noise(4,2,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),5,'filled','HandleVisibility','off')
        p4.MarkerFaceColor = colors(4,:);
        p4.MarkerEdgeColor = colors(4,:);
        p4.MarkerFaceAlpha = 0.5;
        p4.LineWidth = 1;
        
        hold on   
        p41=scatter(squeeze(roi_wm_projection_onto_Stable_noise(4,1,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_noise(4,2,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),5,'filled')
        p41.MarkerFaceColor = colors(4,:);
        p41.MarkerEdgeColor = colors(4,:);
        p41.LineWidth = 1;
        
        if numBins==6
            hold on   

            p52=plot(squeeze(roi_wm_projection_onto_Stable_noise(5,1,:,a)),squeeze(roi_wm_projection_onto_Stable_noise(5,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
            hold on
            p5=scatter(squeeze(roi_wm_projection_onto_Stable_noise(5,1,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_noise(5,2,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),5,'filled','HandleVisibility','off')
            p5.MarkerFaceColor = colors(5,:);
            p5.MarkerEdgeColor = colors(5,:);
            p5.MarkerFaceAlpha = 0.5;
            p5.LineWidth = 1;
            
            hold on   
            p51=scatter(squeeze(roi_wm_projection_onto_Stable_noise(5,1,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_noise(5,2,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),5,'filled')
            p51.MarkerFaceColor = colors(5,:);
            p51.MarkerEdgeColor = colors(5,:);
            p51.LineWidth = 1;
           hold on   

            p62=plot(squeeze(roi_wm_projection_onto_Stable_noise(6,1,:,a)),squeeze(roi_wm_projection_onto_Stable_noise(6,2,:,a)),"LineWidth",1,"Color",[0,0,0,0.2],'HandleVisibility','off')
            hold on
            p6=scatter(squeeze(roi_wm_projection_onto_Stable_noise(6,1,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),squeeze(roi_wm_projection_onto_Stable_noise(6,2,1:size(roi_wm_projection_onto_Stable_noise,3)/2,a)),5,'filled','HandleVisibility','off')
            p6.MarkerFaceColor = colors(6,:);
            p6.MarkerEdgeColor = colors(6,:);
            p6.MarkerFaceAlpha = 0.5;
            p6.LineWidth = 1;
            
            hold on   
            p61=scatter(squeeze(roi_wm_projection_onto_Stable_noise(6,1,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),squeeze(roi_wm_projection_onto_Stable_noise(6,2,size(roi_wm_projection_onto_Stable_noise,3)/2:end,a)),5,'filled')
            p61.MarkerFaceColor = colors(6,:);
            p61.MarkerEdgeColor = colors(6,:);
            p61.LineWidth = 1;
        end

        box off
        lgnd{1} = sprintf(num2str(binsVal(1)));
        lgnd{2} = sprintf(num2str(binsVal(2)));
        lgnd{3} = sprintf(num2str(binsVal(3)))
        lgnd{4} = sprintf(num2str(binsVal(4)))

        if numBins==6
            lgnd{5} = sprintf(num2str(binsVal(5)))
            lgnd{6} = sprintf(num2str(binsVal(6)))
        end

        axis square
        xlim([-10,10])
        ylim([-10,10])        
        xticks([-10 0 10])
        yticks([-10 0 10])
        xlabel(strjoin([string(round(perc_explained_noise(1,tPoint,rr))),"%"],""))
        ylabel(strjoin([string(round(perc_explained_noise(2,tPoint,rr))),"%"],""))


        lh=legend(lgnd,'Location','southeast');
        lpos = lh.Position;
        lpos(1) = lpos(1) +0.1;
        lh.Position = lpos;
        title(my_rois_plus{a})
    end
    
    saveas(fig2,fullfile(resultsSave,['scatter_',timePointName{tPoint},'.svg']))
    close(fig2)
    



%% PCA time - 3d plot
    projectionCombo = cat(5,roi_wm_projection_onto_Stable,roi_wm_projection_onto_Stable_orientation,roi_wm_projection_onto_Stable_noise);
    nameCombo = {'NoDist','OrDist','NoiseDist'};

    plotTimePointsZ = [2,5,8,11,14,17];

    for a=1:length(my_rois_plus)
        for b=1:3
            fig = figure(1)
            
            set(fig, 'PaperUnits', 'inches');
            x_width=1 ;y_width=1;
            set(fig, 'PaperPosition', [0 0 x_width y_width]); %
            fig.Renderer = 'painters';

            p1=plot3(squeeze(projectionCombo(1,1,:,a,b)),squeeze(projectionCombo(1,2,:,a,b)),plotTimePointsZ)
            p1.LineWidth = 0.5;  
            p1.Color = colors(1,:);
            hold on

            p2=plot3(squeeze(projectionCombo(2,1,:,a,b)),squeeze(projectionCombo(2,2,:,a,b)),plotTimePointsZ)
            p2.LineWidth = 0.5;  
            p2.Color = colors(2,:);
            hold on

            p3=plot3(squeeze(projectionCombo(3,1,:,a,b)),squeeze(projectionCombo(3,2,:,a,b)),plotTimePointsZ)
            p3.LineWidth = 0.5;  
            p3.Color = colors(3,:);
            hold on

            p4=plot3(squeeze(projectionCombo(4,1,:,a,b)),squeeze(projectionCombo(4,2,:,a,b)),plotTimePointsZ)
            p4.LineWidth = 0.5;  
            p4.Color = colors(4,:);
            hold on


            if numBins==6
                p5=plot3(squeeze(projectionCombo(5,1,:,a,b)),squeeze(projectionCombo(5,2,:,a,b)),plotTimePointsZ)
                p5.LineWidth = 0.5;  
                p5.Color = colors(5,:);
                hold on
    
                p6=plot3(squeeze(projectionCombo(6,1,:,a,b)),squeeze(projectionCombo(6,2,:,a,b)),plotTimePointsZ)
                p6.LineWidth = 0.5;  
                p6.Color = colors(6,:);
                hold on
            end

            rectangle('Position',[-10 -10 20 20],'FaceColor',[.7 .7 .7])
            
            zlim([0 lengthDelay+1])  
            xlim([-10, 10])
            ylim([-10, 10])
            grid off
            xticks([-10 0 10])
            yticks([-10 0 10])
            
            if b==1
                xlabel(strjoin([string(round(perc_explained_no(1,tPoint,rr))),"%"],""))
                ylabel(strjoin([string(round(perc_explained_no(2,tPoint,rr))),"%"],""))
            elseif b==2
                xlabel(strjoin([string(round(perc_explained_or(1,tPoint,rr))),"%"],""))
                ylabel(strjoin([string(round(perc_explained_or(2,tPoint,rr))),"%"],""))
            elseif b==3
                xlabel(strjoin([string(round(perc_explained_noise(1,tPoint,rr))),"%"],""))
                ylabel(strjoin([string(round(perc_explained_noise(2,tPoint,rr))),"%"],""))
            end
            
            zticks([6 12 18])
            ax = gca;
            ax.FontSize = 6; 

            title([my_rois_plus{a},' ',timePointName{tPoint},' ',nameCombo{b}])
            
            saveas(fig,fullfile(resultsSave,['threeDplot_PCA_time_',my_rois_plus{a},timePointName{tPoint},nameCombo{b},'.svg']))    
            close(fig)
        end
    end

end

end