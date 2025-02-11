function simulationNoiseAddedRealData(trialType)

rng(111)
path2data = "../derivates";
experiment="E1";

results_new = fullfile(path2data,"Simulations");
mkdir(results_new)

sampleFiles = dir(fullfile('../raw_data',experiment,'SampleFile_*'));
timeCourseFiles = dir(fullfile('../raw_data',experiment,'TimeCourses_*'));

my_rois_plus = {'V1'};

trs_to_recon = [-1 1:22];  % 1 tr before stim onset, 22 tr's after
num_trs_to_recon = length(trs_to_recon);
num_trials = 108*3;

wm_fcaMatrix = cell(numel(sampleFiles),numel(my_rois_plus));
indexVox = load(fullfile(path2data,'correlation_analysis_GLM',experiment,'fwhm_30_shuffle_0_IPS','FullResults_correlationAnalysis.mat'),'sub_voxSelect');

SNRValues_to_sample = [0.5:0.5:3, flip(0.5:0.5:3)];


for sub = 1:6
    disp(num2str(sub))
    load(fullfile(sampleFiles(sub).folder,sampleFiles(sub).name));
    load(fullfile(timeCourseFiles(sub).folder,timeCourseFiles(sub).name));


    TimeCoursesLocalizer = [TimeCoursesLocalizer(1:5),{cat(2,TimeCoursesLocalizer{6:9})},TimeCoursesLocalizer(10:11)];
    TimeCoursesMain = [TimeCoursesMain(1:5),{cat(2,TimeCoursesMain{6:9})},TimeCoursesMain(10:11)];
        
    AvgResponsesLocalizer = [AvgResponsesLocalizer(1:5),{cat(2,AvgResponsesLocalizer{6:9})},AvgResponsesLocalizer(10:11)];
    AvgResponsesMain = [AvgResponsesMain(1:5),{cat(2,AvgResponsesMain{6:9})},AvgResponsesMain(10:11)];
    
    wmTrialIndex = InfoMain(:,3)==0;
    wmTrials = InfoMain(InfoMain(:,3)==0,2);

    wmTrialIndex_noise = isnan(InfoMain(:,3));
    wmTrialIndex_orDist = InfoMain(:,3)>0;

    wmTrials_noise = InfoMain(isnan(InfoMain(:,3)),2);
    wmTrials_orDist = InfoMain(InfoMain(:,3)>0,2);
   
   
    for roiInit = 1:numel(TimeCoursesMain)            
        TimeCoursesMain(roiInit)={TimeCoursesMain{roiInit}(:,indexVox.sub_voxSelect{sub,roiInit})};
    end
   
    wmTrials = reshape(wmTrials,length(wmTrials)/3,3)*2;
    wmTrials_noise = reshape(wmTrials_noise,length(wmTrials_noise)/3,3)*2;
    wmTrials_orDist = reshape(wmTrials_orDist,length(wmTrials_orDist)/3,3)*2;

    for roi = 1:length(my_rois_plus)

        my_TC = TimeCoursesMain{roi};

        rsp_wm = NaN(num_trs_to_recon,size(my_TC,2),num_trials);
        for trial = 1:num_trials
            rsp_wm(:,:,trial) = my_TC(TrialStartInd_main(trial)+trs_to_recon(1):TrialStartInd_main(trial)+trs_to_recon(end)-1,:);
        end

        if strcmp(trialType,"noDist")
            wm_noDist = rsp_wm(:,:,wmTrialIndex);
            wm_noDist = mean(wm_noDist(8:12,:,:),1);
            wm_full = zeros(length(SNRValues_to_sample), size(wm_noDist,2), size(wm_noDist,3));
            
            for n = 1:length(SNRValues_to_sample)
                randRange_train = -1+rand(size(wm_noDist))*(1+1);
                wm_full(n,:,:) = wm_noDist + SNRValues_to_sample(n)*randRange_train;
            end

            [wm_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_full,wm_full,wmTrials(:),wmTrials(:));               

            
        elseif strcmp(trialType, "orDist")
            wm_orDist = rsp_wm(:,:,wmTrialIndex_orDist);
            wm_orDist = mean(wm_orDist(8:12,:,:),1);
            wm_full = zeros(length(SNRValues_to_sample), size(wm_orDist,2), size(wm_orDist,3));
        
            for n = 1:length(SNRValues_to_sample)
                randRange_train = -1+rand(size(wm_orDist))*(1+1);
                wm_full(n,:,:) = wm_orDist + SNRValues_to_sample(n)*randRange_train;
            end

            [wm_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_full,wm_full,wmTrials_orDist(:),wmTrials_orDist(:));               
        
        
        elseif strcmp(trialType, "noiseDist")
            wm_noiseDist = rsp_wm(:,:,wmTrialIndex_noise);
            wm_noiseDist = mean(wm_noiseDist(8:12,:,:),1);
            wm_full = zeros(length(SNRValues_to_sample), size(wm_noiseDist,2), size(wm_noiseDist,3));            
            
            for n = 1:length(SNRValues_to_sample)
                randRange_train = -1+rand(size(wm_noiseDist))*(1+1);
                wm_full(n,:,:) = wm_noiseDist + SNRValues_to_sample(n)*randRange_train;
            end

            [wm_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_full,wm_full,wmTrials_noise(:),wmTrials_noise(:));                       
        end
            
    end
end

wmNoDist_reliable = cell(6,1);

for f=1:6
    wmNoDist_reliable(f,:) = cellfun(@(x) mean(x(:,:,:),3), wm_fcaMatrix(f,:), 'UniformOutput', false);
end

combo = cat(3, wmNoDist_reliable{:});


save(fullfile(results_new,'simulationRealData_'+trialType+'.mat'),'combo');

mean_combo = mean(combo,3);
Cluster=permutationTest_cluster_psvr(combo-50);
Dynamics=dynamicCoding_signInversion(combo);

colorsMap = cbrewer('div', 'RdYlBu', 64);
colorsMap = flipud(colorsMap); 
colorsMap(colorsMap<0) = 0;
colorsMap(colorsMap>1)=1;


fig = figure(2)
imagesc(mean_combo)
runBoundary(Cluster)
runBoundary(Dynamics,'b')
set(gca, 'ydir', 'normal');
axis square;
xticks(1:length(SNRValues_to_sample))
xticklabels(split(num2str(SNRValues_to_sample)))
yticks(1:length(SNRValues_to_sample))
yticklabels(split(num2str(SNRValues_to_sample)))
xlabel("Noise amplitude - training data")
ylabel("Noise amplitude - test data")
handles = colorbar;
caxis([45 75])
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Label.String = '% Decoding accuracy';
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Label.String = '% Decoding accuracy';
handles.Label.FontSize = 5;
drawnow;
colormap(colorsMap)
saveas(fig,fullfile(results_new,'SNR_comp_'+trialType+'.svg'))    

end