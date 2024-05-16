function pSVR_temporal(experiment,subject,voxSelect,expandIPS)

global path2data

rng(subject)

if voxSelect==1
    nameAdd = 'ReliableVox';
end

results_new = fullfile(path2data,'pSVR_analysis',experiment,'temporal',nameAdd);
mkdir(results_new)

results_new_individ = fullfile(path2data,'pSVR_analysis',experiment,'temporal',nameAdd,'individ');
mkdir(results_new_individ)

sampleFiles = dir(fullfile('../raw_data',experiment,'SampleFile_*'));
timeCourseFiles = dir(fullfile('../raw_data',experiment,'TimeCourses_*'));

my_rois_plus = {'V1','V2','V3','V3AB','V4','IPS','LO1','LO2'};

trs_to_recon = [-1 1:22];  % 1 tr before stim onset, 22 tr's after
num_trs_to_recon = length(trs_to_recon);
num_trials = 108*3;

wmnoDist_fcaMatrix = cell(numel(sampleFiles),numel(my_rois_plus));
wmOr_fcaMatrix = cell(numel(sampleFiles),numel(my_rois_plus));
wmNoise_fcaMatrix= cell(numel(sampleFiles),numel(my_rois_plus));

wmOronMem_fcaMatrix = cell(numel(sampleFiles),numel(my_rois_plus));
wmNoiseonMem_fcaMatrix= cell(numel(sampleFiles),numel(my_rois_plus));
wmDist_fcaMatrix = cell(numel(sampleFiles),numel(my_rois_plus));

wmDist_fcaMatrix_distOnDist = cell(numel(sampleFiles),numel(my_rois_plus));
wmDist_fcaMatrix_genDisttoOr = cell(numel(sampleFiles),numel(my_rois_plus));


if expandIPS==1
    indexVox = load(fullfile(path2data,'correlation_analysis_GLM',experiment,'fwhm_30_shuffle_0_IPS','FullResults_correlationAnalysis.mat'),'sub_voxSelect');
else
    error("Not implemented")
end

for sub = subject:subject

    disp(['temporal gen pSVR analysis for ',experiment,' subject ', num2str(sub)]);
    % load data and such
    load(fullfile(sampleFiles(sub).folder,sampleFiles(sub).name));
    load(fullfile(timeCourseFiles(sub).folder,timeCourseFiles(sub).name));

    if strcmp('E2',experiment)
        TimeCoursesLocalizer = TimeCoursesLocalizer2;
        AvgResponsesLocalizer = AvgResponsesLocalizer2;
    end

    if expandIPS==1
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

    wmTrials = reshape(wmTrials,length(wmTrials)/3,3)*2;
    wmTrials_noise = reshape(wmTrials_noise,length(wmTrials_noise)/3,3)*2;
    wmTrials_orDist = reshape(wmTrials_orDist,length(wmTrials_orDist)/3,3)*2;
    distractor_WM = reshape(distractor_WM,length(distractor_WM)/3,3)*2;
    
    for roi = 1:length(my_rois_plus)

        % grab the time course in this roi (localized voxels only) as well
        % as the avg response on localizer trials
        my_TC = TimeCoursesMain{roi};

        % Make a matrix of num_trs_to_recon by num_voxels by num_trials
        rsp_wm = NaN(num_trs_to_recon,size(my_TC,2),num_trials);
        for trial = 1:num_trials
            rsp_wm(:,:,trial) = my_TC(TrialStartInd_main(trial)+trs_to_recon(1):TrialStartInd_main(trial)+trs_to_recon(end)-1,:);
        end

        wm_noDist = rsp_wm(:,:,wmTrialIndex);
        wm_orDist = rsp_wm(:,:,wmTrialIndex_orDist);
        wm_noiseDist = rsp_wm(:,:,wmTrialIndex_noise);
                
        [wmnoDist_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_noDist,wm_noDist,wmTrials(:),wmTrials(:));              
 
        [wmOr_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_orDist,wm_orDist,wmTrials_orDist(:),wmTrials_orDist(:));
        
        [wmNoise_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_noiseDist,wm_noiseDist,wmTrials_noise(:),wmTrials_noise(:));
         

        
        [wmOronMem_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_noDist,wm_orDist,wmTrials(:),wmTrials_orDist(:));
         
        [wmNoiseonMem_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_noDist,wm_noiseDist,wmTrials(:),wmTrials_noise(:));

        [wmOronMem_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_noDist,wm_orDist,wmTrials(:),wmTrials_orDist(:));        
        [wmNoiseonMem_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_noDist,wm_noiseDist,wmTrials(:),wmTrials_noise(:));



        [wmDist_fcaMatrix(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_orDist,wm_orDist,wmTrials_orDist(:),distractor_WM(:));

        [wmDist_fcaMatrix_distOnDist(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_orDist,wm_orDist,distractor_WM(:),distractor_WM(:));
        [wmDist_fcaMatrix_genDisttoOr(sub,roi),~] = runPSVR_gen_temporal_runwise(wm_orDist,wm_orDist,distractor_WM(:),wmTrials_orDist(:));

    end
end

save(fullfile(results_new,['temporal_runwise_PSVM_',num2str(subject),'.mat']),'wmnoDist_fcaMatrix','wmOr_fcaMatrix','wmNoise_fcaMatrix', 'wmOronMem_fcaMatrix', 'wmNoiseonMem_fcaMatrix','wmDist_fcaMatrix')
%save(fullfile(results_new,['temporal_Extra_NoDistOnDist_runwise_PSVM_',num2str(subject),'.mat']),'wmDist_fcaMatrix_distOnDist', 'wmDist_fcaMatrix_genDisttoOr')

end