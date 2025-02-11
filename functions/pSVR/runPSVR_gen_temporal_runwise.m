function [genMatrix_fca, genMatrix_bfca] = runPSVR_gen_temporal_runwise(data_train, data_test,labels_train,labels_test)
% You need: 'data' (preprocessed), in the form of n_trials x n_voxels x n_runs x n_TRs
%           'labels', in degrees btw. 0-360°, in the format n_trials x n_runs
%           a 3D binary ROI 'mask', with the dimensions of the original .nii
%           image
%
% Output:   sin_pred, cos_pred, ang_pred: cell arrays with the sine, cosine
%           and angular predictions, respectively
%           true_labels: cell array with the true labels used for
%           reconstruction
   
    %% Prepare necessary variables
    numTRs = size(data_train,1);
    warning("TR changed to 0 instead of 2")
    TR_gen = 0; %three TRs being meaned
    
    numTRs = numTRs-TR_gen;
    
    testIndex = zeros(size(data_train,3),3);
    testIndex(1:size(data_train,3)/3,1) = 1;
    testIndex(size(data_train,3)/3+1:size(data_train,3)/3*2,2) = 1;
    testIndex(size(data_train,3)/3*2+1:size(data_train,3),3) = 1;
    
    trainIndex = ~testIndex;
    testIndex = ~trainIndex;

    genMatrix_fca = zeros([numTRs,numTRs,3]);
    genMatrix_bfca = zeros([numTRs,numTRs,3]);
    
    for TR_test = 1:numTRs
        for TR_train =1:numTRs
            
            bfca_session = [];
            fca_session = [];
            
            for session = 1:3
            % get data dimensions
                data_train_c = data_train(TR_train:TR_train+TR_gen,:,trainIndex(:,session));
                data_train_c = mean(data_train_c,1);
                data_train_c = squeeze(data_train_c)';
                data_test_c = data_test(TR_test:TR_test+TR_gen,:,testIndex(:,session));
                data_test_c = mean(data_test_c,1);
                data_test_c = squeeze(data_test_c)';

                data = [data_train_c;data_test_c];

                % prepare labels
                labels_train_session = labels_train(trainIndex(:,session));
                labels_test_session = labels_test(testIndex(:,session));
                labels = [labels_train_session;labels_test_session];                                        % vecotrize labels
                nan_ind = isnan(labels);                                                    % find inidces of missing trials
                labels(nan_ind) = [];                                                       % remove missing trials
                label_rad = deg2rad(labels)-pi;                                             % turn to rad and shift to range -pi:pi
                label_sin = sin(label_rad);                                                 % extract sin and cos components
                label_cos = cos(label_rad);
                label = mat2cell([label_sin, label_cos], ones(numel(label_rad),1));         % turn to cell array for multitarget regression

                % prepare n_files and chunk for TDT
                chunk = [repmat(1,1,numel(labels_train_session)),repmat(2,1,numel(labels_test_session))]';


                %% Run reconstruction via TDT

                % grab data of curent TR and bring into right format for TDT
                yt = data;

                % reassign passed_data and results to avoid looping issues
                passed_data = [];
                results = [];

                % set up cfg for TDT
                cfg = decoding_defaults;
                cfg.results.overwrite = 1;
                cfg.analysis = 'ROI';
                cfg.multitarget = 1;
                cfg.decoding.method = 'regression';
                cfg.decoding.train.classification.model_parameters = '-s 4 -t 2 -c 1 -n 0.5 -b 0 -q';   % libsvm SVR parameters
                cfg.decoding.software = 'libsvm_multitarget';
                cfg.results.output = {'predicted_labels_multitarget'};
                cfg.plot_selected_voxels = 0;
                cfg.plot_design = 0;
                cfg.results.write = 0;
                cfg.scale.method = 'min0max1';                                          % rescaling to range [0 1]
                cfg.scale.estimation = 'across';                                        % scaling estimated on training data and applied to test data

                % fill passed_data
                passed_data.data = yt;

                [passed_data,cfg] = fill_passed_data(passed_data, cfg, label, chunk);

                % create design (and plot once if required)
                cfg.design = make_design_cv(cfg);

                % peform decoding
                [results, cfg, passed_data] = decoding(cfg, passed_data);

                % get sin and cos predictions from result structure
                sin_y = results.predicted_labels_multitarget.output.model{1}.predicted_labels(numel(labels_train_session)+1:end);
                cos_y = results.predicted_labels_multitarget.output.model{2}.predicted_labels(numel(labels_train_session)+1:end);

                % reconstruct predicted angular label using four-quadran arctangent
                ang_y = atan2(sin_y,cos_y);

                % get true labels
                true_labels=label_rad(numel(labels_train_session)+1:end);

                genMatrix_bfca(TR_test,TR_train,session) = bal_norm_circ_resp_dev(ang_y, true_labels, 'trapz').*100;
                genMatrix_fca(TR_test,TR_train,session) = avg_norm_circ_resp_dev(ang_y, true_labels).*100;
            end                
        end
    end
    
    genMatrix_bfca = {genMatrix_bfca};
    genMatrix_fca = {genMatrix_fca};
end