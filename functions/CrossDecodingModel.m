rng(13080)
path2data = "../derivates";
results_new = fullfile(path2data,"Simulations");
mkdir(results_new)

orientations = [1:180]*2;
orientations_sin = sin(deg2rad(orientations));
orientations_cos = cos(deg2rad(orientations));

num_neurons = 100;

dataset_expand = [];

for i = 0:num_neurons-1
    shifted_cos = circshift(orientations_cos, i);
    shifted_sin = circshift(orientations_sin, i);
    
    dataset_expand = [dataset_expand; shifted_cos; shifted_sin];
end

%% Cross-SNR pSVR 
numIter = 7;
SNR_values = [0.1:0.1:5,10,20];
SNR_results = zeros(numIter,length(SNR_values));
num_rows_to_select = 108;
total_rows = size(dataset_expand, 2);

SNRValues_to_sample = [0.5:0.5:10,20,50];
crossSNR_decodingMatrix  = zeros(length(SNRValues_to_sample),length(SNRValues_to_sample),numIter);

for iter = 1:numIter
    
    random_indices = randperm(total_rows, num_rows_to_select);
    selected_trials = dataset_expand(:,random_indices);
    orientations_selected = orientations(random_indices);
    
    num_rows = size(selected_trials, 1);
    random_order = randperm(num_rows);
    permuted_dataset = selected_trials(random_order, :);
    
    disp(iter)
    
    for noise_now_train = 1:length(SNRValues_to_sample)
        for noise_now_test = 1:length(SNRValues_to_sample)
    
            randRange_train = -1+rand(size(permuted_dataset))*(1+1);
            randRange_test = -1+rand(size(permuted_dataset))*(1+1);
        
            data_train = zscore(permuted_dataset + SNRValues_to_sample(noise_now_train)*randRange_train,[],2);
            data_test = zscore(permuted_dataset + SNRValues_to_sample(noise_now_test)*randRange_test,[],2);
        
            data = [data_train,data_test]';
        
            % prepare labels
            labels_train_session = orientations_selected';
            labels_test_session = orientations_selected';
            
            labels = [labels_train_session;labels_test_session];                                        % vecotrize labels
        
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
            cfg.verbose = 0;
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
        
            crossSNR_decodingMatrix(noise_now_train,noise_now_test,iter) = avg_norm_circ_resp_dev(ang_y, true_labels).*100;
        end
    end                
end

save(fullfile(results_new,['crossSNR_simulation.mat']),'crossSNR_decodingMatrix');

Dynamics=dynamicCoding_signInversion(crossSNR_decodingMatrix);


colorsMap = cbrewer('div', 'RdYlBu', 64);
colorsMap = flipud(colorsMap); 
colorsMap(colorsMap<0) = 0;
colorsMap(colorsMap>1)=1;


fig = figure(2)
imagesc(mean(crossSNR_decodingMatrix,3))
%runBoundary(Dynamics,'b')
set(gca, 'ydir', 'normal');
axis square;
xticks(1:length(SNRValues_to_sample))
xticklabels(split(num2str(SNRValues_to_sample)))
yticks(1:length(SNRValues_to_sample))
yticklabels(split(num2str(SNRValues_to_sample)))
xlabel("Noise amplitude - training data")
ylabel("Noise amplitude - test data")
handles = colorbar;
caxis([45 100])
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Label.String = '% Decoding accuracy';
handles.TickDirection = 'out';
handles.Box = 'off';
handles.Label.String = '% Decoding accuracy';
handles.Label.FontSize = 5;
drawnow;
colormap(colorsMap)
saveas(fig,fullfile(results_new,"Cross-SNR comparison.svg"))    



