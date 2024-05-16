function [orientationRespROIs,orienationSNR] = estimateTuningCurveGLM(responseDist,trials,bins,isTable)
    orienationSNR = cell(1,size(responseDist,2));
    orientationRespROIs = cell(1,size(responseDist,2));
    X = zeros(length(trials),length(bins));
    for i =1:length(bins)
        X(trials==bins(i),i) =1;
    end
        
    
    if isTable==1
        for column = 1:size(responseDist,2)
            Y = responseDist.(column);
            V = eye(size(Y,1));
            B = ME_GLM(Y,X,V); % trial-wise response, orientation on a given trial, covariance matrix (no covar- eye)
            [s, SNR] = ME_GLM_SNR(Y, X, V, B);
            orienationSNR(column) = {SNR'};        
            orientationRespROIs(column) = {B'};
        end
    else
        Y = responseDist;
        V = eye(size(Y,1));
        B = ME_GLM(Y,X,V); % trial-wise response, orientation on a given trial, covariance matrix (no covar- eye)
        [s, SNR] = ME_GLM_SNR(Y, X, V, B);
        orienationSNR = SNR';        
        orientationRespROIs = B';
    end

    
end



