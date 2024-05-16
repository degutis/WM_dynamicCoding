function current_cell = extractTrialsMean_splitHalf(delayPeriodSelected,index,rsp_wm,wmTrials_binned,binsVal,nameCV,numberCV,seed)

    rng(seed)

    current_cell = cell(2,1);

    wm = rsp_wm(delayPeriodSelected,:,index);    
    
    wm_perm = permute(wm,[1,3,2]);
    wm_partition = cvpartition(wmTrials_binned,nameCV,numberCV,'Stratify',true);
    ind_wm_train = training(wm_partition,1);
    ind_wm_test = test(wm_partition,1); 

    currentWM_full = zeros(length(binsVal),size(wm_perm,3),size(wm_perm,1));

    for ind=1:2
        if ind==1
            indexCurrent_wm=ind_wm_train;
        elseif ind==2
            indexCurrent_wm=ind_wm_test;
        end

        for tt=1:size(wm,1)

            currentWM = squeeze(wm_perm(tt,indexCurrent_wm,:));
            try
                currentWM = [mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(1),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(2),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(3),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(4),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(5),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(6),:))];
            catch
                currentWM = [mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(1),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(2),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(3),:));...
                    mean(currentWM(wmTrials_binned(indexCurrent_wm)==binsVal(4),:))];
            end

            currentWM_full(:,:,tt) = currentWM;

        end 
        current_cell(ind) = {currentWM_full}; 
    end
