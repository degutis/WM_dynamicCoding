function [Angles,currentCond] = runPCAangleSplit(WMcell,sub,roi,timePointCell,timePointIndex)
    
    cond1_same = WMcell{sub,roi,1};
    cond2_same = WMcell{sub,roi,2};
    
    Angles = zeros(length(timePointIndex),1);
    for t=1:length(timePointIndex)
    
        time1_train = mean(cond1_same(:,:,timePointCell{t,1}),3);
        time1_train = time1_train - mean(time1_train,1);

        time1_test = mean(cond2_same(:,:,timePointCell{t,1}),3);
        time1_test = time1_test - mean(time1_test,1);

        time2_train = mean(cond1_same(:,:,timePointCell{t,2}),3);
        time2_train = time2_train - mean(time2_train,1);
        
        time2_test = mean(cond2_same(:,:,timePointCell{t,2}),3);
        time2_test = time2_test - mean(time2_test,1);

        currentCond = zeros(4,1);
        
        for cond=1:4
    
            if cond==1
                combo = [time1_train;time1_test];
            elseif cond==2
                combo = [time2_train;time2_test];
            elseif cond==3
                combo = [time1_train;time2_test];                   
            elseif cond==4
                combo = [time2_train;time1_test];  
            end
    
            [~,W1]=pca(combo(1:size(combo,1)/2,:)','NumComponents',2);
            [~,W2]=pca(combo(size(combo,1)/2+1:end,:)','NumComponents',2);
    
            Theta = rad2deg(subspace(W1,W2));
    
            currentCond(cond) = Theta;
        end
    
       Angles(t) = mean([[currentCond(3)-currentCond(1)],[currentCond(4)-currentCond(2)]]);

    end