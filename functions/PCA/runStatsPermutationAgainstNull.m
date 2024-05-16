function [p_value]=runStatsPermutationAgainstNull(X,Y,timePointIndices,greaterThan,half,seed)

if half==0

    numRois = size(X,2);
    numSub=size(X,1);
    numIter = 10000;
    
    p_value= zeros(size(X,3),size(X,4), numRois);
    
    X = mean(mean(X,6),5);
    Y = mean(mean(Y,5),4);
    
    numTests=length(timePointIndices);
    
    for roi=1:numRois
        
        nonNull_rows = zeros(numSub,numTests);
        nonNull = squeeze(X(:,roi,:,:));
    
        for t = 1:length(timePointIndices)
            nonNull_rows(:,t) = nonNull(:,timePointIndices{t}(1),timePointIndices{t}(2));
        end
    
        null = squeeze(Y(:,roi,:));
        null_current = zeros(numSub,numTests);
    
        pVal = zeros(length(timePointIndices),1);
    
        for t = 1:length(timePointIndices)
            if greaterThan==1
                null_current(:,t)=max(null(:,timePointIndices{t}(:))')';
            else
                null_current(:,t)=min(null(:,timePointIndices{t}(:))')';
            end
    
            stat = mean(nonNull_rows(:,t)-null_current(:,t));
    
            rng(seed+t*roi)
            
            bTestStat = zeros(numIter,1);
            for b=1:numIter
                bFull = [nonNull_rows(:,t);null_current(:,t)];
                permBFull = bFull(randperm(length(bFull)));
    %             [~,~,~,stats] = ttest(permBFull(1:numSub),permBFull(numSub+1:end));
    %             bTestStat(b) = stats.tstat; 
                bTestStat(b) = mean(permBFull(1:numSub) - permBFull(numSub+1:end));
            end
    
            if greaterThan==1
                pVal(t) = sum(bTestStat>stat)/numIter;
            else
                pVal(t) = sum(bTestStat<stat)/numIter;
            end            
        end
    
        [~,~,adj_p] = fdr_bh(pVal,0.05,'pdep','yes');
    
        for t3 = 1:length(timePointIndices)  
            p_value(timePointIndices{t3}(1),timePointIndices{t3}(2),roi) = adj_p(t3);
        end
    end

elseif half==1
    
    numRois = size(X,2);
    numSub=size(X,1);
    numIter = size(Y,4);
    
    p_value= zeros(size(X,3),size(X,4), numRois);
    Y = mean(Y,5);    
    numTests=length(timePointIndices);

    for roi=1:numRois
        
        nonNull_rows = zeros(numSub,numTests);
        nonNull = squeeze(X(:,roi,:,:));
    
        for t = 1:length(timePointIndices)
            nonNull_rows(:,t) = nonNull(:,timePointIndices{t}(1),timePointIndices{t}(2));
        end
    
        null = squeeze(Y(:,roi,:,:));
        null_current = zeros(numSub,numTests,numIter);
    
        pVal = zeros(length(timePointIndices),1);
    
        for t = 1:length(timePointIndices)
            if greaterThan==1
                idx = mean(mean(null(:,timePointIndices{t}(1),:),3))>mean(mean(null(:,timePointIndices{t}(2),:),3));
            else
                idx = mean(mean(null(:,timePointIndices{t}(1),:),3))<mean(mean(null(:,timePointIndices{t}(2),:),3));
            end
            if idx==0
                null_current(:,t,:) = null(:,timePointIndices{t}(2),:);
            else
                null_current(:,t,:) = null(:,timePointIndices{t}(1),:);
            end
       
            stat = zeros(numIter,1);
            for i=1:numIter
                stat(i) = mean(nonNull_rows(:,t)-null_current(:,t,i));
            end
        
            if greaterThan==1
                pVal(t) = sum(0>=stat)/numIter;
            else
                pVal(t) = sum(0<=stat)/numIter;
            end            
        end
    
        [~,~,adj_p] = fdr_bh(pVal,0.05,'pdep','yes');
    
        for t3 = 1:length(timePointIndices)  
            p_value(timePointIndices{t3}(1),timePointIndices{t3}(2),roi) = adj_p(t3);
        end
    end

elseif half==2

    numRois = size(X,2);
    numSub=size(X,1);
    numIter = 10000;
    
    p_value= zeros(size(X,3),size(X,4), numRois);
    
    X = mean(mean(X,6),5);
    Y = mean(mean(Y,6),5);
    
    numTests=length(timePointIndices);
    
    for roi=1:numRois
        
        nonNull_rows = zeros(numSub,numTests);
        nonNull = squeeze(X(:,roi,:,:));
    
        for t = 1:length(timePointIndices)
            nonNull_rows(:,t) = nonNull(:,timePointIndices{t}(1),timePointIndices{t}(2));
        end
    
        nonNull2_rows = zeros(numSub,numTests);
        nonNull2 = squeeze(Y(:,roi,:,:));
    
        for t = 1:length(timePointIndices)
            nonNull_rows(:,t) = nonNull(:,timePointIndices{t}(1),timePointIndices{t}(2));
            nonNull2_rows(:,t) = nonNull2(:,timePointIndices{t}(1),timePointIndices{t}(2));
        end
    
        pVal = zeros(length(timePointIndices),1);
    
        for t = 1:length(timePointIndices)
    
            stat = mean(nonNull_rows(:,t)-nonNull2_rows(:,t));
    
            rng(seed+t*roi)
            
            bTestStat = zeros(numIter,1);
            for b=1:numIter
                bFull = [nonNull_rows(:,t);nonNull2_rows(:,t)];
                permBFull = bFull(randperm(length(bFull)));
                bTestStat(b) = mean(permBFull(1:numSub) - permBFull(numSub+1:end));
            end
    
            pVal1 = sum(bTestStat>stat)/numIter;
            pVal2 = sum(bTestStat<stat)/numIter;

            pVal(t) = min(pVal1,pVal2);
        end
    
        [~,~,adj_p] = fdr_bh(pVal,0.025,'pdep','yes');
    
        for t3 = 1:length(timePointIndices)  
            p_value(timePointIndices{t3}(1),timePointIndices{t3}(2),roi) = adj_p(t3);
        end
    end


elseif half==3

    numRois = size(X,2);
    numSub=size(X,1);
    numIter = 10000;
    maxMatrix = max(max(cat(1,timePointIndices{:})));
    
    p_value= zeros(maxMatrix,maxMatrix, numRois,3);
           
    for roi=1:numRois
                
        pVal = zeros(length(timePointIndices),3);
        
        for tType = 1:size(pVal,2)
            for t = 1:size(pVal,1)
        
                currentValues = X(:,roi,(tType-1)*size(pVal,1)+t); 
                stat = mean(currentValues);

                rng(seed+t*roi*tType)

                bTestStat = zeros(numIter,1);
                for b=1:numIter
                    bTestStat(b) = mean(currentValues.*datasample([-1,1],numSub,'Replace',true)');
                end

                if greaterThan==1
                    pVal(t,tType) = sum(bTestStat>stat)/numIter;
                else
                    pVal(t,tType) = sum(bTestStat<stat)/numIter;
                end            
            end
        end
    
        [~,~,adj_p] = fdr_bh(pVal,0.05,'pdep','yes');
    
        for type2 = 1:size(pVal,2)
            for t3 = 1:length(timePointIndices)  
                p_value(timePointIndices{t3}(1),timePointIndices{t3}(2),roi,type2) = adj_p(t3,type2);
            end
        end
    end

elseif half==4

    numRois = size(X,2);
    numSub=size(X,1);
    numIter = 10000;
        
    pVal = zeros(numRois,1);
    for roi=1:numRois
             
        currentValues = X(:,roi); 
        stat = mean(currentValues);
                
        rng(seed+roi)
            
        bTestStat = zeros(numIter,1);
        for b=1:numIter
            bTestStat(b) = mean(currentValues.*datasample([-1,1],numSub,'Replace',true)');
        end
    
        if greaterThan==1
            pVal(roi) = sum(bTestStat>stat)/numIter;
        else
            pVal(roi) = sum(bTestStat<stat)/numIter;
        end            
    end
       
    [~,~,adj_p] = fdr_bh(pVal,0.05,'pdep','yes');
    
    p_value = adj_p; 
    
end





