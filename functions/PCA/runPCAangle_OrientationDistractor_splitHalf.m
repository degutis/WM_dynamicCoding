function Angle = runPCAangle_OrientationDistractor_splitHalf(OrCell,DistCell,sub,roi,iter)


    session1_wm = OrCell{sub,roi,1,iter};
    session2_wm = OrCell{sub,roi,2,iter};
    
    session1_dist = DistCell{sub,roi,1,iter};
    session2_dist = DistCell{sub,roi,2,iter};
       
    wm_session1 = mean(session1_wm,3);
    wm_session1 = wm_session1 - mean(wm_session1,1);

    wm_session2 = mean(session2_wm,3);
    wm_session2 = wm_session2 - mean(wm_session2,1);

    dist_session1 = mean(session1_dist,3);
    dist_session1 = dist_session1 - mean(dist_session1,1);

    dist_session2 = mean(session2_dist,3);
    dist_session2 = dist_session2 - mean(dist_session2,1);
    
    currentCond = zeros(4,1);

    for cond=1:4

        % within subspaces

        if cond==1
            combo = [wm_session1;wm_session2];
        elseif cond==2
            combo = [dist_session1;dist_session2];

        % across subspaces    

        elseif cond==3
            combo = [wm_session1;dist_session2];  
        elseif cond==4
            combo = [wm_session2;dist_session1];  
        end

        [~,W1]=pca(combo(1:size(combo,1)/2,:)','NumComponents',2);
        [~,W2]=pca(combo(size(combo,1)/2+1:end,:)','NumComponents',2);

        Theta = rad2deg(subspace(W1,W2));

        currentCond(cond) = Theta;
    end

    Angle = mean(currentCond(3:4)) - mean(currentCond(1:2));
        
end
