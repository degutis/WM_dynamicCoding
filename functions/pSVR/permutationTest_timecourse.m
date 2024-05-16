function cluster_above_chance = permutationTest_timecourse(M)

    rng(round(M(1,1)*M(1,1)))

    timeDim = size(M,2);
    M_shuffle = zeros(timeDim,1000);
    numSubs = 1:size(M,1);

    PlusMin = [-1 1];
    
    for perm=1:1000    
        
        matrixPerm = zeros(timeDim,1);
    
        for aa=1:timeDim
            matrixPerm(aa) = M(randsample(numSubs,1),aa)*randsample(PlusMin,1);
        end
        
        M_shuffle(:,perm) = matrixPerm;
    end
    
    [clusters_sup, p_values_sup, ~, ~] = permutest(M', M_shuffle, 0, 0.05, 1000, 0);

    cluster_above_chance = zeros(timeDim,1);
    for c=1:length(clusters_sup)
        if p_values_sup(c)<0.05
            cluster_above_chance(clusters_sup{c})=1;
        end
    end
    
end