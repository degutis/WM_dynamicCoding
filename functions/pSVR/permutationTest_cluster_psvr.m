function cluster_above_chance = permutationTest_cluster_psvr(M)

    rng(abs(round(M(1,1,1)*M(1,2,1))+1))

    timeDim = size(M,1);
    M_shuffle = zeros(timeDim,timeDim,1000);
    numSubs = 1:size(M,3);

    PlusMin = [-1 1];
    
    for perm=1:1000    
        
        matrixPerm = zeros(timeDim,timeDim);
    
        for aa=1:timeDim
            for bb=1:timeDim
                matrixPerm(aa,bb) = M(aa,bb,randsample(numSubs,1))*randsample(PlusMin,1);
            end
        end

        
        M_shuffle(:,:,perm) = matrixPerm;
    end
    
    [clusters_sup, p_values_sup, ~, ~] = permutest(M, M_shuffle, 0, 0.05, 1000, 0);

    cluster_above_chance = zeros(timeDim,timeDim);
    for c=1:length(clusters_sup)
        if p_values_sup(c)<0.01
            cluster_above_chance(clusters_sup{c})=1;
        end
    end
    cluster_above_chance = reshape(logical(cluster_above_chance),[timeDim,timeDim]);
    
end