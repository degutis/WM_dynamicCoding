function cluster_above_chance = dynamicCoding_signInversion(input_matrix)

rng(abs(round(input_matrix(1,5,1))))

matrixSize = size(input_matrix);
timeDim = matrixSize(1);
test1=zeros(matrixSize);
test2=zeros(matrixSize);
PlusMin = [-1,1];
for sub=1:matrixSize(3)
    for a= 1:matrixSize(1)
        for b=1:matrixSize(2)
            test1(a,b,sub)=input_matrix(a,a,sub)-input_matrix(a,b,sub); 
            test2(a,b,sub)=input_matrix(b,b,sub)-input_matrix(a,b,sub);
        end
    end
end

matrixShuffled1 = zeros([matrixSize(1),matrixSize(2),1000]);
matrixShuffled2 = zeros([matrixSize(1),matrixSize(2),1000]);
numSubs=1:matrixSize(3);

for perm=1:1000

    matrixPerm1 = zeros(timeDim,timeDim);
    matrixPerm2 = zeros(timeDim,timeDim);

    for aa=1:timeDim
        for bb=1:timeDim
            matrixPerm1(aa,bb) = test1(aa,bb,randsample(numSubs,1))*randsample(PlusMin,1);
            matrixPerm2(aa,bb) = test2(aa,bb,randsample(numSubs,1))*randsample(PlusMin,1);
        end
    end
        
    matrixShuffled1(:,:,perm)=matrixPerm1;
    matrixShuffled2(:,:,perm)=matrixPerm2;
end

[clusters_test1, p_values_test1, ~, ~] = permutest(test1, matrixShuffled1, 0, 0.05, 1000, 0);
cluster_above_chance_test1 = zeros([matrixSize(1),matrixSize(2)]);
for c=1:length(clusters_test1)
    if p_values_test1(c)<0.01
        cluster_above_chance_test1(clusters_test1{c})=1;
    end
end
cluster_above_chance_test1 = reshape(logical(cluster_above_chance_test1),[size(input_matrix,1),size(input_matrix,2)]);


[clusters_test2, p_values_test2, ~, ~] = permutest(test2, matrixShuffled2, 0, 0.05, 1000, 0);
cluster_above_chance_test2 = zeros([matrixSize(1),matrixSize(2)]);
for c=1:length(clusters_test2)
    if p_values_test2(c)<0.01
        cluster_above_chance_test2(clusters_test2{c})=1;
    end
end
cluster_above_chance_test2 = reshape(logical(cluster_above_chance_test2),[size(input_matrix,1),size(input_matrix,2)]);


cluster_above_chance = cluster_above_chance_test1+cluster_above_chance_test2;
cluster_above_chance(cluster_above_chance<2)=0;
cluster_above_chance(cluster_above_chance==2)=1;

end