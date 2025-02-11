function outputMatrix = dynamicsMultiply(originalMatrix)
    
    outputMatrix = zeros(size(originalMatrix));
    
    for roi=1:size(outputMatrix,3)
        newMatrix = [];
        diagonalElements = diag(originalMatrix(:,:,roi));
        zeroDiagonalIndices = find(diagonalElements == 0);  
        newMatrix = ones(size(originalMatrix(:,:,roi)));
        newMatrix(zeroDiagonalIndices, :) = 0; 
        newMatrix(:, zeroDiagonalIndices) = 0; 
        newMatrix(1:4,:)= 0; % too weak of a signal in the start
        newMatrix(:,1:4) = 0;
        outputMatrix(:,:,roi) = newMatrix;
    end;
