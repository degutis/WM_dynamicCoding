function output = dynamicism(inputM)

timeDim = size(inputM,1);

output = zeros(timeDim,1);

for t=1:timeDim    
    output(t) = (1/(2*timeDim))*(sum(inputM(:,t))+sum(inputM(t,:)));
end

end