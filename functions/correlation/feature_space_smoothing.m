function new_data = feature_space_smoothing(data, labels, fwhm)
% function new_data = feature_space_smoothing(data, labels, fwhm)
% 
% Moves over a continuous feature space and computes the average of the
% signal of neighboring features, where neighboring features are weighted
% using a Gaussian kernel. 
%
% INPUT:
%   data: up to 4D array, where the first dim are the features and the
%   second dim is the signal. Dims 3-4 can be e.g. runs or TRs
%   labels: FxR array, where F is features and R is runs
%   fwhm: full-width half maximum (in degrees) to calculate the Gaussian
%   Kernel
%
% OUTPUT:
%   new_data: up to 4D array with processed data
%
% Simon Weber, BCCN Berlin, 15.01.2021
% edit 13.08.2021, SW: deal with NaNs in labels, labels do not have to be
% sorted
%
% TO DO:
% - add sorting of data and lables (if they are unsorted), so that unsorted
% data and labels can be used as input (output would be sorted)

% check if fwhm is positive
if fwhm < 0
    error('''fwhm'' is negative, please use positive value.');
elseif fwhm == 0
    warning('''fwhm'' is zero, no feature-space smoothing applied, continue.');
    new_data = data;
    return;
end

% fprintf('Running feature space smoothing with FWHM of %d ...\n', fwhm);

x = 1:360;
x = x-numel(x)/2;

sigma = fwhm/(2*sqrt(2*log(2)));            % calculate sigma from fwhm

y = normpdf(x,0,sigma);                     % compute pdf values at x
y = y.*1/max(y);                            % scale to 1
% figure;plot(x,y)
y = circshift(y,numel(y)/2+1);          % shift so index 1 is central 1

sz = size(data,1);

% preallocate new_data
new_data = nan(size(data));                % initialize with NaNs in case of NaN labels

for d4 = 1:size(data,4)         % in case data has up to 4 dimensions (run & TRs)
for d3 = 1:size(data,3)         % loop through runs
    
    curr_labels = round(labels(:,d3));
    curr_labels(curr_labels==0) = 1;            % account for zeros
    
    % take only valid labels (exclude NaNs)
    valid = find(~isnan(curr_labels));

    for i = 1:sz
        
        if isnan(curr_labels(i))
            continue;
        end
        
        curr_y = circshift(y, curr_labels(i)-1);        % shift to current index
        
        w = curr_y(curr_labels(valid));                          % take weightings
%         w = curr_y(curr_labels); %%%

        % scale weightings so that sum of weightings is always 1
        w = w./sum(w);

        % assign averaged, weighted data 
        new_data(i,:,d3,d4) = nansum(data(valid,:,d3,d4).*w');
    end

end
end