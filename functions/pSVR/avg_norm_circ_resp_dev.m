function [Acc, Accs] = avg_norm_circ_resp_dev(xr, x)
% _
% Average Normalized Circular Response Deviation
% FORMAT Acc = avg_norm_circ_resp_dev(xr, x)
% 
%     xr   - an n x v matrix of reconstructions (n observations, v channels)
%     x    - an n x 1 vector of angles, in radians, i.e. 0 <= x_i <= 2 pi
% 
%     Acc  - a  1 x v vector of avg. norm. circ. resp. dev. accuracy values
%     Accs - an n x v matrix of trial-wise (not averaged) accuracy values
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 26/08/2019, 16:45
%  Last edit: 22/01/2020, 14:19


% quantify accuracy
%-------------------------------------------------------------------------%
dx      = xr - repmat(x,[1 size(xr,2)]);
dxi     = abs(dx) > pi;
dx(dxi) = -1*sign(dx(dxi)) * 2*pi + dx(dxi);

% Compute accuracy
%-------------------------------------------------------------------------%
Accs = (pi - abs(dx))./pi;
Acc  = nanmean(Accs, 1);