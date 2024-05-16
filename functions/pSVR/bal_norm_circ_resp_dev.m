function [BPrec, Bw] = bal_norm_circ_resp_dev(xr, x, meth)
% _
% Balanced Normalized Circular Response Deviation
% FORMAT [BPrec, Bw] = bal_norm_circ_resp_dev(xr, x, meth)
% 
%     xr    - an n x v matrix of reconstructions (n observations, v channels)
%     x     - an n x 1 vector of angles, in radians, i.e. 0 <= x_i <= 2 pi
%     meth  - method to use for numerical integration ('trapz' or 'rect')
% 
%     BPrec - a  1 x v vector of bal. norm. circ. resp. dev. precision values
%     Bw    - an n x v matrix of weights used to balance trial-wise precisions
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 21/01/2020, 17:48
%  Last edit: 04/03/2020, 12:26


% specify method
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(meth), meth = 'rect'; end;

% quantify precision
%-------------------------------------------------------------------------%
xn      = ~isnan(x);
dx      = xr(xn,:) - repmat(x(xn),[1 size(xr,2)]);
dxi     = abs(dx) > pi;
dx(dxi) = -1*sign(dx(dxi)) * 2*pi + dx(dxi);

% compute precisions
%-------------------------------------------------------------------------%
BPrecs = (pi - abs(dx))./pi;

% sort precisions
%-------------------------------------------------------------------------%
[xs, is] = sort(x(xn));
BPrecs   = BPrecs(is,:);

% integrate numerically
%-------------------------------------------------------------------------%
if strcmp(meth,'trapz')
    xs     = [xs; (2*pi)+min(xs)];
    BPrecs = [BPrecs; BPrecs(1,:)];
    BPrec  = 1/(2*pi) * trapz(xs, BPrecs, 1);
end;
if strcmp(meth,'rect')
    xs    = [max(xs)-(2*pi); xs; (2*pi)+min(xs)];
    ds    = diff(xs);
    Bw    = (ds(1:end-1)/2 + ds(2:end)/2)./(2*pi);
    BPrec = sum(repmat(Bw,[1 size(xr,2)]).*BPrecs);
    Bw(is)= Bw;
end;