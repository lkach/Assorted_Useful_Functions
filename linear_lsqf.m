% [X,H,Y_fit] = linear_lsqf(Y,t,BASES)
% 
% Linear least squares. Uses the simple H\X operator. Designed for ease of
% entry for entering basis functions.
% 
% IN:
%   Y = data column vector
%   t = independent variables (e.g. time), same shape as Y
%   BASES = cell array, where each entry is a character string of the
%       desired basis functions in terms of t, e.g. 'sin(2*pi*t/12)'. The
%       exception to this formatting rule is the constant basis, which may
%       be written as '1', which is converted to 'ones(size(t))'.
%   Plot = (optional) if given, then a plot with Y and Y_fit will be
%       provided
% 
% OUT:
%   X = fit parameters, column vector of length = length(BASES)
%   H = basis matrix, also called the regressor matrix
%   Y_fit = modeled data = H*X

%%
function varargout = linear_lsqf(Y,t,BASES,varargin)

% Quality of life change: let '1' be shorthand for 'ones(size(t))'
for ii = 1:length(BASES)
    if strcmp(BASES{ii},'1')
        BASES{ii} = 'ones(size(t))';
    else
    end
end

if isrow(Y); Y = Y'; else; end
if isrow(t); t = t'; else; end
if length(Y)==length(t); else; error('dependent and independent variable vectors must be the same length'); end

H_ = nan(length(t),length(BASES));
for ii=1:length(BASES)
    eval(['H_(:,ii) = ',BASES{ii},';']);
end

t_nans = t; Y_nans = Y;
t = t(isfinite(Y)); Y = Y(isfinite(Y));

H = nan(length(t),length(BASES));

for ii=1:length(BASES)
    eval(['H(:,ii) = ',BASES{ii},';']);
end

X = H\Y;

Y_fit = H_*X;

t = t_nans;
Y = Y_nans;

if nargout == 1
    varargout{1} = X;
elseif nargout == 2
    varargout{1} = X;
    varargout{2} = H;
elseif nargout == 3
    varargout{1} = X;
    varargout{2} = H;
    varargout{3} = Y_fit;
end

if nargin > 3
    figure('color',[1 1 1])
    plot(t,Y,'b.-');hold on
    plot(t,Y_fit,'r.-')
    legend('data','fit')
else
end

end

%% Example

% DATA = 0.1*randn(1000,1) + sin(2*pi*[0:999]'/90 + 2*pi*rand); DATA(100:150) = nan;
% linear_lsqf(DATA, [0:999], {'1','sin(2*pi*t/90)','cos(2*pi*t/90)'},'');
