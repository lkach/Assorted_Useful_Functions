% [OUT N_data] = area_average_regrid_angular(IN, X_IN, Y_IN, X_OUT, Y_OUT, binDX, binDY)
% 
% Takes an input "IN" with corresponding grid points "X_IN" and "Y_IN" and
% given a new grid "X_OUT" and "Y_OUT" and bin size (see below), and gives
% "OUT", which is "IN" gridded from the appropriate area. This is meant
% only for angular values, i.e. ones which wrap around some interval. E.g.
% for angles in [0 360], the angular mean of [3 359] should be 1, not 181.
% This function reflects that. Use "area_average_regrid" for non-angular
% quantities.
% 
% INT   = interval, two-element vector with min. and max. of what is
%           considered the value at which the angle wraps back
%           (e.g. [0 2*pi])
% IN    = data, NxM
% X_IN  = x-grid on which "IN" is defined, NxM
% Y_IN  = y-grid on which "IN" is defined, NxM
% X_OUT = x-grid on which "OUT" is defined, nxm
% Y_OUT = y-grid on which "OUT" is defined, nxm
% binDX = (OPTIONAL) one half the x-side of the search box (bin) for
%         averaging to build "IN" (default 0.5*mean-dx-of-X_OUT), units of
%         X's and Y's, nxm
% binDY = (OPTIONAL) one half the y-side of the search box (bin) for
%         averaging to build "IN" (default 0.5*mean-dy-of-Y_OUT), units of
%         X's and Y's, nxm
% 
% OUT   = The gridded product after area averaging, nxm
% N_data= (OPTIONAL) the number of data in "IN" that were used to calculate
%         the corresponding value in "OUT", size nxm
% 
% If only binDX is given and not binDY, then bins are assumed to be square
% and thus binDY = binDX.

function [OUT, varargout] = area_average_regrid_angular(INT, IN, X_IN, Y_IN, X_OUT, Y_OUT, varargin)
%% Display error nessages if there are incorrectly sized inputs 
if sum(size(X_IN)==size(Y_IN)) == 2
else
    error('X_IN and Y_IN must be the same size.')
end

if sum(size(X_OUT)==size(Y_OUT)) == 2
else
    error('X_OUT and Y_OUT must be the same size.')
end

if nargin == 6
    binDX = abs(0.5*nanmean(nanmean(diff(X_OUT,1,2)))*ones(size(X_OUT)));
    binDY = abs(0.5*nanmean(nanmean(diff(Y_OUT,1,1)))*ones(size(Y_OUT)));
elseif nargin == 7
    binDX = varargin{1};
    binDY = binDX;
elseif nargin == 8
    binDX = varargin{1};
    binDY = varargin{2};
else
    error('Only 5, 6, or 7 inputs are accepted.')
    help area_average_regrid
end

if (sum(size(X_OUT)==size(binDX)) == 2) && (sum(size(binDX)==size(binDY)) == 2)
else
    error('binDX and binDY must be the same size as X_OUT and Y_OUT.')
end

%% Average and interpolate

OUT = nan(size(X_OUT));
for yi = 1:size(X_OUT,1)
    for xi = 1:size(X_OUT,2)
        % Note that the < vs >= is arbitrary
        sampled_points = IN(Y_IN < Y_OUT(yi,xi) + binDY(yi,xi) & Y_IN >= Y_OUT(yi,xi) - binDY(yi,xi) & ...
                            X_IN < X_OUT(yi,xi) + binDX(yi,xi) & X_IN >= X_OUT(yi,xi) - binDX(yi,xi));
        sampled_points = reshape(sampled_points,numel(sampled_points),1);
        nanmean_ij = Angular_Mean(sampled_points,INT);
        OUT(yi,xi) = nanmean_ij;
        N_data(yi,xi) = sum(isfinite(sampled_points));
    end
end

if nargout == 2
    varargout{1} = N_data;
else
end

end

%% Angular mean function
function M = Angular_Mean(D,I)

if isrow(D)
    D = D';
else
end

if sum(D > max(I) | D < min(I))
    error('All values in D must be in the interval I.')
else
end

dI = max(I) - min(I);
D_ = D - min(I);
D_ = D_*(2*pi/dI) - pi;
D_ = exp(1i*D_);

M = nanmean(D_);
M = angle(M) + pi;
M = M*(dI/(2*pi));
M = M + min(I);

end

%% Example:

% % % Compare the following:

% area_average_regrid(repmat([3 3 3 3 359 359 359 359],5,1),...
%                     [1:8;1:8;1:8;1:8;1:8],[1:5;1:5;1:5;1:5;1:5;1:5;1:5;1:5]',...
%                     [1:7;1:7;1:7;1:7;1:7]+0.5,[1:5;1:5;1:5;1:5;1:5;1:5;1:5]'+0.5,...
%                     ones(5,7),ones(5,7))
% 
% area_average_regrid_angular([0,360],repmat([3 3 3 3 359 359 359 359],5,1),...
%                             [1:8;1:8;1:8;1:8;1:8],[1:5;1:5;1:5;1:5;1:5;1:5;1:5;1:5]',...
%                             [1:7;1:7;1:7;1:7;1:7] + 0.5,[1:5;1:5;1:5;1:5;1:5;1:5;1:5]' + 0.5,...
%                             ones(5,7),ones(5,7))

% % % Which respectively output different results (the latter being
% % % desired, because the angular mean of [3 359] should be 1, not 181):

%      3     3     3   181   359   359   359
%      3     3     3   181   359   359   359
%      3     3     3   181   359   359   359
%      3     3     3   181   359   359   359
%      3     3     3   181   359   359   359
% 
%      3     3     3     1   359   359   359
%      3     3     3     1   359   359   359
%      3     3     3     1   359   359   359
%      3     3     3     1   359   359   359
%      3     3     3     1   359   359   359
