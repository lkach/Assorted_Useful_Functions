% [OUT N_data] = area_average_regrid(IN, X_IN, Y_IN, X_OUT, Y_OUT, binDX, binDY)
% 
% Takes an input "IN" with corresponding grid points "X_IN" and "Y_IN" and
% given a new grid "X_OUT" and "Y_OUT" and bin size (see below), and gives
% "OUT", which is "IN" gridded from the appropriate area.
% 
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

function [OUT, varargout] = area_average_regrid(IN, X_IN, Y_IN, X_OUT, Y_OUT, varargin)
%% Display error nessages if there are incorrectly sized inputs 
if sum(size(X_IN)==size(Y_IN)) == 2
else
    error('X_IN and Y_IN must be the same size.')
end

if sum(size(X_OUT)==size(Y_OUT)) == 2
else
    error('X_OUT and Y_OUT must be the same size.')
end

if nargin == 5
    binDX = abs(0.5*nanmean(nanmean(diff(X_OUT,1,2)))*ones(size(X_OUT)));
    binDY = abs(0.5*nanmean(nanmean(diff(Y_OUT,1,1)))*ones(size(Y_OUT)));
elseif nargin == 6
    binDX = varargin{1};
    binDY = binDX;
elseif nargin == 7
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
        nanmean_ij = nanmean(sampled_points);
        OUT(yi,xi) = nanmean_ij;
        N_data(yi,xi) = sum(isfinite(sampled_points));
    end
end

if nargout == 2
    varargout{1} = N_data;
else
end

end
%% Example

% xx = -10:0.1:10;
% yy = xx;
% sx = 3;
% sy = 1;
% 
% % IN = ZZ
% [XI,YI] = meshgrid(xx,yy);
% ZI = exp(-(XI.^2)/(2*sx^2) - (YI.^2)/(2*sy^2)) + ...
%      0.5*exp(-((XI - 3).^2)/(2*sx^2) - ((YI - 3).^2)/(2*sy^2)) + ...
%      2*exp(-((XI - 0).^2)/(.1*2*sx^2) - ((YI + 7).^2)/(2*2*sy^2));
% ZI(100:105,100:105) = NaN;
% 
% [XO,YO] = meshgrid([-9:9] + 0.5, [-9:9] + 0.5);
% % OUT = ZO
% [ZO,NO] = area_average_regrid(ZI, XI, YI, XO, YO);
% 
% figure
% subplot(121)
% contourf(XI,YI,ZI)
% subplot(122)
% contourf(XO,YO,ZO)
% 
% figure
% subplot(131)
% imagesc(ZI)
% subplot(132)
% imagesc(ZO)
% subplot(133)
% imagesc(NO)

