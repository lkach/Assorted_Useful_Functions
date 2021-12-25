% Gives the angular mean of a vector of values, all of which must be in the
% interval of interest.
% 
% IN:   D = data, vector; all elements in D must lie in I
% IN:   I = interval, two-element vector with min. and max. of what is
%           considered the value at which the angle wraps back
%           (e.g. [0 2*pi])
% 
% OUT:  M = angular mean of D
% OUT:  C = (optional) centered version of D, i.e. with angular mean of
%           mean(I). This will have the same size as D

function [M,varargout] = angular_mean(D,I)

if isrow(D)
    D = D';
    transposeC = 1;
else
    transposeC = 0;
end

if sum(D > max(I) | D < min(I))
    error('All values in D must be in the interval I.')
else
end

dI = max(I) - min(I);
D_ = D - min(I);
D_ = D*(2*pi/dI);
D_ = exp(1i*D_);

M = nanmean(D_);
M = angle(M);
M = M*(dI/(2*pi));
M = M + min(I);

if nargout == 2
    % lkaddpath
    % warning('Requires custom function "interval_wrap".')
    C = interval_wrap(D - M,[min(I) max(I)] - mean(I));
    if transposeC
        C = C';
    else
    end
    varargout{1} = C;
elseif nargout > 2
    error('Too many outputs.')
else
end

end
