% Cov = slope_to_cov(Slope, length)
% 
% IN:   Slope   = spectral slope (may be positive or negative)
% IN:   Length  = number of autocovariance coefficients to be given
% 
% OUT:  AutoCov = autocovariance coefficients corresponding with the
%                 spectral slope of "Slope" via the Wiener-Khinchin theorem
% OUT:  Lags    = lags corresponding to the elements in "AutoCov"

function [AutoCov, Lags] = slope_to_autocov(Slope, Length)

df = 1/(2*Length);
f = [df:df:0.5]';
S = f.^Slope;
S_ = [S; flip(S(2:(end-1)))];
AutoCov = ifft(S_);
Lags = [0:(length(AutoCov)-1)]';

% Truncate:
AutoCov = AutoCov(1:Length);
Lags = Lags(1:Length);

% figure;plot(Lags,AutoCov,'.-');hold on;plot(Lags(1:Length),AutoCov(1:Length),'o--')

end
