% cov_gaps_window
% 
% Not to be confused with MATLAB's default "nancov", which does basically
% the same thing. This (cov_gaps_window) exists because I know exactly what went
% into it.
% 
% Possible inputs:
% 
% IN:   X      = Nx1 or 1xN vector of data, evenly spaced and may have
%                NaNfor missing data.
% IN:   Y      = Same size as X, this is for cases where you want the
%                covariance between two different vectors, not an
%                autocovariance of X with itself.
% IN:   maxlag = The maximum number of lags considered. The deafult is
%                N - 1. For very long X and Y, this needs to be much
%                smaller than N - 1, otherwise it could take a long time to
%                do the calculation (a warning will be issued so that the
%                user knows if they need to force stop).
% IN:   Window = The window by which the raw covariance is multiplied; this
%                is useful for suppressing spurious correlation at large
%                lags for which there are few pairs whose products are
%                averaged. Default value is 'rectwin'.
%                This input may be a string (e.g. 'hanning', 'triang',
%                etc.)
% 
% Possible outputs:
% 
% OUT:  cov_xy     = The covariance between X and Y (or covariance of X).
%                    Length = maxlag + 1
% OUT:  lags       = The corresponding lags. Length = maxlag + 1
% OUT:  cov_std_xy = Standard deviation of points that made each element of
%                    autocov. Length = maxlag + 1
% OUT:  n_cov_xy   = 2-column matrix (Length = maxlag + 1):
%                    - First is how many pairs were finite for calculating
%                      the corresponding element in cov_xy
%                    - Second is how many pairs would have been used if
%                      there were no gaps (NaN's).
%                    This is given so that n_cov_xy(:,1)./n_cov_xy(:,2) is
%                    another metric for how good your estimate of cov_xy is
%                    given that X and/or Y have gaps.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                       POSSIBLE INPUT ARRANGEMENTS:
% {X}
% 
% {X, Y}
% {X, maxlag}
% 
% {X, maxlag, Window}
% {X, Y,      maxlag}
% {X, Y,      Window}
% 
% {X, Y, maxlag, Window}
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                       POSSIBLE OUTPUT ARRANGEMENTS:
% {cov_xy}
% 
% {cov_xy, lags}
% 
% {cov_xy, lags, cov_std_xy}
% 
% {cov_xy, lags, cov_std_xy, n_cov_xy}
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


function varargout = cov_gaps_window(varargin)%(varargin):

if nargin == 1 % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    X = varargin{1};
    Y = varargin{1};
    maxlag = length(X) - 1;
    Window = 'rectwin';
elseif nargin == 2 % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    X = varargin{1};
    if length(varargin{2}) == 1
        Y = varargin{1};
        maxlag = varargin{2};
        Window = 'rectwin';
    elseif size(varargin{2}) == size(varargin{1})
        Y = varargin{2};
        maxlag = length(X) - 1;
        Window = 'rectwin';
    else
        error('Inputs are incorrectly formatted.')
    end
elseif nargin == 3 % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    X = varargin{1};
    if length(varargin{2}) == 1
        Y = varargin{1};
        maxlag = varargin{2};
        Window = varargin{3};
    elseif size(varargin{2}) == size(varargin{1})
        Y = varargin{2};
        if length(varargin{3}) == 1
            maxlag = varargin{3};
            Window = 'rectwin';
        elseif ischar(varargin{3})
            maxlag = length(X) - 1;
            Window = varargin{3};
        else
            error('Inputs are incorrectly formatted.')
        end
    else
        error('Inputs are incorrectly formatted.')
    end
elseif nargin == 4 % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    X = varargin{1};
    Y = varargin{2};
    maxlag = varargin{3};
    Window = varargin{4};
else
    error('Inputs are incorrectly formatted.')
end % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if isrow(X)
    X = X';
    Y = Y';
elseif iscolumn(X)
else
    error('Inputs are incorrectly formatted.')
end

X(~isfinite(X)) = nan;
Y(~isfinite(Y)) = nan;

N_x = length(X);
X = X - nanmean(X);
Y = Y - nanmean(Y);

cov_xy = zeros(1 + maxlag,1); % the autocovariance
cov_std_xy = zeros(1 + maxlag,1); % the std of the points that went into each element in pre-shaped "autocov_x" (starting with zero-lag)
n_cov_xy = zeros(1 + maxlag,2); % how many pairs were finite and how many pairs could have been finite if perfect
for i=1:(1 + maxlag) % lag = (i-1)*dt
    cov_xy_i = zeros(N_x - i + 1,1);
    for j=1:(N_x - i + 1)
        cov_xy_i(j) = X(j)*Y(j+i-1);
    end
    cov_xy(i) = nanmean(cov_xy_i); % autocov
    cov_std_xy(i) = nanstd(cov_xy_i); % std of points that made each element of autocov
    n_cov_xy(i,1) = sum(isfinite(cov_xy_i)); % how many pairs were finite
    n_cov_xy(i,2) = N_x - i + 1; % how many pairs could have been finite if perfect
end
% NOTE: "cov_std_xy" and "n_cov_xy" are unused, and are
% retained here for troubleshooting and possible inclusion.

full_cov_xy = [cov_xy;flip(cov_xy(2:(end-1)))];
% ^ Used in "eval"; this is also what one would use to calculate the
% spectrum with the FFT, but seeing as I already have a function
% (WKT_spectrum) that does that, I do not have "full_cov_xy" as an output.
% Feel free to modify if desired.
if strcmp(Window,'blackman') || ...
        strcmp(Window,'flattopwin') || ...
        strcmp(Window,'hamming') || ...
        strcmp(Window,'hann') || ...
        strcmp(Window,'hanning') || ...
        strcmp(Window,'blackmanharris')
    Window_ = eval(['fftshift(',Window,'(length(full_cov_xy),''periodic''));']);
else
    Window_ = eval(['fftshift(',Window,'(length(full_cov_xy)));']);
end

cov_xy = Window_(1:length(cov_xy)).*cov_xy;
lags = 0:(maxlag);

if nargout == 0
    ans = cov_xy
elseif nargout == 1
    varargout{1} = cov_xy;
elseif nargout == 2
    varargout{1} = cov_xy;
    varargout{2} = lags;
elseif nargout == 3
    varargout{1} = cov_xy;
    varargout{2} = lags;
    varargout{3} = cov_std_xy;
elseif nargout == 4
    varargout{1} = cov_xy;
    varargout{2} = lags;
    varargout{3} = cov_std_xy;
    varargout{4} = n_cov_xy;
else
    error('Only 1, 2, 3, or 4 outputs are expected.')
end

end
