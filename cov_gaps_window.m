% COV_GAPS_WINDOW
% 
% Possible inputs:
% 
% IN:   X      = Nx1 or 1xN vector of data, evenly spaced and may have
%                NaNfor missing data.
% IN:   Y      = Same size as X, this is for cases where you want the
%                covariance between two different vectors, not an
%                autocovariance of X with itself.
% IN:   maxlag = The maximum number of lags considered. The default is
%                +-(N - 1). For very long X and Y, this needs to be much
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
%                    Length = 2*maxlag + 1
% OUT:  lags       = The corresponding lags. Length = 2*maxlag + 1
% OUT:  cov_std_xy = Standard deviation of points that made each element of
%                    autocov. Length = 2*maxlag + 1
% OUT:  n_cov_xy   = 2-column matrix (Length = 2*maxlag + 1):
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
    elseif ischar(varargin{2})
        Y = varargin{1};
        maxlag = length(X) - 1;
        Window = varargin{2};
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
    Y = conj(Y); % for when the inputs are complex

cov_xy = zeros(1 + 2*maxlag,1); % the autocovariance
cov_std_xy = zeros(1 + 2*maxlag,1); % the std of the points that went into each element in pre-shaped "autocov_x" (starting with zero-lag)
n_cov_xy = zeros(1 + 2*maxlag,2); % how many pairs were finite and how many pairs could have been finite if perfect
% There is certainly a way to code this to go in one loop, but this works:
for i=1:(1 + maxlag) % lag = (i-1)*dt
    cov_xy_i = X(i:end).*Y(1:(end-i+1)); % length N_x - i + 1 (I think)
    cov_xy(i+maxlag) = nanmean(cov_xy_i); % autocov
    cov_std_xy(i+maxlag) = nanstd(cov_xy_i); % std of points that made each element of autocov
    n_cov_xy(i+maxlag,1) = sum(isfinite(cov_xy_i)); % how many pairs were finite
    n_cov_xy(i+maxlag,2) = N_x - i + 1; % how many pairs could have been finite if perfect
end
for i=2:(1 + maxlag) % lag = -(i-1)*dt
    cov_xy_i = X(1:(end-i+1)).*Y(i:end);
    cov_xy(maxlag+2-i) = nanmean(cov_xy_i);
    cov_std_xy(maxlag+2-i) = nanstd(cov_xy_i);
    n_cov_xy(maxlag+2-i,1) = sum(isfinite(cov_xy_i));
    n_cov_xy(maxlag+2-i,2) = N_x - i + 1;
end
% NOTE: "cov_std_xy" and "n_cov_xy" are unused, and are
% retained here for troubleshooting and possible inclusion later.

% "cov_xy" is also what one would use to calculate the
% spectrum with the FFT, but seeing as I already have a function
% (WKT_spectrum) that does that, I do not have "cov_xy" as an output.
% Feel free to modify if desired.
if strcmp(Window,'blackman') || ...
        strcmp(Window,'flattopwin') || ...
        strcmp(Window,'hamming') || ...
        strcmp(Window,'hann') || ...
        strcmp(Window,'hanning') || ...
        strcmp(Window,'blackmanharris')
    Window_ = eval([Window,'(length(cov_xy),''periodic'');']);
else
    Window_ = eval([Window,'(length(cov_xy));']);
end

cov_xy = Window_.*cov_xy;
lags = [(-maxlag):(maxlag)]';

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

%% Demonstration and comparison to xcorr:

% X0 = zeros(100,1);
% for i=2:100; X0(i) = X0(i-1) + 0.9*randn; end
% X1 = X0 + 0.5*randn(100,1); X1 = wshift('1D',X1,10);
% X2 = X0 + 0.3*randn(100,1);
% [C1,lags1] = cov_gaps_window(X1/std(X1),X2/std(X2),'triang');
% [C2,lags2] = xcorr(X1-mean(X1),X2-mean(X2),'coef');
% figure;plot(X1,'.-');hold on;plot(X2,'.-')
% figure;plot(lags1,C1,'.-');hold on;plot(lags2,C2,'.-')

