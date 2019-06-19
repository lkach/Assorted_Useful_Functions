% [Cov, Lags, Cov_std, n_Cov] = cov_gaps(A,B,maxlags)
% 
% Gives the true cross covariance between two equally sized vectors. This
% does not use any tricks to make the calculation fast, but rather it uses
% the definition of covariance:
% 
% C_ab(dt) = <A(t)B(t + dt)>
% 
% Where both A and B are de-meaned vectors.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% IN:   A =         Nx1 vector (may be 1xN, but outputs will always be columns).
%                   This may have NaN or Inf for undefined values, but what is key
%                   is that the data (whether finite or not) be even spaced.
% IN:   B =         Vector of the same size as A. If A = B, then this calculates
%                   autocovariance.
% IN:   maxlags =   (Optional, scalar) Maximum lag considered for calculating Cov.
%                   The default value is length(A) - 1.
%                   Note that Lags(end) = maxlags.
% 
% OUT:  Cov =       maxlags x 1 vector, covariance between A and B. See the
%                   equation for C_ab above.
% OUT:  Lags =      (Optional, maxlags x 1) The corresponding lags of the
%                   entries in Cov.
% OUT:  Cov_std =   (Optional, maxlags x 1) The STD of the values that went
%                   into calculating each entry in Cov.
% OUT:  n_Cov =     (Optional, maxlags x 2) The first column is the number of
%                   pairs A(t)B(t + dt) for each dt which were finite (i.e.
%                   neither A(t) nor B(t + dt) were NaN or Inf). The second
%                   column is the number of pairs that would have been used if
%                   there were no NaN's or Inf's.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [Cov, varargout] = cov_gaps(A,B,varargin)

if size(A) ~= size(B)
    error('A and B must be the same size.')
else
end

if isrow(A)
    A = A'; B = B';
else
end

N = length(A);

A = A - nanmean(A);
B = B - nanmean(B);

if nargin == 2
    maxlags = N - 1;
elseif nargin == 3
    maxlags = varargin{1};
else
    error('Incorrect number of inputs; please see documentation.')
end

Cov = zeros(1 + maxlags,1); % initialize the crosscovariance
Lags = zeros(1 + maxlags,1); % initialize the corresponding lags
Cov_std = zeros(1 + maxlags,1); % the std of the points that went into each element in pre-shaped "Cov" (starting with zero-lag)
n_Cov = zeros(1 + maxlags,2); % how many pairs were finite and how many pairs could have been finite if perfect
for i=1:(1 + maxlags) % lag = (i-1)*dt
    Crosscov_ABi = zeros(N - i + 1,1);
    Crosscov_BAi = zeros(N - i + 1,1);
    for j=1:(N - i + 1)
        Crosscov_ABi(j) = A(j)*B(j+i-1);
        Crosscov_BAi(j) = B(j)*A(j+i-1);
    end
    Crosscov_i = [Crosscov_ABi ; Crosscov_BAi];
    Lags(i) = i - 1;
    Cov(i) = nanmean(Crosscov_i); % autocov
    Cov_std(i) = nanstd(Crosscov_i); % std of points that made each element of autocov
    n_Cov(i,1) = sum(isfinite(Crosscov_i)); % how many pairs were finite
    n_Cov(i,2) = 2*(N - i + 1); % how many pairs could have been finite if perfect
end

% % Use "Cov_unfolded" (others are for completeness/clarity) for
% % estimating (cross)spectra by taking the fft.
% Cov_unfolded = [Cov;flip(Cov(2:(end-1)))];
% Cov_std_unfolded = [Cov_std;flip(Cov_std(2:(end-1)))];
% Lags_unfolded = [Lags;flip(Lags(2:(end-1)))];

if nargout == 2
    varargout{1} = Lags;
elseif nargout == 3
    varargout{1} = Lags;
    varargout{2} = Cov_std;
elseif nargout == 4
    varargout{1} = Lags;
    varargout{2} = Cov_std;
    varargout{3} = n_Cov;
elseif nargout == 0
else
    error('The number of outputs must be 1, 2, 3, or 4.')
end

% % This appears to taper just the right amount, at least for white noise:
% foo = randn(1000,1);
% [Cov, Lags, Cov_std, n_Cov] = cov_gaps(foo,foo);
% plot(Lags,Cov.*(n_Cov(:,1)/n_Cov(1,1)).^.5,'.:')

end
