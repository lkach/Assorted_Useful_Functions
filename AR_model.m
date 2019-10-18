% [Phi, varNoise, S, f_vec] = AR_model(D, p)
% 
% This function uses the Yule-Walker equations to solve for the
% coefficients of an autoregressive model. Additionally, it gives a
% spectral estimate for the process based on the AR model.
% 
% Gives the coefficients (phi_1, phi_2, ... phi_p in Phi) of the
% Autoregressive model of order p of the data  vector D. Written such that
% data with gaps may be used. Typically p << length(D).
% 
% % IN:     D = Data vector (may have gaps, expressed as NaN)
% % IN:     p = order of the AR process to be modeled (i.e. AR(p)). If you
%               want the spectral estimates to be calculated from fewer
%               coefficients than order p (maybe because you want to model
%               higher lags but not have the correspondingly higher
%               frequencies in your spectral estimate), then make this a
%               two-element vector, where the first element is the AR order
%               and the second (must be < the first) is the maximum order
%               coefficient you want contributing to the spectral estimate
%               "S", which is itself an optional output.
% 
% NOTE: "Phi" is always output. Two requested outputs means that the second
%       output is "varNoise." Finally, four requested outputs results in
%       "S" and "f_vec" being assigned. Requesting three or five+ outputs
%       will result in an error message.
% 
% % OUT:    Phi = Vector of the AR coefficients, where phi_i = Phi(i), and
%                 i = 1, ..., p
% % OUT:    varNoise = variance of the forcing noise terms (which should
%                      be of zero mean
% % OUT:    S = power spectrum of D (it seems to be scaled such that
%               sum(Smean)*f_vec(1) = var(Data) approximately, which is
%               typical (here LH is actually slightly less than RHS))
% % OUT:    f_vec = vector of linear frequencies corresponding to elements
%                   in S. Units of cycles per [time unit of dt], where dt =
%                   time step on elements in D (assumed to be evenly
%                   spaced).
% 
% NOTE: If sum(Phi) = 1 (or approximately so), then the coefficients are
% close to their correct values but always a little off the mark.

% Luke Kachelein, 2019

function varargout = AR_model(D, p)

if iscolumn(D)
elseif isrow(D)
    D = D';
else
    error('D must be a vector (column or row).')
end

if length(p) == 1
    p_spectral = p;
elseif length(p) == 2 && p(1)>=p(2)
    p_spectral = p(2);
    p = p(1);
else
    error('Second argument "p" must be either a scalar or a two-element vector, with the first element >+ the second.')
end

D = D - nanmean(D); % de-mean D
N = length(D);
N_nonan = sum(isfinite(D));

% Calculate autocorrelation coefficients:
% First define the variance, which is needed for normalizing:

            varD = nansum(D.^2)/N_nonan; % for agreement
            
            % Now define the vector r, which is p+1 elements long:
            r = nan(p,1);
            for i = 2:(p+1)
                DD_lag = D(i:end).*D(1:(end - i + 1));
                r(i-1) = nanmean(DD_lag)/varD;
            end
            r0 = nanmean(D.*D)/varD;

% Error message for if r(1) ~= 1
if r0 == 1
else
    warning(['Something might be wrong; the 0th autocorrelation coefficient should be 1, but it''s actually ',r0,'.'])
end

% Now define the matrix R, which is inverted to solve for the AR
% coefficients:
R = zeros(p,p);
rr = [flip(r(1:(end-1))); r0; r(1:(end-1))];
                         % Backwards r (except r(1)) appended to r, used to
                         % build R. Note that I shaved off the last element
                         % in r because it's not used in building R.

rr_columns = (rr*rr'*ones(2*p-1))'/sum(rr);
R = spdiags(rr_columns,-(p-1):(p-1),R) + R;
                      % adding the R at the end while it's still all zeros
                      % is to revert R back to non-sparse, which is less
                      % memory intensive for non-sparse matrices like this
                      % (I really just wish that "diag" had the same
                      % functionality that "spdiags" does.

Phi = R\r; % alternatively (but maybe more slowly) "Phi = inv(R)*r;"

if nargout == 1
    varargout{1} = Phi;
elseif nargout == 2

    % Define the variance of the noise (random forcing) term, which is
    % solved on page 204 (Ch. 10.3.2 of the book cited in the "Spectral
    % Estimation" section below):
    varNoise = varD*(1 - sum(Phi.*r));
    
    varargout = cell(1,nargout);
    varargout{1} = Phi;
    varargout{2} = varNoise;
elseif nargout == 4
%% Spectral Estimation

% Reference:
% Statistical Analysis in Climate Research
% Hans von Storch, Francis W. Zwiers
% Chapter IV.12.3 (p.279)

% "auto-regressive spectral estimation"

% Define "varNoise" for this condition:
varNoise = varD*(1 - sum(Phi.*r));

df = 1/N;
fNy = 0.5;
f_vec = (df:df:fNy)';

foo = zeros(length(f_vec),p_spectral); % initialize
for k = 1:p_spectral
    foo(:,k) = Phi(k)*exp(-2*pi*1i*k*f_vec);
end
S = varNoise./(abs(1 - sum(foo,2)).^2);

% I have found that I need to multiply this by 2 to make it agree with the
% power spectrum calculated by other means (e.g. "Welch's method", SIO's
% go-to). This is probably because this (above) is the theoretical power
% spectrum, whereas the one that we get is double of half the fft (so we
% always take half the fft and double it to preserve variance in practive
% due to the ft's symmetry for real data). Because we're looking at half of
% the "possible" frequencies, we should double this S:
S = 2*S;

    varargout{1} = Phi;
    varargout{2} = varNoise;
    varargout{3} = S;
    varargout{4} = f_vec;
else
    error('Only 1, 2, or 4 outputs may be requested. Please read the ducumentation.')
end

end
