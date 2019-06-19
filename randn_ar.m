% randn_ar is like randn, but gives an autoregressive process instead of a
% white one.
% 
% X_t = \sum_{i=1}^{p}(a_i X_{t-i}) + \varepsilon_t
% 
% [X, S] = randn_ar(N, A)
% 
% IN:   N = Prescribed length of X
% IN:   A = (Optional) The AR coefficients (default A = 1). Length of A is
%           the order of the AR process (p), where A(i) is a_i
% 
% OUT:  X = Nx1 time series for this AR process.
% OUT:  S = (Optional) floor(N/2)x2 matrix: first column is frequencies,
%           second column is the true spectrum of X, calculated from:
% 
%       S(f) = 1/abs(1 - \sum_{k=1}^p (a_k exp(-i 2 pi f k)) )^2
% 
%       Formally, the spectrum should be scaled by the noise variance, but
%       this is not done here.

function varargout = randn_ar(varargin)

if nargin == 1
    N = varargin{1};
    a = 1;
    p = 1;
elseif nargin == 2
    N = varargin{1};
    a = varargin{2};
    p = length(a);
    if isrow(a); a = a'; else; end
else
    error('Only one or two inputs are expected. Please see documentation.')
end

Out = randn(N+p,1);

for i = (1+p):(N+p)
    Out(i) = Out(i) + sum( flip(a).*Out( (i-p):(i-1) ) );
end

Out = Out((1+p):end);

if nargout == 1
    varargout{1} = Out;
elseif nargout == 2
    df = 1/N;
    f = (df:df:0.5)';
    k = 1:p;
    Spec = 1./abs( (1 - sum(exp(-1i*2*pi*(f*k) ).*a', 2) ) ).^2;
    varargout{1} = Out;
    varargout{2} = [f,Spec];
else
    error('Only one or two outputs are expected. Please see documentation.')
end

end

%% Demonstration (requires "nanspectrum")

% close all
% [foo,sss] = randn_ar(10003,[.6,-.5,0,0,0,.3,-0.4]);
% figure;plot(foo)
% 
% figure;loglog(sss(:,1),sss(:,2),'.-');hold on
% nanspectrum(foo, 1, 'step', 7, '.-', 1, 0, 'hanning');
