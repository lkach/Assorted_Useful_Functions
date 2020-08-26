% Out = AR_make(phi,N)
% 
% Given AR(p) coefficients phi (p-long vector), get a random autoregressive
% process of order p, N elements long.
% 
% IN:   phi = AR(p) coefficients, where phi(1) is the first lag coefficient
%             and phi(n) = n'th lag coefficient. Length p.
% IN:   N   = length of desired time series built as an AR(p) process with
%             AR coefficients phi.
% IN:   Init = (Optional) Initialization time series of length p from which
%             "OUT" will be constructed based on "phi". Default is white
%             noise of variance 1.
% 
% OUT:  Out = length N AR(p) process with AR coefficients phi. This will
%             have a variance of 1 and mean of 0.

function Out = AR_make(phi,N,varargin)

p = length(phi);

if nargin == 3
    Init = varargin{1};
    if isrow(Init); Init = Init'; else; end
    if length(Init) ~= p; error('"Init" must have the same length as "phi".'); else; end
elseif nargin == 2
    Init = randn(p, 1);
else
    error('Wrong number of inputs.')
end

if isrow(phi)
    phi = phi';
else
end

rand_white = randn(N + p, 1);
rand_white = rand_white/std(rand_white); % rand_white should already have
                                         % var = 1, but this forces it.

Out = rand_white;
Out(1:p) = Init;

for i = (p+1):length(Out)
    Out(i) = sum( Out((i-p):(i-1)).*flip(phi) ) + rand_white(i);
end

Out = Out((p+1):end);

Out = Out/std(Out);
Out = Out - mean(Out);

end
