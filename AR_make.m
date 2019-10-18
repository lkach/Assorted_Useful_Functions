% Out = AR_make(phi,N)
% 
% Given AR(p) coefficients phi (p-long vector), get a random autoregressive
% process of order p, N elements long.
% 
% IN:   phi = AR(p) coefficients, where phi(1) is the first lag coefficient
%             and phi(n) = n'th lag coefficient. Length p.
% IN:   N   = length of desired time series built as an AR(p) process with
%             AR coefficients phi.
% 
% OUT:  Out = length N AR(p) process with AR coefficients phi. This will
%             have a variance of 1 and mean of 0.

function Out = AR_make(phi,N)

if isrow(phi)
    phi = phi';
else
end

p = length(phi);

rand_white = randn(N + p, 1);
rand_white = rand_white/std(rand_white); % rand_white should already have
                                         % var = 1, but this forces it.

Out = rand_white;

for i = (p+1):length(Out)
    Out(i) = sum( Out((i-p):(i-1)).*flip(phi) ) + rand_white(i);
end

Out = Out((p+1):end);

end
