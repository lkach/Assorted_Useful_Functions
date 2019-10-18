% [S, f, C, tau] = AR_Spec(AR_coef,dt,T)
% 
% Spectrum of an AR process:
% 
% S(f) = \sigma^2/|1 - \sum_{k=1}^{p}(\phi_k \exp(-i 2\pi f k) )|^2
% 
% where \sigma^2 = 1 is assumed (scale it on your own if necessary).
% 
% IN:   AR_coef     Pretty self-explanatory, p-long vector, first entry is
%                   the coefficient for lag 1*dt
% IN:   dt          Time step (used in calculating f)
% IN:   T           Length of record (used in calculating f)
% 
% OUT:  S           Power spectrum (variance units)
% OUT:  f           Corresponding frequencies (natural, not angular)
% OUT:  C           Autocovariance (ifft of S via WKTheorem)
% OUT:  tau         Lags corresponding to entries in C
% 
% Outputs after "S" are optional.

function varargout = AR_Spec(AR_coef,dt,T)

if isrow(AR_coef)
elseif iscolumn(AR_coef)
    AR_coef = AR_coef';
else
    error('"AR_coef" is formatted incorrectly, please make it a vector.')
end

p = length(AR_coef);

df = 1/T;
fNy = 1/(2*dt);

f = (df:df:fNy)';

F = repmat(f,1,p);

S = 1./abs(1 - sum(exp(-2*pi*1i*F.*[1:p]).*AR_coef,2)).^2;

if nargout == 1
    varargout{1} = S;
elseif nargout == 2
    varargout{1} = S;
    varargout{2} = f;
elseif nargout == 3
    varargout{1} = S;
    varargout{2} = f;
    SS = [S;flip(S(2:(end-1)))];
    C = ifft(SS);
    C = C/C(1);
    C = C(1:length(f));
    varargout{3} = C;
elseif nargout == 4
    varargout{1} = S;
    varargout{2} = f;
    SS = [S;flip(S(2:(end-1)))];
    C = ifft(SS);
    C = C/C(1);
    C = C(1:length(f));
    varargout{3} = C;
    tau = 0:(length(C) - 1);
    varargout{4} = tau;
else
    error('1, 2, 3, or 4 outputs are expected.')
end

end

%% Demonstration

% ARC = [.8,-.2];
% foo = randn_ar(1000,ARC); foo = foo/std(foo);
% figure;plot(foo)
% [S,f,C,tau] = AR_Spec(ARC,1,length(foo));
% figure;semilogy(f,S,'.-')
% figure;plot(real(C),'.-');hold on;plot(imag(C))
%        plot(flip(fftshift(xcorr(foo,'coeff'))),'.-r')
% figure;semilogy(abs(fft(foo)).^2 / length(foo),'.-k');hold on;semilogy(S,'.-')
