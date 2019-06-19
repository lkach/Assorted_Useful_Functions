% [C_xx, C_yy, C_xy, Lags] = plot_cov_gaps(x,y,maxlags)
% 
% Plots and gives the autocovariances and covariance of two variables.
% 
% Requires "cov_gaps.m"
% 
% IN:   x = Nx1 or 1xN vector, first variable
% IN:   y = second variable, same size as x
% IN:   maxlags = (Optional) either an integer (or blank), in which case
%           all covariances are calculated up to the same maximum lag, or
%           it is given as a 3x1 or 1x3 vector of integers, which are the
%           respective maxlags of C_xx, C_yy, and C_xy
% 
% OUT:  C_xx = Autocovariance of x
% OUT:  C_yy = Autocovariance of y
% OUT:  C_xx = Covariance of x with y

function [C_xx, C_yy, C_xy] = plot_cov_gaps(x,y,varargin)


if nargin == 2
    [C_xx, Lags_xx] = cov_gaps(x,x);
    C_yy = cov_gaps(y,y);
    C_xy = cov_gaps(x,y);
    Lags_yy = Lags_xx;
    Lags_xy = Lags_xx;
elseif nargin == 3
    if length(varargin{1}) == 1
        maxlags = varargin{1}*[1 1 1];
    elseif length(varargin{1}) == 3
        maxlags = varargin{1};
    else
        error('Please format the third argument correctly: positive integer or 3-element vector of integers.')
    end
    
    [C_xx, Lags_xx] = cov_gaps(x,x,maxlags(1));
    [C_yy, Lags_yy] = cov_gaps(y,y,maxlags(2));
    [C_xy, Lags_xy] = cov_gaps(x,y,maxlags(3));
else
    error('Incorrect number of inputs, see documentation.')
end

figure
plot(Lags_xx,C_xx,'.:');hold on
plot(Lags_yy,C_yy,'.:')
plot(Lags_xy,C_xy,'.:')
legend('C_x_x','C_y_y','C_x_y')
xlabel('Lags')
ylabel('Covariance')

end
