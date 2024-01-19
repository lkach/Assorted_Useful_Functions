% [X,J,Y_fit,ConvergenceRecord] = nonlinear_lsqf(Y,t,BASIS,PARAMETERS_0,dBASES,Tol,   MaxIt,Plot)
% 
% Nonlinear least squares. Uses the simple H\X operator. Designed for ease
% of entry for entering basis functions. This is designed to not need the
% MATLAB symbolic math toolbox, which would make it easier by removing the
% requirement that "dBASES" be given.
% This function uses the Gauss?Newton algorithm.
% 
% IN:
%   Y = data column vector
%   t = independent variables (e.g. time), same shape as Y
%   BASIS = cell array where the first element is a single basis function
%       in terms of "t" and any number of parameters, e.g. "a" and "f". The
%       following elements list the names of the parameters, e.g.:
%       {'a*sin(2*pi*f*t)','a','f'}
%       Note that this is linear in "a" but not in "f".
%   PARAMETERS_0 = initial guesses for the parameters, vector of length =
%       length(BASIS)-1. Matching the order of parameters as given in BASIS
%   dBASES = cell array, where each entry is a character string of the
%       partial derivative of the desired basis functions in terms of t,
%       with respect to each parameter. E.g. if the model function is
%       'a*sin(2*pi*f*t)', where amplitude a and frequency f are the
%       parameters, then dBASES =
%       {'sin(2*pi*f*t)','2*pi*t*a*cos(2*pi*f*t)'}
%   Tol = tolerancefor convergence, defined here as the sum of the absolute
%       values of the difference between two iterations of the parameters
%       divided by the absolute value of the earlier iteration, e.g.
%       abs(a2 - a1)/abs(a1) + abs(f2 - f1)/abs(f1)
%   MaxIt = (optional) maximum number of iterations (default 100)
%   Plot = (optional) if given, then a plot with Y and Y_fit will be
%       provided
% 
% OUT:
%   X = fit parameters, column vector of length = length(BASIS)-1 = length(dBASES)
%   J = Jacobian matrix
%   Y_fit = modeled data = H*X
%   ConvergenceRecord = each value of convergence before "Tol" is reached
%       or "MaxIt" iterations are performed, whichever is sooner. Variable
%       length vector.
%   ResidualVar_Record = same size as "ConvergenceRecord", lists the values
%       of the variance of the residuals at each iteration
%   X_record = matrix of variable size containing the parameter values at
%       each iteration, including the initial guess in the first column


%%
function varargout = nonlinear_lsqf(Y,t,BASIS,PARAMETERS_0,dBASES,Tol,varargin)

if isrow(Y); Y = Y'; else; end
if isrow(t); t = t'; else; end
if isrow(PARAMETERS_0); PARAMETERS_0 = PARAMETERS_0'; else; end
if length(Y)==length(t); else; error('Dependent and independent variable vectors must be the same length.'); end
t_nans = t; Y_nans = Y;
t = t(isfinite(Y)); Y = Y(isfinite(Y));
t = t(isfinite(t)); Y = Y(isfinite(t));
if length(BASIS)==(length(PARAMETERS_0)+1); else; error('Incorrect number of parameter initial guesses given.'); end
if nargin > 6
    MaxIt = varargin{1};
else
    MaxIt = 100;
end

J = nan(length(t),length(dBASES));

% Quality of life change: let '1' be shorthand for 'ones(size(t))'
% This shouldn't generally be used (because it's not a nonlinear function
% of a parameter) but it's kept in just in case.
for ii = 1:length(dBASES)
    if strcmp(dBASES{ii},'1')
        dBASES{ii} = 'ones(size(t))';
    else
    end
end

PARAMETERS_i = PARAMETERS_0;
X = nan(size(PARAMETERS_0));
ConvergenceRecord = [];
if nargin > 7 || nargout > 4
    X_record = PARAMETERS_0;
    ResidualVar_Record = [];
else
end
ItNum = 0;

% First condition before OR ensures that we take one step, second
% condition checks tolerance
while (( ~isfinite(sum( abs(X - PARAMETERS_i)./abs(PARAMETERS_i) )) ||  (sum( abs(X - PARAMETERS_i)./abs(PARAMETERS_i) ) > Tol)) && ItNum <= MaxIt )
    
    if sum(isfinite(X)); PARAMETERS_i = X; else; end
    
    for ii=1:length(PARAMETERS_i) % assigns numerical values to arbitrarily named given parameters
        eval([BASIS{ii+1},' = ',num2str(PARAMETERS_i(ii)),';']);
    end
    for ii=1:length(dBASES)
        eval(['J(:,ii) = ',dBASES{ii},';']);
    end
    
    % Directly inverting J'*J shouldn't be a problem because
    % rows(J) >> columns(J)
    eval(['R = Y - (',BASIS{1},');']);
    
    % X = PARAMETERS_i + inv(J'*J)*J'*R; % This is equivalent to the line below in terms of output and time performance
    X = PARAMETERS_i + (J'*J)\J'*R;
    ConvergenceRecord = [ ConvergenceRecord; sum( abs(X - PARAMETERS_i)./abs(PARAMETERS_i) ) ];
    if nargin > 7 || nargout > 4
        X_record = [X_record,X];
        ResidualVar_Record = [ResidualVar_Record,var(R(isfinite(R)))];
    else
    end
    ItNum = ItNum + 1;
end

t = t_nans;
Y = Y_nans;
eval(['Y_fit = ',BASIS{1},';']);

if nargout == 1
    varargout{1} = X;
elseif nargout == 2
    varargout{1} = X;
    varargout{2} = J;
elseif nargout == 3
    varargout{1} = X;
    varargout{2} = J;
    varargout{3} = Y_fit;
elseif nargout == 4
    varargout{1} = X;
    varargout{2} = J;
    varargout{3} = Y_fit;
    varargout{4} = ConvergenceRecord;
elseif nargout == 5
    varargout{1} = X;
    varargout{2} = J;
    varargout{3} = Y_fit;
    varargout{4} = ConvergenceRecord;
    varargout{5} = ResidualVar_Record;
elseif nargout == 6
    varargout{1} = X;
    varargout{2} = J;
    varargout{3} = Y_fit;
    varargout{4} = ConvergenceRecord;
    varargout{5} = ResidualVar_Record;
    varargout{6} = X_record;
end


if nargin > 7
    figure('color',[1 1 1])
    subplot(3,2,[1 3])
    plot(t,Y,'b.-');hold on
    plot(t,Y_fit,'r.-')
    legend({'data','fit'},'FontSize',14)
    xlabel('Independent variable'); ylabel('Dependent variable')
    
    subplot(3,2,[5])
    plot(ResidualVar_Record,'k.-')
    xlabel('Iteration'); ylabel('Residual Variance')
    
    subplot(3,2,[2 4])
    plot([0:(size(X_record,2)-1)]',[X_record./repmat(PARAMETERS_0,1,size(X_record,2))]','.-')
    xlabel('Iteration'); ylabel('X_n./X_0')
    LegendInput = cell(size(X_record,1),1);
    for ii=1:size(X_record,1)
        LegendInput{ii} = ['X(',num2str(ii),')_n / X(',num2str(ii),')_0'];
    end
    legend(LegendInput,'FontSize',14)
    
    subplot(3,2,[6])
    plot(ConvergenceRecord,'k.-');hold on
    plot([0 length(ConvergenceRecord)],Tol*[1 1],':k')
    if max(ConvergenceRecord)>0
        set(gca,'YLim',[0 1.1*max(ConvergenceRecord)])
    else
        set(gca,'YLim',[0 1])
    end
    xlabel('Iteration'); ylabel('[X_n_+_1 - X_n]./[X_n]')
else
end

end

%% Examples

% close all
% TIME = [0:0.1:9.9]'; DATA = 0.2*randn(100,1) + 7*sin(2*pi*0.44*TIME);
% [X,J,Y_fit,ConvergenceRecord,~,X_record] = nonlinear_lsqf(DATA,TIME,{'a*sin(2*pi*f*t)','a','f'},[6,0.40],{'sin(2*pi*f*t)','2*pi*t*a.*cos(2*pi*f*t)'},0.001,   100,[]); disp(X)
% figure; plot(X_record(1,:),X_record(2,:),'.-'); hold on; plot(X_record(1,end),X_record(2,end),'r*')
% 
% %%
% close all
% TIME = [0:0.1:9.9]'; DATA = 0.2*randn(100,1) + exp(-[(TIME - 5)/(3)].^2);
% [X,J,Y_fit,ConvergenceRecord,~,X_record] = nonlinear_lsqf(DATA,TIME,{'exp(-[(t - b)/(a)].^2)','a','b'},[2.5,4.5],...
%     {'(2*(b - t).^2.*exp(-(b - t).^2 /a^2))/a^3','(2*(t - b).*exp(-(b - t).^2 /a^2))/a^2'},...
%     0.001,   100,[]); disp(X)
% figure; plot(X_record(1,:),X_record(2,:),'.-'); hold on; plot(X_record(1,end),X_record(2,end),'r*')
% 
% %%
% close all
% TIME = [0:0.1:9.9]'; DATA = 0.2*randn(100,1) + exp([(TIME - 5)/(3)].^2);
% [X,J,Y_fit,ConvergenceRecord,~,X_record] = nonlinear_lsqf(DATA,TIME,{'exp([(t - b)/(a)].^2)','a','b'},[2.5,4.5],...[3,5]
%     {'-(2*[(b - t).^2].*exp([(b - t)/a].^2))/[a^3]','((2*b - 2*t).*exp([(b - t).^2]/[a^2]))/[a^2]'},...
%     0.001,   100,[]); disp(X)
% figure; plot(X_record(1,:),X_record(2,:),'.-'); hold on; plot(X_record(1,end),X_record(2,end),'r*')
% %%
