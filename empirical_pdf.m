% Y_pop = empirical_pdf(PDF_x, y, B_N, lims);
% Y_pop = empirical_pdf(PDF_x, y, B_N);
% 
% Empirical PDF of y = y(x), where x is a value with a known or assumed
% PDF.
% 
%                                   IN:
% 
% PDF_x = String or cell of strings.
%         If string: functional form of PDF(x), e.g. 'exp(-x.^2)' or
%         If cell:   first element is like above but with constants, e.g. 'exp(-(x - m).^2)'
%                    second etc. element(s) is/are the constants, e.g. 'm = 2'
%         Alternative cell form:
%                    first element is 'gaussian' or 'g', which stands in for 'exp(-(x - m).^2 / (2*s^2))'
%                    second element is 'm = M' where M is explicitly a number, the mean.
%                    third element is 's = S' where S explicitly a number, the standard deviation = sqrt(variance).
%         NOTE: x is to be treated as a vector in this argument, i.e.
%         'x.^2' is correct syntax but 'x^2' is not.
% y     = y(x) in funtional form, e.g. 'x.^2'. Adhere to element-wise
%         operations on x.
% B_N   = Two-element vector:
%         B_N(1) = B, number of bins between lims(1) and lims(2) (see
%                  below) in which discretization will be set.
%         B_N(2) = N, average number of realizations fo x that will be in
%                  each bin.
%         This means that x will be a vector of B*N elements with a
%         histogram that approximates the analytical form of PDF_x.
% lims  = Upper and lower limits of the interval on which you want this
%         function to look. OPTIONAL if "PDF_x" is a cell with the first
%         element either 'gaussian' or 'g', in which case the interval is
%         [m - 3s, m + 3s]
% 
% 
%                                   OUT:
% 
% Y_pop = Vector of length (approximately) B*N. The result of running a
%         population of x's that follow PDF_x through y(x) to get a
%         population of y's. From this can be obtained the desired
%         statistical parameters.
% 

function varargout = empirical_pdf(PDF_x, y, B_N, varargin)

if iscell(PDF_x) && length(PDF_x) > 1
    pdf_x = PDF_x{1};
    for ii = 2:(length(PDF_x))
        eval([PDF_x{ii},';'])
    end
    if strcmp(pdf_x,'gaussian') || strcmp(pdf_x,'g')
        pdf_x = 'exp(-(x - m).^2 / (2*s^2))';
    else
    end
elseif ischar(PDF_x)
    pdf_x = PDF_x;
elseif iscell(PDF_x) && length(PDF_x) == 1
    pdf_x = PDF_x{1};
else
    error('Incorrectly formatted "PDF_x" input.')
end

if nargin == 4
    lims = varargin{1};
elseif nargin == 3 && strcmp(PDF_x{1},'gaussian')
    lims = [m - 3*s, m + 3*s];
elseif nargin == 3 && strcmp(PDF_x{1},'g')
    lims = [m - 3*s, m + 3*s];
else
    error('Input error: either wrong number of inputs or a formatting error.')
end

B = B_N(1);
N = B_N(2);

dx = (lims(2) - lims(1))/B;
x = (lims(1)):(dx):(lims(2));

eval(['PDF_x_discrete = ',pdf_x,';'])

% Normalize:
PDF_x_discrete = round(B*N*PDF_x_discrete/(sum(PDF_x_discrete)));

% Build population of x's (X) that will match the PDF:
X = [];
X = [X, (x(1)):(dx/PDF_x_discrete(1)):(x(1) + dx/2)];
for jj = 2:(length(PDF_x_discrete) - 1)
    X = [X, (x(jj)):(dx/PDF_x_discrete(jj)):(x(jj) + dx)];
end
X = [X, (x(end) - dx/2):(dx/PDF_x_discrete(end)):(x(end))];

% This part is confusing, because we are overwriting "x" with the data in
% "X" so that the input "y" makes sense (I didn't want to force the user to
% distinguish between the variable x and population X):
foo = x;
x = X;

eval(['Y_pop = ',y,';'])
% figure;histogram(Y_pop,'Normalization','pdf')

% Possibly unnecessary, otherwise switch "varargout" at beginning of function with "PDF_y":
if nargout == 1
    varargout{1} = Y_pop;
elseif nargout == 0
    disp(Y_pop)
else
    error('Wrong number of outputs (should be zero or one).')
end

end
