% Y_pop = empirical_pdf(PDF_x, B_N, y, lims);
% Y_pop = empirical_pdf(PDF_x, B_N, y,);
% 
% Empirical PDF of y = y(x), where x is a value with a known or assumed
% PDF.
% 
%                                   IN:
% 
% PDF_x = A 'string' or {cell of 'strings'} or {cell of {cells of 'strings'}}.
%         If string: functional form of PDF(x), e.g. 'exp(-x.^2)' or
%         If cell of strings:
%                    -first element is like above but with constants, e.g. 'exp(-(x - m).^2)'
%                    -second etc. element(s) is/are the constants, e.g. 'm = 2'
%         Alternative cell form:
%                    -first element is 'gaussian', 'g', 'Gaussian', or 'G',
%                     which stands in for 'exp(-(x - m).^2 / (2*s^2))'
%                    -second element is 'm = M', where M is explicitly a
%                     number, the mean.
%                    -third element is 's = S' where S explicitly a number,
%                     the standard deviation = sqrt(variance).
%         If cell of cells:
%                    Each cell is formatted like the cell-of-strings input
%                    format, but for a different variable, e.g.:
%                    { {'exp(-(x1 - m1).^2)', 'm1 = 2'}, ...
%                      {'exp(-(x2 - m2).^2)', 'm2 = -1'} }
%                    Unsurprisingly, the input "y" should have x1 and x2 in
%                    it.
%         NOTE: x is to be treated as a vector in this argument, i.e.
%         'x.^2' is correct syntax but 'x^2' is not.
% 
% B_N   = A two-element vector:
%         B_N(1) = B, number of bins between lims(1) and lims(2) (see
%                  below) in which discretization will be set.
%         B_N(2) = N, average number of realizations of x that will be in
%                  each bin.
%         This means that x will be a vector of B*N elements with a
%         histogram that approximates the analytical form of PDF_x.
%         If PDF_x is a cell of cells, then B_N applies to each of x1, x2,
%         etc.
% 
% y     = y(x) in funtional form, e.g. 'x.^2'. Adhere to element-wise
%         operations on x (or x1, x2, etc. for multivariable functions).
% 
% lims  = Upper and lower limits of the interval on which you want this
%         function to look. OPTIONAL if "PDF_x" is a cell with the first
%         element either 'gaussian' or 'g', in which case the interval is
%         [m - 3s, m + 3s].
%         If "PDF_x" is a cell of cells, even if they're Gaussian
%         functions, the intervals must be specified like above as elements
%         of yet another cell (corresponding to each variable described in
%         the sub-cells of "PDF_x"), because accounting for all the x_n's
%         and their intervals of interest is too difficult in the case that
%         at least one of x_n does not follow a Gaussian distribution.
% 
% 
%                                   OUT:
% 
% Y_pop = Vector of length B*N (approximately, due to rounding). The result
%         of running a population of x's that follow PDF_x through y(x) to
%         get a population of y's. From this can be obtained the desired
%         statistical parameters.
%         This variable is deterministic for y(x), but will necessarily
%         have some randomness for y(x1,x2,...).
% 

function varargout = empirical_pdf(PDF_x, B_N, y, varargin)

if iscell(PDF_x)
else
    PDF_x = {PDF_x};
end

if iscell(PDF_x{1})
    % Warning: this is going to break a lot of rules of good programming
    % practices for the sake of convenience (heavy use of the "eval"
    % function). If this breaks your computer, you did something wrong
    % because it definitely shouldn't be that bad.
    if iscell(PDF_x{1})
    else
        error('Incorrectly formatted input "PDF_x".')
    end
    
    if nargin ~= 4
        error('Input variable "lims" is required for multivariable PDF''s.')
    else
        lims_cell = varargin{1};
    end
    
    B = B_N(1);
    N = B_N(2);
    
    %%
    % Loop through all the variables (1:length(PDF_x)), heavily relying on
    % the "eval" function:
    
    for n = 1:length(PDF_x)
        
        pdf_xn = PDF_x{n}{1};
        for ii = 2:(length(PDF_x{n}))
            eval([PDF_x{n}{ii},';'])
        end
        if strcmp(pdf_xn,'gaussian') || strcmp(pdf_xn,'g') || strcmp(pdf_xn,'Gaussian') || strcmp(pdf_xn,'G')
            pdf_xn = ['exp(-(x',num2str(n),' - m',num2str(n),').^2 / (2*s',num2str(n),'^2))'];
        else
        end
        
        dxn = (lims_cell{n}(2) - lims_cell{n}(1))/B;
        xn = (lims_cell{n}(1)):(dxn):(lims_cell{n}(2));
        eval(['x',num2str(n),' = xn;'])
        
        PDF_x_discreten = eval([pdf_xn,';']);
        
        % Normalize:
        PDF_x_discreten = round(B*N*PDF_x_discreten/(nansum(PDF_x_discreten)));
        
        % Build population of x's (X) that will match the PDF:
        Xn = [];
        Xn = [Xn, (xn(1)):(dxn/PDF_x_discreten(1)):(xn(1) + dxn/2)];
        for jj = 2:(length(PDF_x_discreten) - 1)
            Xn = [Xn, (xn(jj)):(dxn/PDF_x_discreten(jj)):(xn(jj) + dxn)];
        end
        Xn = [Xn, (xn(end) - dxn/2):(dxn/PDF_x_discreten(end)):(xn(end))];
        
        % This part is confusing, because we are overwriting "x" with the data in
        % "X" so that the input "y" makes sense (I didn't want to force the user to
        % distinguish between the variable x and population X):
        eval(['x',num2str(n),' = Xn;'])
        
    end
    
    % We should now have x1, x2, x3, etc., which will have almost the same
    % length, but there were rounding errors. Therefore, we need to
    % randomly remove elements from all bu the shortest x_n until they all
    % have the same length as the shortest one (necessary in calculating
    % y).
    
    x_Population_sizes = [];
    for n = 1:length(PDF_x)
        x_Population_sizes(n) = eval(['length(x',num2str(n),');']);
    end
    Min_Pop = min(x_Population_sizes);
    
    % Randomly rearrange all the vectors and cut them off at the "Min_Pop"
    % index (that way the ones you remove are random and the population of
    % x1, x2, etc. are out of order and hence uncorrelated for calculating
    % the population of y; this is a double whammy):
    for n = 1:length(PDF_x)
        eval(['x',num2str(n),' = x',num2str(n),'(randperm(length(x',num2str(n),')));']);
        eval(['x',num2str(n),' = x',num2str(n),'(1:Min_Pop);']);
    end
    
    % y is defined after the if-statement
    %%
else
    if length(PDF_x) > 1
        pdf_x = PDF_x{1};
        for ii = 2:(length(PDF_x))
            eval([PDF_x{ii},';'])
        end
        if strcmp(pdf_x,'gaussian') || strcmp(pdf_x,'g') || strcmp(pdf_x,'Gaussian') || strcmp(pdf_x,'G')
            pdf_x = 'exp(-(x - m).^2 / (2*s^2))';
        else
        end
    elseif ischar(PDF_x{1})
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
    
    PDF_x_discrete = eval(pdf_x);
    
    % Normalize:
    PDF_x_discrete = round(B*N*PDF_x_discrete/(nansum(PDF_x_discrete)));
    
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
    x = X;
    % figure;histogram(x,'Normalization','pdf');title('histogram of x')
end

%%

eval(['Y_pop = ',y,';'])
% figure;histogram(Y_pop,'Normalization','pdf');title('histogram of y')

% Possibly unnecessary, otherwise switch "varargout" at beginning of function with "PDF_y":
if nargout == 1
    varargout{1} = Y_pop;
elseif nargout == 0
    disp(Y_pop)
else
    error('Wrong number of outputs (should be zero or one).')
end

end

% % EXAMPLE:
% Y_pop = empirical_pdf({{'g','m1 = 10','s1 = 5'},{'g','m2 = 10','s2 = 5'}}, [1000,20], 'atan2(x1,x2)', {10+[-15 15],10+[-15 15]});
