% ncreadall(NC_file) loads all of the variables in a .nc and names them in
% MATLAB the same as they are named in the file. WARNING: There is no way
% to tell if the .nc file is very large (i.e. if its contents would fill of
% the memory of the MATLAB workspace). Therefore, be aware of the size of
% the file before using this function.
% 
% IN:   NC_file  = string of the .nc file name, e.g. 'data.nc'
%       EXCLUDED = (optional) cell of strings of names of variables to be
%                  excluded
% 
% OUT:  var_struct = OPTIONAL, this is a structure with the variables
%                    stored in it. If you choose not to have an output
%                    for ncreadall, then the variables in NC_file simply
%                    go to your workspace separately, rather than being
%                    bundled into a structure.

% Luke Kachelein - February 21, 2018
% 
% For use to meet course requirements and perform academic research at
% degree-granting institutions only. It is not available for government,
% commercial, or other organizational use, unless granted permission by the
% original author.

function var_struct = ncreadall(NC_file,varargin)

if nargin == 1
    EXCLUDED = {''};
elseif nargin == 2
    EXCLUDED = varargin{1};
else
    error('One or two inputs expected.')
end

if sum(NC_file((end-2):end) == '.nc') == 3
else
    help ncreadall
    error('NCREADALL only accepts .nc files. You must terminate the string argument "NC_file" with ".nc". Read the instructions.')
end

INFO = ncinfo(NC_file);

variables_struct = INFO.Variables; % this will be a struct
N = length(variables_struct); % number of variables

if length(EXCLUDED) == 1 && strcmp(EXCLUDED{1},'')
    variables_cell = cell(1,N);
else
    variables_cell = {};
end

n_var = 1;
for n = 1:N
    % Exclude variables if request:
    if ~sum(strcmp(variables_struct(n).Name,EXCLUDED)) % 0 if the n'th variable is in EXCLUDED
        variables_cell{n_var} = variables_struct(n).Name;
        n_var = n_var + 1;
    else
    end
end

% Redefine N after eliminating possible undesired variables:
N = length(variables_cell);

% If no output is expected, just define variables in workspace
if nargout == 0
    for n = 1:N
        eval([variables_cell{n},' = ncread(''',NC_file,''',''',variables_cell{n},''');']);
        eval(['assignin(''base'',''',variables_cell{n},''',',variables_cell{n},');']);
    end

% If one is expected, define a structure
elseif nargout == 1
    var_struct = struct;
    for n = 1:N
        eval(['var_struct.',variables_cell{n},' = ncread(''',NC_file,''',''',variables_cell{n},''');']);
    end

else
    help ncreadall
    error('NCREADALL allows only one or zero output arguments. Read the instructions.')
end

end
