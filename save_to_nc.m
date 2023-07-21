% Save one new variable into a netCDF file in one step. This can overwrite
% any pre-existing variables of the same name given the right option is
% selected.
% 
% IN:
%   FILE = string, full file path of .nc file
%   NAME = string, name of the new variable to be added
%   DATA = data of any type, usually vectors of double or single precision
%   DIM_METADATA = an annoying variable structured as a cell with one or
%          two elements. If one, the element is the name of a known
%          pre-existing dimension in the file. If two, the first is the new
%          dimension's name, the second is its size. NOTE: see the optional
%          variable "MULTIPLE_DIMS" for allowing multiple dimensions for
%          the same variable.
%   METADATA = arbitrarily long cell comprised of length-2 cells of
%          attribute-value pairs, e.g. {'Units','m/s'}
% 
% OPTIONAL INPUTS:
%   OVERWRITE = (default "false") prints the system command needed to
%          delete the variable in question if it already exists. MATLAB has
%          the native "system" function but that does not seem to work with
%          commands that are installed by the user, only native system
%          commands. This variable may be either "true" or "false" (or 1 or 0).
%   SKIP = (default "false") if NAME is already a variable in FILE, then
%          this option, if set to "true", skips the execution of the rest
%          of this script. Leaving this option "false" will lead to an
%          error if the variable already exists. If the variable does not
%          exist, then "save_to_nc" will save it regardless of "SKIP".
%   MULTIPLE_DIMS = (default "false") if true, this allows multiple inputs
%          for the "DIM_METADATA" input, and is necessary for variables
%          with multiple dimensions. In this case, "DIM_METADATA" is
%          structured as a cell of cells, where each cell is either
%          one-element long with a string naming a pre-existing variable,
%          or it is a string-number pair for a new dimension and its size.
% 
% 
% Examples of the more complicated inputs:
% 
% DIM_METADATA = {'LAT',180} or {'LAT'} if you know 'LAT' already exists in
% the file.
% METADATA = {{'Attribute1','Value1'},{'Attribute2','Value2'},...}

function save_to_nc(FILE,NAME,DATA,DIM_METADATA,METADATA,varargin)
% Open the existing netCDF file in write mode

if nargin == 5
    OVERWRITE = false;
    SKIP = false;
    MULTIPLE_DIMS = false;
elseif nargin == 6
    OVERWRITE = varargin{1};
    SKIP = false;
    MULTIPLE_DIMS = false;
elseif nargin == 7
    OVERWRITE = varargin{1};
    SKIP = varargin{2};
    MULTIPLE_DIMS = false;
elseif nargin == 8
    OVERWRITE = varargin{1};
    SKIP = varargin{2};
    MULTIPLE_DIMS = varargin{3};

else
    error(['The function save_to_nc takes 5, 6, or 7 input arguments, while here ' num2str(nargin) ' were given.'])
end

if OVERWRITE
    disp(' ')
    disp(['The following is the bash command for deleting the variable "' NAME '" from the file.'])
    disp(['This requires the nco package to be installed.'])
    disp(' ')
    disp(['ncks -x -v ' NAME ' ' FILE ' ' FILE]);
    disp(' ')
else
end

NCID = netcdf.open(FILE, 'WRITE');

% Get all dimensions already present in the netCDF file:
PRIOR_DIM_IDS = netcdf.inqDimIDs(NCID);
NUM_PRIOR_DIMS = length(PRIOR_DIM_IDS);
PRIOR_DIM_NAMES = cell(NUM_PRIOR_DIMS,1);
PRIOR_DIM_SIZE  = cell(NUM_PRIOR_DIMS,1);
for ii = 1:NUM_PRIOR_DIMS
    [PRIOR_DIM_NAMES{ii},PRIOR_DIM_SIZE{ii}] = netcdf.inqDim(NCID,PRIOR_DIM_IDS(ii));
end

% Get names of existing variables:
EXVARIDS = netcdf.inqVarIDs(NCID); EXVARNAMES = cell(size(EXVARIDS));
for ii = 1:length(EXVARIDS)
    EXVARNAMES{ii} = netcdf.inqVar(NCID,EXVARIDS(ii));
end
% Identify the index of the variable with the same name as NAME
DELVARIND = [];
for ii = 1:length(EXVARIDS)
    if strcmp(NAME,EXVARNAMES{ii})
        DELVARIND = [DELVARIND ; EXVARIDS(ii)];
    else
    end
end

if ~isempty(DELVARIND) && SKIP
warning(['The variable "' NAME '" already exists and was requested to be skipped and not saved over.'])
else



% Define DIMIDS
if MULTIPLE_DIMS
    DIMIDS = nan(length(DIM_METADATA),1);
    for ii = 1:length(DIM_METADATA)
        IS_EXISTING_DIM = zeros(size(PRIOR_DIM_NAMES));
        for jj = 1:length(PRIOR_DIM_NAMES)
            IS_EXISTING_DIM(jj) = strcmp(PRIOR_DIM_NAMES{jj}, DIM_METADATA{ii}{1});
        end
        IS_EXISTING_DIM = sum(IS_EXISTING_DIM);
        if IS_EXISTING_DIM
            DIMIDS(ii) = netcdf.inqDimID(NCID,DIM_METADATA{ii}{1}); % pre-existing dimension
        else
            DIMIDS(ii) = netcdf.defDim(NCID,DIM_METADATA{ii}{1},DIM_METADATA{ii}{2}); % add a new dimension
        end
    end
else
    IS_EXISTING_DIM = zeros(size(PRIOR_DIM_NAMES));
    for jj = 1:length(PRIOR_DIM_NAMES)
        IS_EXISTING_DIM(jj) = strcmp(PRIOR_DIM_NAMES{jj}, DIM_METADATA{1});
    end
    IS_EXISTING_DIM = sum(IS_EXISTING_DIM);
    if IS_EXISTING_DIM
        DIMIDS = netcdf.inqDimID(NCID,DIM_METADATA{1}); % pre-existing dimension
    else
        DIMIDS = netcdf.defDim(NCID,DIM_METADATA{1},DIM_METADATA{2}); % add a new dimension
    end
end

% Create the new variable
VARID = netcdf.defVar(NCID, NAME, class(DATA), DIMIDS);

% Add metadata attributes to the new variable
for ii = 1:length(METADATA)
    netcdf.putAtt(NCID, VARID, METADATA{ii}{1}, METADATA{ii}{2});
end

% Assign values to the new variable if desired
netcdf.putVar(NCID, VARID, DATA);

% Close the netCDF file
netcdf.close(NCID);

end
end

% save_to_nc('/full/path/test.nc',...
%            'new_variable',[1 2 3 4 5]',{'new_dim',5},{{'units','donuts'},{'importance','minimal'}})
% 
% Adds this new variable:
% 
% new_variable
%        Size:       5x1
%        Dimensions: new_dim
%        Datatype:   double
%        Attributes:
%                    units      = 'donuts'
%                    importance = 'minimal'
