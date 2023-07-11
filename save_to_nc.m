% Save one new variable into a netCDF file in one step
% 
% IN:
%   FILE = string, full file path of .nc file
%   NAME = string, name of the new variable to be added
%   DATA = data of any type, usually vectors of double or single precision
%   DIM_METADATA = an annoying variable structured as a cell with one or
%          two elements. If one, the element is the name of a known
%          pre-existing dimension in the file. If two, the first is the new
%          dimension's name, the second is its size (scalar, it seems we
%          can only write vectors here)
%   METADATA = arbitrarily long cell comprised of length-2 cells of
%          attribute-value pairs, e.g. {'Units','m/s'}
% 
% 
% 
% Examples of the more complicated inputs:
% 
% DIM_METADATA = {'LAT',180} or {'LAT'} if you know 'LAT' already exists in
% the file.
% METADATA = {{'Attribute1','Value1'},{'Attribute2','Value2'},...}

function save_to_nc(FILE,NAME,DATA,DIM_METADATA,METADATA)
% Open the existing netCDF file in write mode
NCID = netcdf.open(FILE, 'WRITE');

% Get all dimensions already present in the netCDF file:
PRIOR_DIM_IDS = netcdf.inqDimIDs(NCID);
NUM_PRIOR_DIMS = length(PRIOR_DIM_IDS);
PRIOR_DIM_NAMES = cell(NUM_PRIOR_DIMS,1);
PRIOR_DIM_SIZE  = cell(NUM_PRIOR_DIMS,1);
for ii = 1:NUM_PRIOR_DIMS
    [PRIOR_DIM_NAMES{ii},PRIOR_DIM_SIZE{ii}] = netcdf.inqDim(NCID,PRIOR_DIM_IDS(ii));
end

% Define DIMIDS
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

% Create the new variable
VARID = netcdf.defVar(NCID, NAME, class(DATA), DIMIDS);

% Add metadata attributes to the new variable
for ii = 1:length(METADATA)
    netcdf.putAtt(NCID, VARID, METADATA{ii}{1}, METADATA{ii}{2});
end

% Assign values to the new variable if desired
whos NCID VARID DATA %$
netcdf.putVar(NCID, VARID, DATA);

% Close the netCDF file
netcdf.close(NCID);

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
