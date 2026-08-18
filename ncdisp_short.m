% NCDISP_SHORT(NC_FILE)
% 
% Takes the string of a netCDF file's path as input and outputs a simple
% display of the variables names and dimensions in two columns. No
% overwhelming metadata to sort through, just the essential values you want
% to know most of the time.

function ncdisp_short(NC_FILE)
NCINFO = ncinfo(NC_FILE);
disp([{NCINFO.Variables.Name}',{NCINFO.Variables.Size}'])
end