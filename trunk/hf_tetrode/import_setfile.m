function [setfile] = import_setfile(filename)
%IMPORT_SETFILE Import all information from axona .set file as a struct.
%
%   The returned argument is a struct, where
%   each field name is set similar to variable name in the .set-file, 
%   and for now, only numeric values are returned as either int or double.
%   The remaining will be strings.
%
% Example, getting field name trial_time from a .set-file:
%   [setfile] = import_setfile('path/to/file.set');
%   setfile.trial_time

fprintf(1, ['reading setfile ', filename, '\n']);
fileID = fopen(filename, 'r');

if fileID == -1
    msg = ['Could not open the set file! Make sure the filname and',...
        'path are correct.'];
    error('hf_tetrode:import_setfile:fopen_error', msg);
end

dataArray = textscan(fileID, '%s', 'Delimiter', '');

setfile = {};

setfile.filename = filename;

for ii = 1:length(dataArray{1})
    x = strsplit(dataArray{1}{ii});
    field = x{1};
    value = strjoin(x(2:end));
    if sum(isletter(value)) > 0 || isnan(str2double(value))
        setfile.(field) = value;
    else
        setfile.(field) = str2double(value);
    end
end
