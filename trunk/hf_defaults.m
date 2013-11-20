% This is file hf_defaults. Put in configurations into this file
%
%
%

hf_path = fileparts(mfilename('fullpath'));
hf_path = strrep(hf_path, '\', '\\');

if isempty(regexp(path, [hf_path pathsep '|' hf_path '$'], 'once'))
    msg = 'hafting-fyhn-tetrode is not yet on your MATLAB path, adding %s';
    warning('hf_defaults:not_in_path', msg, hf_path);
    addpath(hf_path);
end
addpath(fullfile(hf_path, 'hf_tetrode'));
rf_path = fullfile(hf_path, 'rf_pack');
if isdir(rf_path)
    addpath(rf_path);
end