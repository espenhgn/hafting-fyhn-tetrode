function out = input_wrapper(input_file)
%INPUT_WRAPPER    Wraps around get_tetrode_data

fid = fopen(input_file, 'r');
if fid == -1
    msg = ['Could not open the input file! Make sure the filname and',...
        'path are correct.'];
    error('hf_tetrode:input_wrapper:fopen_error', msg);
end

while 1
end