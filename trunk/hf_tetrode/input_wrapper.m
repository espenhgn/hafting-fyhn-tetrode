function out = input_wrapper(input_file)
%INPUT_WRAPPER    Wraps around *get_tetrode_data()*
% S = INPUT_WRAPPER(input_file) calls *get_tetrode_data()* with arguments
% stated in input_file. The input file is a plain text document in this
% format:
%     session <path_to_session_file_#1>
%     session <path_to_session_file_#2>
%     synch <path_to_synch_file_#1>
%     synch <path_to_synch_file_#2>
%     cut <path_to_cut_file_#1>
%     cut <path_to_cut_file_#1>
% 
% If ALL synch and/or ALL cut files follow standard naming conventions,
% these files can be omitted. 
% If all session files have been cut together in one cut file, this cut
% file only needs to be specified once.
% In addition other optional arguments to *get_tetrode_data()* can be
% specified. Write the name and value of the argument separated by a blank
% space.

sessionfiles = {};
stimfiles = {};
cutfiles = {};

ncomment = 0;

fid = fopen(input_file, 'r');
if fid == -1
    msg = ['Could not open the input file! Make sure the filname and',...
        'path are correct.'];
    error('hf_tetrode:input_wrapper:fopen_error', msg);
end

while ~feof(fid)
    sline = fgetl(fid);
    contents = regexp(sline, ' ', 'split');
    if regexp(contents{1}, '%', 'start')==1
        ncomment=ncomment+1;
    elseif length(contents)~=2
        nspace = length(contents)-1;
        msg = ['Format of this line in %s is wrong:\n %s\n',... 
            'I allow only one space in each line, but I see %d spaces\n'];
        error('hf_tetrode:input_wrapper:space_error', msg, input_file,...
            sline, nspace)
    else
        [var, val] = contents{:};
        switch lower(var)
            case 'session'    
                sessionfiles{end+1} = val; %#ok
            case 'synch'
                stimfiles{end+1} = val; %#ok
            case 'cut'
                cutfiles{end+1} = val; %#ok
            case 't_delay'
                t_delay = str2num(val); %#ok                
            case 't_blank'
                t_blank = str2num(val); %#ok
            case 'export_csv'
                export_csv = val;
            case 'export_args'
                eval([var, '=', val, ';']);
            case 'verbose'
                eval([var, '=', val, ';']);
            otherwise
                msg = ['Input variable %s in input file %s is unknown.',...
                    ' Check spelling!'];
                error('hf_tetrode:input_wrapper:variable_error', msg, var,...
                    input_file)
        end
    end
end

fclose(fid);
args = {'stimfiles', stimfiles, 'cutfiles', cutfiles, ...
    't_delay', t_delay, 't_blank', t_blank, ...
    'export_csv', export_csv, 'export_args', export_args, ...
    'verbose', verbose};
out = get_tetrode_data(sessionfiles, args{:});