function write_csv(in_struct, ofile, varargin)
%WRITE_CSV Writes a cvs file that libreoffice/excel can open
%
% WRITE_CSV(S, ofile) writes OSI and total spike count for each
% orientation and each cell to file *ofile*. If *ofile* has no
% extension, '.csv' will be added as extension.
%
% The required arguments can be followed by parameter/value pairs:
%     'expanded' : boolean or int (default true)
%         If true (default), fourier OSI will also be written to file
%     'separator' : str (default ';')
%         Each entry in *ofile* is separated by this string.
%         libreoffice/excel uses this to separate data into different
%         cells.
%
% Version 0.1 September 2012, Eivind Skjønsberg Norheim
% Version 0.2 January 2013, ESN:
%     Changes:
%         - Changed name from my_csvwrite_expanded() to write_csv()
%         - New parameter/value pair 'expanded'
%         - Function checks if path to *ofile* exists, and if not creates
%           the whole path
%

% Parsing input
p = inputParser();
def_expanded = true;
def_separator = ';';
addRequired(p, 'in_struct', @isstruct)
addRequired(p, 'ofile', @isstr)
addParamValue(p, 'expanded', def_expanded, @is_my_logical)
addParamValue(p, 'separator', def_separator, @isstr)

parse(p, in_struct, ofile, varargin{:})

sep = p.Results.separator;
expanded = p.Results.expanded;

% if no extension add csv as extension
[pathstr, ~, ext] = fileparts(ofile);
if isempty(ext)
    ofile=[ofile, '.csv'];
end

% Check if file *ofile* exists,
reply = '';
while exist(ofile, 'file') && ~strcmpi(reply, 'w')  
    pstr = sprintf(['\nFile %s exists.\nOver[w]rite, [r]ename or [q]uit?',...
        ' w/r/q: '], ofile);    
%    while isempty(reply)
    reply = input(pstr, 's');
    if strcmpi(reply, 'w')
        fprintf(1, 'Ok, I\''ll overwrite.\n');
    elseif strcmpi(reply, 'r')
        rename_str = sprintf('New name: ');
        fname_new = input(rename_str, 's');
        ofile = fullfile(pathstr, [fname_new, '.csv']);
    elseif strcmpi(reply, 'q')
        fprintf('Ok, aborting.\n');
        return
    else
        fprintf('Not a valid option, try again\n\n');
        reply = '';
    end
%    end
end

orientations = in_struct.orientations;
norient = length(orientations);

% Create directory if it doesn't exist and open file for writing
if ~exist(pathstr, 'dir')
    fprintf(1, 'Directory %s does not exist. Creating directory\n', pathstr);
    mkdir(pathstr);
end
fid = fopen(ofile, 'w');

% Write header
sline = ['Cell #', sep, 'OSI'];
if expanded
    sline = [sline, sep, 'fOSI'];
end   

%    sprintf(';%g', orientations), '\n'];
for iorient=1:norient
    orientation = orientations(iorient);    
    sline = [sline, sep, num2str(orientation)]; %#ok
end
sline = [sline ,'\n'];

slines = {['Session', sep, in_struct.sessionfile, '\n'], ...
    ['Cut', sep, in_struct.cutfile, '\n'], ...
    ['Stim synch', sep, in_struct.stimfile, '\n\n'], ...
    sline};

for iline=1:length(slines)
    sline = slines{iline};
    % If the file didn't open successfully, this is where the error occurs
    try
        fprintf(fid, sline, 'char');
    catch ME
        error('el_phys:write_csv:write_header_error', ME.message)
        return
    end
end

for icell=1:length(in_struct.cells)
    this = in_struct.cells(icell);
    sline = [num2str(icell), sep, num2str(this.OSI)];
    if expanded
        sline = [sline, sep, num2str(this.fOSI)]; %#ok
    end
    for iorient=1:norient
        nspikes = this.nspikes(iorient);
        sline = [sline, sep, num2str(nspikes)]; %#ok
    end
    % Changing decimal delimiter from '.' to ',' unless *sep* has been set
    % to ','
    if ~strcmp(sep, ',')
        sline(sline=='.')=',';
    end
    try
        fprintf(fid, [sline, '\n'], 'char');
    catch ME
        error('el_phys:write_csv:write_bulk_error', ME.message)
        return
    end
end

fclose(fid);
