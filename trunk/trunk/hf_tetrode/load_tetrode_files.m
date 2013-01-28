function out = load_tetrode_files(sessionfiles, varargin)
%LOAD_TETRODE_FILES    load files used for tetrode analysis
% S = LOAD_TETRODE_FILES(sessionfiles) where *sessionfiles* is cell array
% containing (relative or absolute)  paths to files. Files in *sessionfile*
% are tetrode data files created by AXONA.
% 
% Required arguments can be followed by parameter/value pairs:
%     'cutfiles' : cell (default {})
%         array containing paths to cut files. If naming convention from
%         Tint is followed, and there's one cut file for each session  file
%         in *sessionfiles*, this can be omitted. Specify this if (i) other
%         than standard naming convention is used, or if (ii) joint cutfile
%         is available.
%     'stimfiles' : cell (default {})
%         Absolute or relative path to stimulus log created by older
%         versions of PsychStimController. The stimulus log contains at
%         least 'logorientation', and for newer verisions also 't0'.
%     't_blank' : double array (default zeros(size(sessionfiles)))
%         Length of grey period between stimuli in seconds. Overridden if
%         stimfiles contain values for start of gray periods or if there's
%         a .inp file present
%     't_delay' : double array (default zeros(size(sessionfiles)))
%         Time delay between button press and stimulus onset in seconds.
%         Overridden if stimfiles contain values for start of stimuli (and
%         not the end of stimuli as it does in early versions). Also
%         overridden by .inp files.
%     'verbose' : boolean or int (default false)
%         Toggles messages on/off. false means no messages, true means 
%         print messages to standard output. Integer values are also
%         accepted: 0 and negative values means false, positive values
%         means true.

% This file takes all file paths, checks for existence (catch exceptions),
% loads files, checks that the sessionfiles and cut match, 
% return timestamps, clusters, ncells, logtimes,  and orientations. 
% ***Return it all in either struct or cell format***

% Default cut scenario: Each sessionfile has it's own cutfile and own .inp
% file. -> check that *cut_for* is correct separately
% 1st cut scenario: Each sessionfile has a separate cutfile, with
% non-standard name -> same chack as above
% 2nd cut scenario: All sessionfiles has the same cutfile -> Only include
% the sessionfiles that are contained in cutfile
% NB! Cutfile only contains parts of sessionfile's name. Check that number
% of spikes is the same as number of cluster indices

% stimfile scenarios: sync is either specified in stimfile or in default
% .inp file with same name as sessionfile.

% Parse input
p = inputParser();
def_cutfiles = {};
def_stimfile = {};
def_tblank = zeros(size(sessionfiles));
def_tdelay = zeros(size(sessionfiles));
def_verbose = false;

addRequired(p, 'sessionfiles', @iscell)
addParamValue(p, 'cutfiles', def_cutfiles, @iscell)
addParamValue(p, 'stimfiles', def_stimfile, @iscell)
addParamValue(p, 'verbose', def_verbose, @is_my_logical)
addParamValue(p, 't_blank', def_tblank, @isnumeric)
addParamValue(p, 't_delay', def_tdelay, @isnumeric)
parse(p, sessionfiles, varargin{:})

stimfiles = p.Results.stimfiles;
cutfiles = p.Results.cutfiles;
t_blank = p.Results.t_blank;
t_delay = p.Results.t_delay;
verbose = p.Results.verbose;

% Initialize some stuff
nsess = length(sessionfiles);
nspikes = 0;
timestamps = cell(nsess, 1);
fnames = cell(nsess, 1);
sync_data = struct();

% Try to load spikes from all sessionfiles and put into cell array,
% and load syncronization data from stimfiles or default .inp file
if ~isempty(stimfiles) && length(stimfiles)~=nsess        
    msg = 'Same number of sessionfiles and stimfiles needed';        
    error('el_phys:load_tetrode_files:input_error', msg)
end

for isess=1:nsess
    %%%%%%%%%    Get spikes    %%%%%%%
    sessionfile = sessionfiles{isess};
    [pathstr, fname, ~] = fileparts(sessionfile);
    try
         tstamps = getspikes(sessionfile);         
    catch ME
        error('el_phys:load_tetrode_files:getspikes_error', ME.message)
        out = -1;
        return
    end
    fnames{isess} = fname;
    nspikes = nspikes + length(tstamps);
    timestamps{isess} = tstamps;
    %%%%%%%%     Read stimfile/.inp file     %%%%%%%%%
    % stimfiles is either empty (use default file) or has same length as
    % sessionfiles
    if ~isempty(stimfiles)
        stimfile = stimfiles{isess};
    else
        % Default .inp file has same name as sessionfile and is located in
        % same folder
        stimfile = fullfile(pathstr, [fname, '.inp']);
        stimfiles{isess} = stimfile;
    end
    try
        [sync_data(isess).logtime, sync_data(isess).logtime_blank, ...
            sync_data(isess).orients] = get_sync(stimfile,...
            't_blank', t_blank(isess),'t_delay', t_delay(isess),...
            'verbose', verbose);
    catch ME
        error('el_phys:load_tetrode_files:get_sync_error', ME.message)
        out = -1;
        return
    end
end

%%%%%    Load cut data    %%%%
clusts = cell(nsess,1);
cut_tmp = {};

if isempty(cutfiles) || length(cutfiles)==nsess
    for isess=1:nsess
        sessionfile = sessionfiles{isess};
        [pathstr, fname, ext] = fileparts(sessionfile);
        % Use default name if not stated explicitly
        if isempty(cutfiles)
            cutfile = fullfile(pathstr, [fname, '_', ext(2:end), '.cut']);
            cut_tmp{isess} = cutfile;
        else
            cutfile = cutfiles{isess};
        end
        try
            [clust, cut_for] = getcut(cutfile);
        catch ME
            error('el_phys:load_tetrode_files:single_getcut_error', ...
                ME.message)
            out = -1;
            return
        end
        % Check if cut is for this and only this session file. We can't
        % know if the cut is for the right tetrode (except if standard
        % naming has been used)
        if length(cut_for)~=1
            msg = ['Cut file %s is cut for several files. I don''t', ...
                ' know how to extract single cut from this yet'];
            error('el_phys:load_tetrode_files:cut_for_several', ...
                msg, cutfile)
        elseif ~strcmp(fname, cut_for{1})
            msg = 'Cut file %s not cut for session file %s';
            error('el_phys:load_tetrode_files:not_cut_for', msg, ...
                cutfile, sessionfile)
        elseif length(clust)~=length(timestamps{isess})
            msg = ['Number of spikes from cutfile %s and session file', ...
                ' %s do not match. (%d vs %d). Are you sure this is ',...
                'right cut for this session file?'];
            error('el_phys:load_tetrode_files:spike_count_error', msg, ...
                cutfile, sessionfile, length(clust),...
                length(timestamps{isess}))
        else
            % If none of exceptions are thrown then keep the cut data
            clusts{isess} = clust;
        end
    end
    if isempty(cutfiles)
        cutfiles = cut_tmp;
    end
elseif length(cutfiles)==1 && nsess>1
    % This is the case with joint cut file.
    cutfile = cutfiles{1};
    try
        [joint_clust, cut_for] = getcut(cutfile);
    catch ME
        error('el_phys:load_tetrode_files:joint_getcut_error', ME.message)
        out = -1;
        return
    end
    ncut = length(cut_for);
    % Check that number of cuts in cut file is same as number of session
    % files, and that total number of cut entries is same as total
    % number of spikes.
    if ncut~=nsess
        msg = ['Number of session files is %d, but cut file %s is cut ',...
            'for %d files'];
        error('el_phys:load_tetrode_files:cut_count_error', msg, nsess,...
            cutfile, ncut)    
    elseif length(joint_clust)~=nspikes
        msg = ['Number of spikes from cutfile %s and provided session',...
            ' files do not match. (%d vs %d). Are you sure this is',...
            ' the right cut for these session files?'];
        error('el_phys:load_tetrode_files:joint_spikecnt_error', msg,...
            cutfile, length(joint_clust), nspikes)
    end
    for icut=1:ncut
        % Loop through cuts, find the right session to assign the cut to
        % and number of spikes for this session
        isess = find(strcmp(cut_for{icut}, fnames));
        if isempty(isess)
            msg = ['Cut file has cut for %s, but this is not one of the',...
                ' session files provided as input'];
            error('el_phys:load_tetrode_files:cut_missing_sessfile', ...
                msg, cut_for{icut})
        elseif length(isess)>1
            msg = 'Duplicate session file';
            error('el_phys:load_tetrode_files:duplicate_sessfile', msg)
        end
        nstamps = length(timestamps{isess});
        % Pop the cluster values for this session and put it in clusts
        clusts{isess} = joint_clust(1:nstamps);
        joint_clust(1:nstamps) = [];
    end
end

out = struct('timestamps', {timestamps}, 'sync_data', sync_data,...
    'clusts', {clusts}, 'cutfiles', {cutfiles}, 'stimfiles', {stimfiles});
