function [stim_t, blank_t, orients] = get_sync(log_file, varargin)

%GET_SYNC    Get synchronization data and orientations from log file
% [stim_t, blank_t, orients] = GET_SYNC(log_file) checks if log_file is
% Axona .inp file or .mat file and collects synchronization information
% from file. 
% The file name can be followed by parameter/value pairs:
%     't_blank' : float (default 0)
%         Length of grey period between stimuli in seconds. Overridden if
%         stimfile contains values for start of gray periods
%     't_delay' : float (default 0)
%         Time delay between button press and stimulus onset in seconds.
%         Overridden if stimfile contains values for start of stimuli (and
%         not the end of stimuli as it does in early versions)
%     'verbose' : boolean or int (default false)
%         Toggles messages on/off. false means no messages, true means 
%         print messages to standard output. Integer values are also
%         accepted: 0 and negative values means false, positive values
%         means true.
%
% Born in version 0.2, Jan 2013, Eivind Skjønsberg Norheim

% arg handling
p = inputParser();
def_verbose = false;
def_tblank = 0;
def_tdelay = 0;
addRequired(p, 'log_file', @isstr)
addParamValue(p, 't_blank', def_tblank, @isnumeric)
addParamValue(p, 't_delay', def_tdelay, @isnumeric)
addParamValue(p, 'verbose', def_verbose, @is_my_logical)

parse(p, log_file, varargin{:})

verbose = p.Results.verbose;
% Done handling arguments

stim_t = 0;
blank_t = 0;
orients = 0;

[~, ~, ext] = fileparts(log_file);

switch ext
    case '.inp'
        % Read the .inp file and pick the values where the type 
        % was 'I' for input
        [ts_tmp, type, keyval_tmp] = getinp(log_file);
        ts = ts_tmp(type=='I');
        keyval = keyval_tmp(type=='I');
        
        % Since we strobe the values we only look for double entries (some
        % entries might even be triple or more)
        idxs = double_entries(keyval);
        % The first input is the offset (also the keyval for blank screen)
        offset = double(keyval(idxs(1)));
        % If the last index is single and is a blank, add it to the idxs
        if keyval(end)~=keyval(end-1) && double(keyval(end))==offset
            idxs(end+1) = length(keyval);
        end

        % Only use the indices confirmed to carry data
        all_t = ts(idxs);
        all_kv = double(keyval(idxs)) - offset; % Decode ascii chars
        stim_t = all_t(all_kv>0);
        blank_t = all_t(all_kv==0);
        
        % Find the corresponding orientations. This is hardcoded to be
        % valid for experiments where only orientation is changed, not
        % other parameters like spatial freq, contrast level and so on
        norients = max(all_kv);
        orients = (all_kv(all_kv>0)-1)*360/norients;
    case '.mat'
        % Unpack arguments not needed in the case of .inp files
        t_blank = p.Results.t_blank;
        t_delay = p.Results.t_delay;
        % Serial date, as recorded in log_file, is time in days. Conversion to
        % seconds:
        converse = 60*60*24;
        % Automatically generate name for cut file
        %[pathstr, fname, ext] = fileparts(sessionfile);
        %ext = ext(2:end);
        %cutname = [pathstr, filesep, fname, '_', ext, '.cut'];
        
        % Load stimulus log, saved in *S.logorientation* and *t0*
        S = load(log_file);
        raw_logtime = S.logorientation(:,2) * converse;
        % max_dt = max(diff(raw_logtime));
        mean_dt = mean(diff(raw_logtime));
        orients = S.logorientation(:,3);

        % Get time from start of first orientation stimulus. *logtime* is
        % start of stimuli, *logtime_blank* is start of gray screen (for
        % spontaneous rate calculation) 

        if mean_dt - floor(mean_dt) < 0.5
            duration = floor(mean_dt);
        else
            duration = ceil(mean_dt);
        end

        % Check if end time has been recorded
        [~, ncolumn] = size(S.logorientation);
        if ncolumn==3
            % t0 is approx time of button press, 
            t0 = raw_logtime(1) - t_delay - mean_dt;
            stim_t = [t_delay; raw_logtime - t0];
            blank_t = [0; stim_t(2:end) - t_blank];
            stim_t = stim_t(1:end-1);

            warning('el_phys:get_sync:ancient_warning',...
                ['End of stimulus preiods not available, baseline',...
                ' calculations done with t_blank = %.f and t_delay',...
                ' = %.f.\nThis leads to',...
                ' uncertainty\n'], t_blank, t_delay);
        elseif ncolumn==4
            t0 = S.t0 * converse;
            stim_t = raw_logtime - t0;
            % stim_t = [raw_logtime - t0; raw_logtime(end) - t0 + mean_dt];
            blank_t = [0; S.logorientation(:,4)*converse - t0];
            if verbose
                fprintf(1, ['End of stimulus preiods available, baseline',...
                    ' calculations are quite precise\n']);
            end
        end
end



function idxs = double_entries(ivec)
%DOUBLE_ENTRIES    Find double entries in ivec
% ovec = DOUBLE_ENTRIES(ivec) returns indices where ivec has at least
% two trailing entries with the same value.

idxs = [];
ientry = 1;
while ientry<=length(ivec)-1
    tmp_idx = ientry;
    val = ivec(ientry);
    nentry = 1;
    for icheck=ientry:length(ivec)
        if ivec(icheck)==val
            nentry = nentry + 1;
        else
            ientry=icheck;
            break
        end
    end
    if nentry>1
        idxs = [idxs, tmp_idx];
    end
end
