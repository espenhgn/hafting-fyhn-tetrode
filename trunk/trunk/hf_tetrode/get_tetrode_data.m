function out = get_tetrode_data(sessionfiles, varargin)
%GET_TETRODE_DATA    Load spiketimes with corresponding cut and stimulus.
% S = GET_TETRODE_DATA(sessionfiles) where *sessionfiles* is cell array
% containing (relative or absolute)  paths to files. Files in *sessionfile*
% are tetrode data files created by AXONA. Additional cut files must be 
% present in the same folder as *sessionfiles*, with standard name
% convention from AXONA
%
% The file names can be followed by parameter/value pairs:
%     'cutfiles' : cell (default {})
%         array containing paths to cut files when naming
%         convention is not followed, or if joint cut file is created. 
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
%     'export_csv' : boolean or int (default false)
%         Set this to true or 1 to write some results to a csv file that
%         can be loaded in libreoffice/excel. You can do this afterwards as
%         well by calling the function write_csv().
%     'export_args' : cell (default {})
%         Parameter/Value pairs passed to *write_csv()* (if *export_csv* is set to
%         true).
%         NOTE: Required arguments for *write_csv()* (*in_struct* and
%         *ofile*) are created and passed internally, so arguments in
%         *export_args* should be the parameter/value pairs described in
%         *write_csv()*. for documentation write 'help write_csv'on the
%         command line, or press the link at the bottom of this page.
%     'verbose' : boolean or int (default false)
%         Toggles messages on/off. false means no messages, true means 
%         print messages to standard output. Integer values are also
%         accepted: 0 and negative values means false, positive values
%         means true.
%
%
% Struct *S* has fields *cells*, *orientations*, *duration*, 
% *mean_duration*, *logtime*, *logtime_blank* and *t_blank*. 
% *S.cells* is a struct with length *ncells*, where each entry in the
% struct has field *spike_trains*, which is a cell array with
% *norientations* rows and *ntrials* columns. Each element in
% *spike_trains* is an array containg the spike times for the given
% orientation and trial. *S.orientations* is an array with stimulus
% orientations, *S.duration* is the rounded off stimulus duration and
% *mean_duration* is the mean of the real stimulus durations pulled out of
% *stimfiles*. *logtime* is the time vector for stimulus onset,
% *logtime_blank* is the time vector for gray screen onset and t_blank is
% the duration of gray screen.
%
% Version 0.1, September 2012, Eivind Skjønsberg Norheim
% Version 0.2, Jan 2013, ESN:
%     Changes:
%         - Functionality for using sync info from .inp files.
%         - 'stimfile' is a new optional parameter/value pair, instead of a
%           required argument. If no stimfile is provided, the function
%           will look for an .inp file containing sync info.
%
% See also write_csv


% Argument handling
p = inputParser();
def_cutfiles = {};
def_stimfiles = {};
def_verbose = false;
def_tblank = zeros(size(sessionfiles));
def_tdelay = zeros(size(sessionfiles));
def_writecsv = false;
def_csvargs = {};

addRequired(p, 'sessionfiles', @is_cellorchar)
addParamValue(p, 'cutfiles', def_cutfiles, @is_cellorchar)
addParamValue(p, 'stimfiles', def_stimfiles, @is_cellorchar)
addParamValue(p, 't_delay', def_tdelay, @isnumeric)
addParamValue(p, 't_blank', def_tblank, @isnumeric)
addParamValue(p, 'verbose', def_verbose, @is_my_logical)
addParamValue(p, 'export_csv', def_writecsv, @is_my_logical)
addParamValue(p, 'export_args', def_csvargs, @iscell)

parse(p, sessionfiles, varargin{:})

stimfiles = p.Results.stimfiles;
cutfiles = p.Results.cutfiles;
t_blank = p.Results.t_blank;
t_delay = p.Results.t_delay;
verbose = p.Results.verbose;
% Done with argument handling

% Initialize and load files
out = struct();
nsess = length(sessionfiles);
data_S = load_tetrode_files(sessionfiles, 'cutfiles', cutfiles, ...
    'stimfiles', stimfiles, 't_delay', t_delay, ...
    't_blank', t_blank, 'verbose', verbose);

% Loop through all session files, pick the right data from the load
% structure and put it into the out structure
for isess=1:nsess
    timestamps = data_S.timestamps{isess};
    logtime = data_S.sync_data(isess).logtime;
    logtime_blank = data_S.sync_data(isess).logtime_blank;
    orients = data_S.sync_data(isess).orients;
    clust = data_S.clusts{isess};
    ncells = max(clust);

    % find durations and mean duration of stimulus
    dt = diff([logtime, logtime_blank(2:end)], 1, 2);
    % dt_blank = diff([logtime_blank(1:end-1), logtime], 1, 2);
    mean_dt = mean(dt);
    % mean_dt_blank = mean(dt_blank);

    % Get unique orientations, and number of repetitions of each orientation
    unique_orients = unique(orients);
    norients = length(unique_orients);
    trials_pr_orientation = hist(orients, unique_orients);
    cells = struct();

    if verbose
        fprintf(1, 'Orientations are:\n');
        for iorient=1:length(unique_orients)
            fprintf(1, '%g, # of trials: %g\n', unique_orients(iorient),...
                trials_pr_orientation(iorient));
        end
    end

    % Loop through cells from cut file and assign timestamps to the right
    % cells 
    for icell=1:ncells
        % Update *cells* struct with empty cells for evoked spike trains
        % and spontatneous spike trains
        cells(icell,1).spike_trains = cell(norients,...
            max(trials_pr_orientation));        
        cells(icell,1).spont_trains = ...
            cell(size(cells(icell,1).spike_trains));
        % spike rates are just numbers, so use an array
        cells(icell,1).spike_rates = zeros(norients,...
            max(trials_pr_orientation));

        % Use a counter pr orientation to keep track of which spike trains
        % have been counted
        nspike_trains = ones(length(unique_orients),1);
        nspikes = zeros(size(nspike_trains));
        spont_spikes = 0;
        spont_time = 0;
        % pick timestamps for the right cell
        cell_timestamps = timestamps(clust==icell);

        for itimeint=1:length(logtime)
            % pick start and end times from *logtime* and *logtime_blank*.
            % *logtime_blank* has one more entry than logtime to ensure
            % that we know when the last stimulus period ends.
            % ***stim periods start immediately after blank periods***
            tstart = logtime(itimeint);
            tend = logtime_blank(itimeint+1);
            tstart_base = logtime_blank(itimeint);
            tend_base = tstart;

            % Find the corresponding orientation
            orient = orients(itimeint);
            iorient = find(unique_orients==orient);

            % Append timestamps within stimulus and grey period to *cells*
            % structure
            spike_train = cell_timestamps(cell_timestamps>=tstart &...
                cell_timestamps<tend) - tstart;
            cells(icell).spike_trains{iorient, nspike_trains(iorient)} = ...
                spike_train;
            cells(icell).spont_trains{iorient, nspike_trains(iorient)} = ...
                cell_timestamps(cell_timestamps>=tstart_base&...
                        cell_timestamps<tend_base) - tstart;
            cells(icell).spike_rates(iorient, nspike_trains(iorient)) = ...
                length(spike_train)/(tend - tstart);

            % Update book keeping variables
            spont_spikes = spont_spikes + ...
                length(cells(icell).spont_trains{iorient,...
                                        nspike_trains(iorient)});        
            spont_time = spont_time + (tend_base - tstart_base);
            nspike_trains(iorient) = nspike_trains(iorient) + 1;
            nspikes(iorient) = nspikes(iorient) + length(spike_train);        
        end
        
        cells(icell).spont_rate = spont_spikes/spont_time;
        cells(icell).nspikes = nspikes;
        cells(icell).evoked_rate = mean(mean(cells(icell).spike_rates));
        [cells(icell).OSI, cells(icell).fOSI] = calc_OSI(nspikes);
    end
    % Update out structure with the last cells struct
    % ***Redundant save values, update plot_spiketrains.m to accept a new
    % struct***

    if length(cutfiles)==1 && nsess~=1
        cutfile = GetFullPath(cutfiles{1});
        j_cut = true;
    else
        cutfile = GetFullPath(data_S.cutfiles{isess});
        j_cut = false;
    end
    
    stimfile = GetFullPath(data_S.stimfiles{isess});
    sessionfile = GetFullPath(sessionfiles{isess});
    
    sess_out = struct('cells', cells, ...
        'orientations', unique_orients, ...
        'duration', mean_dt, 'mean_duration', mean_dt, ...
        'logtime', logtime, 'logtime_blank', logtime_blank, ...
        't_blank', t_blank, 't_delay', t_delay, ...
        'sessionfile', sessionfile, ...
        'cutfile', cutfile, 'stimfile', stimfile, 'joint_cut', j_cut, ...
        'trials_pr_orientation', trials_pr_orientation);
    
    out.session(isess) = sess_out;
    
    if p.Results.export_csv
        [pathstr, fname, ext] = fileparts(sessionfile);
        csv_name = fullfile(pathstr, [fname, '_', ext(2:end), '_summary.csv']);
        write_csv(out.session(isess), csv_name, p.Results.export_args{:});
    end
end

if verbose
    fprintf(1, 'Done loading tetrode data\n');
end


function out = is_cellorchar(n)
%ISMYLOGICAL    checks if value is cell or char.

out = 0;
if iscell(n) || ischar(n)
    out=1;
end
