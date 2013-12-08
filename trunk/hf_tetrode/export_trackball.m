function outputdata = export_trackball(session, setfile, trackball)
%EXPORT_TRACKBALL Writes a cvs file that libreoffice/excel can open 
%
% outputdata = export_trackball(session, setfile, trackball)
%
% arguments:
%   'session' : session object from input_wrapper(), aka S.session
%   'setfile' : setfile object from import_setfile()
%   'trackfile' : trackfile object from import_trackball()
%
% outputdata = 
% 
%        trackball: [1x1 struct] = 
%                 speeds: [96x1 double]
%                     Vx: [96x1 double]
%                     Vy: [96x1 double]
%                  dWZdt: [96x1 double] 
%             data: [480x10 double]
%           header: {5x2 cell}
%     columntitles: {1x10 cell}

%% container
outputdata = {};

%% get mean speeds etc. during each trial
durations = session.cells(1).stim_durations;
speeds = zeros(length(session.logtime), 1);
Vx = zeros(length(speeds), 1);
Vy = zeros(length(speeds), 1);
dWZdt = zeros(length(speeds), 1);
for ii = 1:length(speeds)
    tstart = session.logtime(ii);
    tend = tstart + durations(ii);
    tinds = (trackball.time(1:end-1) >= tstart) & ...
            (trackball.time(1:end-1) <= tend);
    speeds(ii) = mean(trackball.speed(tinds));
    Vx(ii) = mean(trackball.Vx(tinds));
    Vy(ii) = mean(trackball.Vy(tinds));
    dWZdt(ii) = mean(trackball.domegaZdt(tinds));
end

%% create struct with these data
outputdata.trackball = {};
outputdata.trackball.speeds = speeds;
outputdata.trackball.Vx = Vx;
outputdata.trackball.Vy = Vy;
outputdata.trackball.dWZdt = dWZdt;


%% prepare output data array
numrows = length(session.data_S.sync_data.orients);
columntitles = {'cell #', 'orientation', ...
    'logtime (s)', 'stim_durations (s)', ...
    'spike_rates (Hz)', 'nspikes', ...
    'speed (m/s)', 'Vx (m/s)', 'Vy (m/s)', 'dWdt (rad/s)'};
output = [];
for ii = 1:length(session.cells)
    data = zeros(numrows, 10);
    data(:, 1) = repmat(ii, numrows, 1);
    data(:, 2) = session.data_S.sync_data.orients;
    data(:, 3) = session.logtime;
    data(:, 4) = reshape(session.cells(ii).stim_durations, 1, []);
    data(:, 5) = reshape(session.cells(ii).spike_rates, 1, []);
    data(:, 6) = cellfun(@(x) numel(x), ...
                        reshape(session.cells(ii).spike_trains, 1, []));
    data(:, 7) = speeds;
    data(:, 8) = Vx;
    data(:, 9) = Vy;
    data(:, 10) = dWZdt;
    
    output = [output; data];
end

%% attach output
outputdata.data = output;

%% file header info
header = [{'Session' session.sessionfile}; ...
    {'Cut', session.cutfile}; ...
    {'Stim synch', session.stimfile}; ...
    {'Set', setfile.filename}; ...
    {'Trackball', trackball.filename}];

%% attach to output
outputdata.header = header;
outputdata.columntitles = columntitles;

%% open and write to file
outputf = strsplit(session.cutfile, '.cut');
outputf = strjoin([outputf(1) '_summary_trackball.csv'], '');
fprintf(1, ['writing output file ', outputf, '\n']);
fid = fopen(outputf, 'w');
DLM = ';';
for ii = 1:length(header)
    fprintf(fid, [strjoin(header(ii,:), ';'), '\n']);
end
fprintf(fid, '\n');
fprintf(fid, [strjoin(columntitles, ';'), '\n']);
fclose(fid);

dlmwrite(outputf, outputdata.data, 'delimiter', DLM, '-append');



