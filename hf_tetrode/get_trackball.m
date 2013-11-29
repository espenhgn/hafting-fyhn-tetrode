function [trackball] = get_trackball(setfile, DPI, diameter, dt)
%GET_TRACKBALL Assess and read trackball file corresponding to axona
%setfile, by matching timestamps
%
%Usage
%   trackball = get_trackball(setfile, DPI, diameter)
%arguments:
%   setfile = import_setfile(setfilename)
%   DPI = 800, DPI setting of trackball mice
%   diameter = 0.198, (m) diameter of trackball
%   dt = 0.01, (s) rate for resampling to regular grid of position data
%
%output structure:
% trackball = 
%      filename: 'trackball_20131013_1014.txt'      file
%            ID: [190051x1 double]      mouse ID, either 1 or 2
%        DYTime: {190051x1 cell}        timestamps column 1
%       DYTime1: {190051x1 cell}        timestamps column 2, not used
%            DX: [190051x1 double]      change in X-readout
%            DY: [190051x1 double]      change in Y-readout
%           DPI: 3600                   DPI setting of computer mice
%           DPm: 3.1496e+04             Dots Per meter
%      diameter: 0.1980                 trackball diameter in units of (m)
%            dt: 0.0100                 time res of regular grid readouts
%            T1: [96779x1 double]       time stamps mouse 1 in (s)
%            T2: [93268x1 double]       time stamps mouse 2 in (s)
%            X1: [96779x1 double]       x-position mouse 1 in (m)
%            Y1: [96779x1 double]       y-position mouse 1 in (m)
%            X2: [93268x1 double]       x-position mouse 2 in (m)
%            Y2: [93268x1 double]       y-position mouse 2 in (m)
%          time: [1x41101 double]       resample time vector
%          ipX1: [1x41101 double]       interpolated x-position mouse 1 (m) 
%          ipY1: [1x41101 double]       interpolated y-position mouse 1 (m)
%          ipX2: [1x41101 double]       interpolated x-position mouse 2 (m)
%          ipY2: [1x41101 double]       interpolated y-position mouse 2 (m)
%        omegaZ: [1x41101 double]       angular position top-down axis (rad)
%         dXdt1: [1x41100 double]       tangential x-velocity mouse 1 (m/s) 
%         dYdt1: [1x41100 double]       tangential y-velocity mouse 1 (m/s)
%         dXdt2: [1x41100 double]       tangential x-velocity mouse 2 (m/s)
%         dYdt2: [1x41100 double]       tangential y-velocity mouse 2 (m/s)
%     domegaZdt: [1x41100 double]       angular velocity (rad/s)
%         speed: [1x41100 double]       scalar speed from tang. components
%
%TODO: 
%   filtering of position data, though they seem quite clean already
%   confirm which axis is forward, sideways, top-down rotation
%   confirm which USB mouse 1 is on back, 2 on side of trackball

%% get date of trial from setfile:
trial_date = strsplit(setfile.trial_date, ', ');
trial_date = strjoin(strsplit(trial_date{end}), '-');
trial_date = datestr(datenum(trial_date), 'yyyymmdd');


%% assess corresponding trackball file:
%we now that trackball must start at or before the timestamp of .set-file
%so we try and open file created same minute
trackballfile = ['trackball_', trial_date, '_', ...
    datestr(setfile.trial_time, 'HHMM'), '.txt'];
fid = fopen(trackballfile, 'r');
fclose(fid);
%if fid was nonexistent, roll back in time, one minute at the time, up to
%one hour
if fid == -1
    for ii = 1:60
        HHMM = datestr(datenum(setfile.trial_time, 'HHMM') ...
            - datenum(['', ii], 'MM'), 'HHMM');
        trackballfile = ['trackball_', trial_date, '_', HHMM, '.txt'];
        fid = fopen(['trackball_', trial_date, '_', HHMM, '.txt'], 'r');
        fclose(fid);
        if fid ~= -1
            break
        end
    end
end
%with some luck, trackballfile should now be the one corresponding to .set
trackballfile = GetFullPath(trackballfile);

%% container for trackball data, set some fields
trackball = {};
trackball.filename = trackballfile;
[trackball.ID, trackball.DYTime, trackball.DYTime1, ...
    trackball.DX, trackball.DY] = import_trackball(trackballfile);

%% USB mouse 
%Logitech spec sheet say optical resolution is up to 3600 DPI,
%but hardware adjustable, Christina think it was set to 800 DPI. 
%Logitech software however, say as low as 200 is possible (software
%tricks?)
trackball.DPI = DPI; %3600
trackball.DPm = trackball.DPI / 2.54 * 100;
trackball.diameter = diameter; %m, measured ball diameter

%dt for regularly spaced, resampled positions and velocity components 
trackball.dt = dt;

%% find valid time stamps
%separate timestamps and coords for each out of two computer mice
%this only need to be semi accurate, as setfile is only reporting times 
%to the nearest second. We will use first occurrence of the same sec.
skipfac = 10; % testing every skipfac timestamp to speed up things
HHMMSS = datestr(trackball.DYTime(1:skipfac:end), 'HH:MM:SS');
startt = datestr(setfile.trial_time, 'HH:MM:SS'); 
endt = datestr(double(datenum(setfile.trial_time, 'HH:MM:SS')) ...
    + double(datenum(num2str(setfile.duration), 'SS')), 'HH:MM:SS');

%% find indices that span recording duration:
startt_ind = find(datenum(HHMMSS) == datenum(startt), 1, 'first');
endt_ind = find(datenum(HHMMSS) == datenum(endt), 1, 'last');
inds = startt_ind*skipfac:endt_ind*skipfac;

%% discard nonneeded trackball datas
trackball.ID = trackball.ID(inds);
trackball.DYTime = trackball.DYTime(inds);
trackball.DYTime1 = trackball.DYTime1(inds);
trackball.DX = trackball.DX(inds);
trackball.DY = trackball.DY(inds);

%for each of the two mice, compute the cumsum of DX and DY (i.e.,
%position).
%To ease interpolation below, times to seconds starting at t=0 for 
%first valid time point in trackball.DYTime
first_tstamp = datenum(trackball.DYTime(1));
inds1 = trackball.ID == 1;
inds2 = trackball.ID == 2;

%% corresponding timestamps, convert to unit of seconds
trackball.T1 = (datenum(trackball.DYTime(inds1)) - first_tstamp)*24*60*60;
trackball.T2 = (datenum(trackball.DYTime(inds2)) - first_tstamp)*24*60*60;

%% position data
% position is cumsum of change in each direction
trackball.X1 = cumsum(trackball.DX(inds1));
trackball.Y1 = cumsum(trackball.DY(inds1));
trackball.X2 = cumsum(trackball.DX(inds2));
trackball.Y2 = cumsum(trackball.DY(inds2));

%first coordinate is position 0
trackball.X1 = trackball.X1 - trackball.X1(1);
trackball.Y1 = trackball.Y1 - trackball.Y1(1);
trackball.X2 = trackball.X2 - trackball.X2(1);
trackball.Y2 = trackball.Y2 - trackball.Y2(1);

%filter out repeated timestamps:
badtinds1 = find(diff(trackball.T1) == 0);
badtinds2 = find(diff(trackball.T2) == 0);

%% correct for bad time inds
trackball.T1(badtinds1) = [];
trackball.T2(badtinds2) = [];
trackball.X1(badtinds1) = [];
trackball.Y1(badtinds1) = [];
trackball.X2(badtinds2) = [];
trackball.Y2(badtinds2) = [];

%% convert to units of m, compute velocity and z-axis rotation (top-down)
trackball.X1 = trackball.X1 / trackball.DPm;
trackball.Y1 = trackball.Y1 / trackball.DPm;
trackball.X2 = trackball.X2 / trackball.DPm;
trackball.Y2 = trackball.Y2 / trackball.DPm;

%% resample signals to timeres dt s.
trackball.time = 0:trackball.dt:setfile.duration;
trackball.ipX1 = ...
    interp1(trackball.T1, trackball.X1, trackball.time, 'linear');
trackball.ipY1 = ...
    interp1(trackball.T1, trackball.Y1, trackball.time, 'linear');
trackball.ipX2 = ...
    interp1(trackball.T2, trackball.X2, trackball.time, 'linear');
trackball.ipY2 = ...
    interp1(trackball.T2, trackball.Y2, trackball.time, 'linear');

%% compute the velocity vector components
trackball.dXdt1 = diff(trackball.ipX1) / trackball.dt;
trackball.dYdt1 = diff(trackball.ipY1) / trackball.dt;
trackball.dXdt2 = diff(trackball.ipX2) / trackball.dt;
trackball.dYdt2 = diff(trackball.ipY2) / trackball.dt;

%% angular position around top-down axis, from the mean of the Y-components
trackball.omegaZ = (trackball.ipY1 + trackball.ipY2) / 2 ...
    * 1. / (pi*trackball.diameter);
trackball.domegaZdt = diff(trackball.omegaZ) / trackball.dt;

%% compute the scalar speed
trackball.speed = sqrt(trackball.dXdt2.^2 + trackball.dXdt1.^2);
