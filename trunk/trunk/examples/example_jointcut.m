% File example_jointcut.m
% September 2012, Eivind Skjønsberg Norheim
% Jan 2013, ESN: Updated to use new versions of functions
% *get_tetrode_data()* and *write_csv()*

% This example uses the 'semi-old' synchronization version, where start
% times for stimulus onsets and stimulus ends were recorded. In addition an
% additional 't0' parameter for the start of the program was recorded.
%
% Give path to data file (and StimLog file if older versions of PSC has
% been used). If .cut (and .inp) files follow standard naming convention
% and lie in the same folder as sessionfiles, it's not necessary to specify
% their location.

% sync_version can be either:
%     - 'ancient', for really old cases, 
%     - 'semi-old' for cases where PsychStimController saves t0 and 4
%       columns in 'logorientation', or
%     - 'new' if proper synchronization has been performed, reported in
%       .inp file (not tested as of Jan-03-2013)

clear
close('all')

% Semi-old version, t0 and blank time recorded
datapath = ['..', filesep, 'data', filesep, 'kombinerte_filer'];
datafiles = {...
    fullfile(datapath, 'DHV_2012060232.1'),...
    fullfile(datapath, 'DHV_2012060231.1'),...
    fullfile(datapath, 'DHV_2012060230.1')...    
    };

stimfile = ['..', filesep, 'data', filesep, 'session_semi-old', filesep,...
    'StimLog-2709-1752.mat'];
for i=1:length(datafiles)
    stimfiles{i} = stimfile; %#ok
end

cutfiles = {fullfile(datapath, 'DHV_2012060230mfcom_1.cut')};
args = {'stimfiles', stimfiles, 'cutfiles', cutfiles};

% *export_args* is a cell containing extra arguments passed to
% *write_csv()* inside *get_tetrode_data()*. Since the struct contatining
% all data and the name of csv file are created inside *get_tetrode_data()*
% only other arguments are defined here.
% NOTE: These arguments are optional parameter/value pairs that can be
% excluded from the argument list.
export_csv = true;
sep = ';'; % Sign used to separate cells in spreadsheet
export_args = {'expanded', true, 'separator', sep};

% NOTE on csv files:
% If file can't be loaded in excel/libreoffice, try changing the separator
% between cells by changing *sep* to something else, for instance ':'
%
% The OSI based on fourier analysis is computed in *get_tetrode_data()* and
% this written to file in the default case. It can be omitted by setting
% the *expanded* argument in *export_args* to false.


% Append to *args*
args = [args, {'export_csv', export_csv, 'export_args', export_args}];

% Set verbose on/off and append to *args*
verbose = 1;
args = [args, {'verbose', verbose}];

% First call a function to pull out data and order it in the right way for
% later use. Write 'help get_tetrode_data' in the Command Window to get
% more info about the function.
% ***Version 0.1***: 'stimfile' is a new optional parameter/value pair,
% instead of a required argument. If no stimfile is provided, the function
% will look for an .inp file containing sync info.

S = get_tetrode_data(datafiles, args{:});

% Make a raster plot and orientation tuning curve for one of the cells.
% Change *icell* to 0 if you want to plot all cells at the same time
icell = 1;
%plot_spiketrains(S, 'icell', icell);

