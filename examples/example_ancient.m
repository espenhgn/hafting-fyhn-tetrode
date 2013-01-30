% File example_ancient.m
% September 2012, Eivind Skjønsberg Norheim
% Jan 2013, ESN: Updated to use new versions of functions
% *get_tetrode_data()* and *write_csv()*

% This example uses the 'ancient' synchrony version, where only end times
% of stimuli were recorded. t_blank and t_delay must be provided.
%
% Give path to data file (and StimLog file if older versions of PSC has
% been used). If .cut (and .inp) files follow standard naming convention
% and lie in the same folder as sessionfiles, it's not necessary to specify
% their location.

clear
close('all')

% Ancient version, no t0, no blank time recorded
stimfile = {['..', filesep, 'data', filesep, 'session_ancient', filesep,...
    'StimLog-0907-1037.mat']};
datafile = {['..', filesep, 'data', filesep, 'session_ancient', filesep,...
    'r1040090712s3.8']};
t_delay = 1;
t_blank = 0.75;
args = {'stimfiles', stimfile, 't_delay', t_delay,...
    't_blank', t_blank};

% Set arguments for writing csv file with results
oname = 'summary.csv';
opath = fullfile('..', 'outfiles', '');
ofile = fullfile(opath, oname);

% *export_args* is a cell containing extra arguments passed to
% *write_csv()* inside *get_tetrode_data()*. Since the struct contatining
% all data is created inside *get_tetrode_data()* only other arguments are
% defined here.
% NOTE: *ofile* is a required argument, others are optional parameter/value
% pairs that can be excluded from the argument list.
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

S = get_tetrode_data(datafile, args{:});

% Make a raster plot and orientation tuning curve for one of the cells.
% Change *icell* to 0 if you want to plot all cells at the same time
icell = 1;
%plot_spiketrains(S, 'icell', icell);

