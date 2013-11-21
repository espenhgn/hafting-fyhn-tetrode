% File example_input.m
% January 2013, Eivind Skjønsberg Norheim
% September 2013, Torkel Hafting: added saveplot

% This example uses an input file listing all session files, synch files
% and cut file. The input file can contain optional arguments for
% get_tetrode_data() as well. Write 'help get_tetrode_data' on the command
% line for details.

clear 
close('all')
addpath(['..', filesep, 'hf_tetrode'])

% The input_file is a text file containing everything needed for
% *get_tetrode_data()*
input_file = 'inputRxt1_ko_6.txt';

S = input_wrapper(input_file);

% Plot
saveplot = 1; % saves jpg files of the figures
isession = 1;
icell = 1; % icell = 0 gir alle celler
plot_spiketrains(S, isession, saveplot, 'icell', icell)
disp('finished');