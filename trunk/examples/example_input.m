% File example_input.m
% January 2013, Eivind Skjønsberg Norheim

% This example uses an input file listing all session files, synch files
% and cut file. The input file can contain optional arguments for
% get_tetrode_data() as well. Write 'help get_tetrode_data' on the command
% line for details.

clear 
close('all')
addpath(['..', filesep, 'hf_tetrode'])

input_file = 'example_input.txt';

S = input_wrapper(input_file);