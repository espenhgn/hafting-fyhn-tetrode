clear;
close ALL;
input_file = 'example_input1267180313.txt';
S = input_wrapper(input_file);
nsess = length(S.session);
icell = 0;
saveplot = 1
for isess=1:nsess
    plot_spiketrains(S, isess, saveplot, 'icell', icell)
end

