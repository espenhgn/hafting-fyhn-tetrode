function plot_spiketrains(in_struct, varargin)
%PLOT_SPIKETRAINS   Plot spiketrains contained in struct
% PLOT_SPIKETRAINS(S) where *S* is a struct obtained from function
% receptive_fields, will plot spike trains for all orientations for all
% cells contained in *S*.
%
% Required arguments can be followed by parameter/value pairs:
%     'icell' : int (default 0)
%         Decide which cell(s) to plot. Default is to plot all cells
%         (corresponding to 0 value)
%
%
% Version 0.1, September 2012, Eivind Skjønsberg Norheim
% Version 0.2, Jan 2013, ESN, no changes
%
% See also receptive_fields.

% Argument handling
p = inputParser();
def_icell = 0;
addRequired(p, 'in_struct', @isstruct)
addOptional(p, 'icell', def_icell, @isnumeric)

parse(p, in_struct, varargin{:})

icell = p.Results.icell;

if icell==0
    icell_start = 1;
    icell_stop = length(in_struct.cells);
else
    icell_start = icell;
    icell_stop = icell;
end
    
for ic=icell_start:icell_stop
    
    cell_struct = in_struct.cells(ic);

    [norient, ntrains] = size(cell_struct.spike_trains);
    m = ceil(sqrt(norient));

    figure();    
    max_time = in_struct.duration;
    t_blank = in_struct.t_blank;
    nspikes_orient = zeros([norient,1]);

    for iorient=1:norient
        subplot(m, m, iorient)
        hold on
        for itrain=1:ntrains
            spikes = cell_struct.spike_trains{iorient, itrain};
            spont_train = cell_struct.spont_trains{iorient, itrain};
            spike_train = [spont_train; spikes];
            nspikes = length(spikes);
            nspikes_orient(iorient) = nspikes_orient(iorient) + nspikes;
            spike_shape = [itrain-0.9, itrain-0.1];
            for ispike=1:nspikes
                spike_time = spike_train(ispike);
                plot([spike_time, spike_time], spike_shape, 'k-',...
                    'LineWidth', 1)
                if spike_time>max_time
                    max_time = spike_time;
                end
            end
        end
        hold off
        axis([-t_blank, max_time - t_blank, 0, itrain])
        if iorient~=m*(m-1)+1
            set(gca, 'XTick', [])
        else
            xlabel('Time, t (s)')
            ylabel('Trials')
        end        
        set(gca, 'YTick', [])
        titlestr = sprintf('%.0f\\circ', in_struct.orientations(iorient));
        title(titlestr)
    end
    
    % Making main title
    ha = axes('Position', [0 0 1 1], 'Xlim', [0 1], 'Ylim', [0 1], ...
        'Box', 'off', 'Visible', 'off', 'Units', 'normalized',...
        'clipping', 'off');
    titlestr = sprintf('Raster for cell # %g', ic);
    text(0.5, 1, titlestr, 'HorizontalAlignment',...
        'center','VerticalAlignment', 'top');
    % Making orientation tuning plot
    spikes_pr_sec = nspikes_orient/(in_struct.duration*ntrains);

    figure();
    plot(in_struct.orientations, spikes_pr_sec, 'k-o')
    axis([0, in_struct.orientations(end), 0, 1.1*max(spikes_pr_sec)])
    xlabel('direction (degrees)');
    ylabel('spikerate (s^{-1})');
    titlestr = sprintf('Orientation tuning curve, cell # %g, OSI = %.3f',...
        ic, cell_struct.OSI);
    title(titlestr)

end

function out = isPosInt(n)

out = 0;

if isinteger(n(1))
    if n(1)>=0
        out=1;
    end
end
