function plot_spiketrains(in_struct, isession, varargin)
%PLOT_SPIKETRAINS   Plot spiketrains contained in struct
% PLOT_SPIKETRAINS(S, isession) plots spike trains for all orientations for all cells
% contained in *S.session(isession)*. *S* is a struct obtained from
% function receptive_fields, and isession is index of the session of
% interest.
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
%def_session = 1;
addRequired(p, 'in_struct', @isstruct)
addRequired(p, 'isession', @isnumeric)
addParamValue(p, 'icell', def_icell, @isnumeric)


parse(p, in_struct, isession, varargin{:})

icell = p.Results.icell;

session = in_struct.session(isession);

if icell==0
    icell_start = 1;
    icell_stop = length(session.cells);
else
    icell_start = icell;
    icell_stop = icell;
end
    
for ic=icell_start:icell_stop
    
    cell_struct = session.cells(ic);

    [norient, ntrains] = size(cell_struct.spike_trains);
    m = ceil(sqrt(norient));

    figure();
    max_time = max(max(cell_struct.stim_durations));
    max_t_blank = mean(mean(cell_struct.blank_durations));
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
        grid on
        hold off
        axis([-max_t_blank, max_time, 0, itrain])
        if iorient~=m*(m-1)+1
            set(gca, 'XTick', [])
        else
            xlabel('Time, t (s)')
            ylabel('Trials')
        end        
        set(gca, 'YTick', [])
        titlestr = sprintf('%.0f\\circ', session.orientations(iorient)); % edit
        title(titlestr)
    end
    
    % Making main title
    ha = axes('Position', [0 0 1 1], 'Xlim', [0 1], 'Ylim', [0 1], ...
        'Box', 'off', 'Visible', 'off', 'Units', 'normalized',...
        'clipping', 'off'); %#ok
    titlestr = sprintf('Session # %g: Raster for cell # %g', isession, ic);
    text(0.5, 1, titlestr, 'HorizontalAlignment',...
        'center','VerticalAlignment', 'top');
    % Making orientation tuning plot
    spike_rates = mean(cell_struct.spike_rates, 2);

    figure();
    plot(session.orientations, spike_rates, 'k-o') % edit
    axis([0, session.orientations(end), 0, 1.1*max(spike_rates)]) % edit
    xlabel('direction (degrees)');
    ylabel('spikerate (s^{-1})');
    titlestr = sprintf(['Session # %g: Orientation tuning curve,',...
        ' cell # %g, OSI = %.3f'],...
        isession, ic, cell_struct.OSI);
    title(titlestr)
end

function out = isPosInt(n) %#ok

out = 0;

if isinteger(n(1))
    if n(1)>=0
        out=1;
    end
end
