function out = rfp_receptive_field(cell_struct, params, varargin)
%RFP_RECEPTIVE_FIELD    Generate receptive fields from drifting bar data
%
% 

% Argument handling
p = inputParser();
def_freq = 60;

addRequired(p, 'cell_struct', @isstruct)
addRequired(p, 'params', @isstruct)
addParamValue(p, 'rf_freq', def_freq, @isnumeric)

parse(p, cell_struct, params, varargin{:})
freq = p.Results.rf_freq;

% Done with argument handling

norients = length(cell_struct(1).nspikes);

% This is picked up from PsychStimController
ScreenSizeDegX = params.SizeX*atan(1/params.ScreenDist)*180/pi;
degPerPix = ScreenSizeDegX/params.PixelsX;
sizeLut = 256;

for icell=1:length(cell_struct)
    rf = zeros(params.PixelsY, params.PixelsX);
    
    for iorient=1:norients
        s_trains_cell = cell_struct(icell).spike_trains(iorient, :);        
        for itrain=1:length(s_trains_cell)
            spike_train = s_trains_cell{itrain};
            duration = cell_struct(icell).stim_durations(iorient, itrain);
            %nframes = floor(duration*freq)
            [img cl] = rfp_generateBars_lut(params.orient(iorient),...
                params.freq(iorient), params.speed(iorient),...
                params.contrast(iorient), params.length(iorient),...
                params.positionX, duration, degPerPix,...
                params.PixelsX, params.PixelsY, freq, BlackIndex(0),...
                WhiteIndex(0), sizeLut);
            %t = 0:(1/freq):duration-1/freq;
            
            rf = rf + rfp_clut2mask(img, cl, freq, spike_train);
        end
        fprintf(1, 'Cell #%d, orient #%d: cumulative idx_sum=%d\n', icell, iorient,...
            sum(sum(rf)));
    end
    cell_struct(icell).receptive_field = rf;
end

out = cell_struct;
