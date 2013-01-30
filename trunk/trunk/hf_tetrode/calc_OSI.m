function [OSI, fourier_OSI] = calc_OSI(nspikes, varargin)
%CALC_OSI   calculate the orientation selectivity index (OSI)
% OSI = CALC_OSI(nspikes), where *nspikes* is a 1D array containing total
% spike counts for one orientation, calculates the OSI as
%
% OSI = (R_pref - R_ortho)/(R_pref + R_ortho)   ,
%
% where R_pref is the direction with most spikes, and R_ortho is 'the'
% orthogonal direction.

% Input parsing
p = inputParser();
default_mode = 'fourier';
expected_modes = {'naive', 'fourier'};
addRequired(p, 'nspikes', @isnumeric);
addParamValue(p,'shape',default_mode,...
    @(x) any(validatestring(x,expected_modes)));
parse(p, nspikes, varargin{:})

% For the moment both OSIs are found, functionality is ready for choosing
% between different modes
%mode = p.results.mode;

% find OSI based on different modes
norients = length(nspikes);
[r_pref, i_pref] = max(nspikes);

if mod(norients, 2)
    % Number of orientations is odd, let's not go there for now
    disp('Not implemented yet, odd number of orientations.')
    r_pref = 1;
    r_ortho = 1;
else
    if mod(norients, 4)
        % Orthogonal direction is not measured, use mean of two closest
        % orientations?
    else
        % Orthogonal orientation is measured
        i_orthos = [i_pref+1/4*norients,...
            i_pref + 3/4*norients];
        for i=1:2
            if i_orthos(i)>norients
                i_orthos(i) = i_orthos(i)-norients;
            end
        end
        r_ortho = mean(nspikes(i_orthos));
    end
end

OSI = (r_pref - r_ortho)/(r_pref + r_ortho);

Nspikes = fft(nspikes);
A0 = abs(Nspikes(1))/norients;
A2 = abs(Nspikes(3))/norients;
fourier_OSI = A2/(A0+A2);