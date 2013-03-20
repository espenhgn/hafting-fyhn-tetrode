function out = rfp_generate_rf(input_file, varargin)
%GENERATE_RF    Generate receptive fields from drifiting bar representations


[S, paramfiles] = input_wrapper(input_file);
 
if ~length(paramfiles)==length(S)
    msg = ['# of sessionfiles and # of parameterfiles does not match.',...
        ' Check input file for errors'];
    error('rf_pack:generate_rf:nparams', msg)
end

for isess=1:length(S)
    params = load(paramfiles{isess});
    cells_internal = rfp_receptive_field(S.session(isess).cells,...
        params, 'rf_freq', 60);
    S.session(isess).cells = cells_internal;
    clear params
end

out = S;