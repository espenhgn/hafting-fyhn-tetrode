function out = generate_rf(S)
% GENERATE_RF Generate receptive fields from drifiting bar representations

[S, paramfiles] = input_wrapper(input_file);
 
if ~length(paramfiles)==length(S)
    msg = ['# of sessionfiles and # of parameterfiles does not match.',...
        ' Check input file for errors'];
    error('rf_pack:generate_rf:nparams', msg)
end

S_internal = rf_update_data(S);


out = S_internal;