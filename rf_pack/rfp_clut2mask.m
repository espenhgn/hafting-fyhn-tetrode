function pixel_mask = rfp_clut2mask(img, cl, freq, spike_train)
% RFP_CLUT2MASK    Find binary mask from image and clut

pixel_mask = zeros(size(img));
[~, ~, nframes] = size(cl);

for ispike=1:length(spike_train)
    t_spike = spike_train(ispike);
    clut_idx = floor(t_spike*freq) + 1;
    if clut_idx<=nframes
        try
            idxs_on = find(cl(:,1,clut_idx)==max(cl(:,1,clut_idx)));
        catch ME
            t_spike
            clut_idx
            size(cl)
            ME.message
            break
        end
        single_mask = ones(size(img)) .* (img>=idxs_on(1)) .* (img<=idxs_on(end));
        pixel_mask = pixel_mask + single_mask;
    end
end
