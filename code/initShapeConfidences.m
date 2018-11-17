function ShapeConfidences = initShapeConfidences(LocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    % Initalize setup values
    confidences = {};
    sigmas = {};
    c_conf = ColorModels.Confidences;
    c_dist = ColorModels.Distances;
    seg_masks = ColorModels.SegmentationMasks;
    
    for window = 1:(length(LocalWindows))
        % Find index wihtin color confidences of the current window
        idx = window;
        % Use index to set up initial values
        f_c = c_conf{idx};
        dists = c_dist{idx};
        sigma = SigmaMin;
        if f_c > fcutoff 
            addVal = A * (f_c - fcutoff)^R;
            sigma = sigma + addVal;
        end
        
        % Calculate shape confidence value
        f_s = 1 - exp(-(dists.^2) ./ (sigma.^2));
        
        confidences{window} = f_s;
        sigmas{window} = sigma;
    end
    
    ShapeConfidences = struct('Confidences', {confidences}, 'SegmentationMasks', seg_masks, 'Sigmas', {sigma});
end
