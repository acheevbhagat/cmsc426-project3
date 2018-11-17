function ShapeConfidences = initShapeConfidences(LocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    % Initalize setup values
    confidences = {};
    sigmas = {};
    c_conf = ColorModels.Confidences;
    c_dist = ColorModels.Distances;
    seg_masks = ColorModels.SegmentationMasks;
    %c_center = ColorConfidences.LocalWindowCenter;
    
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
    %imshow(confidences{10});
    ShapeConfidences = struct('Confidences', {confidences}, 'SegmentationMasks', seg_masks, 'Sigmas', {sigma});
%     sigmas = {};
%     cConf = ColorModels.Confidences;
%     cDist = ColorModels.Distance;
%     cCenter = ColorModels.LocalWindowCenter;
%     for window = 1:(length(LocalWindows))
%         % Find index within color confidences of the current window
%         idx = find(cCenter==window);
%         % Use index to set up initial values
%         fc = cConf(idx);
%         dx = cDist(idx);
%         sigma = SigmaMin;
%         % Adapt sigma value based on fc
%         if fc > cutoff 
%             addVal = A * (fc - fcutoff).^(R);
%             sigma = sigma + addVal;
%         end
%         % Calculate shape confidence value
%         fs = 1 - exp((-(dx.^2)) / (sigma.^2));
%         % Append to confidence vector
%         confidences{window} = fs;
%         sigmas{window} = sigma;
%         si
%     end
%     ShapeConfidences = struct('Confidences', {confidences}, 'Sigmas', {sigmas});
end
