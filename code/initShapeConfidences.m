function ShapeConfidences = initShapeConfidences(LocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    % Initalize setup values
    confidences = {};
  %{
    c_conf = ColorConfidences.Confidences;
    c_dist = ColorConfidences.Distances;
    seg_masks = ColorConfidences.SegmentationMasks;
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
%         
%         mask = zeros([WindowWidth WindowWidth]);
%         for i = 1:WindowWidth + 1
%             for j = 1:WindowWidth + 1
%                 dist_x = dists(i, j);
                
        % Calculate shape confidence value
        f_s = 1 - exp(-(dists.^2) ./ (sigma.^2))
        % Append to mask
%             end
%         end
        confidences{window} = f_s;
    end
    imshow(confidences{10});
    ShapeConfidences = struct('Confidences', {confidences}, 'SegmentationMasks', seg_masks);
    %}
    sigmas = {};
    cConf = ColorModels.Confidences;
    cDist = ColorModels.Distance;
    cCenter = ColorModels.LocalWindowCenter;
    for window = 1:(length(LocalWindows))
        % Find index within color confidences of the current window
        idx = find(cCenter==window);
        % Use index to set up initial values
        fc = cConf(idx);
        dx = cDist(idx);
        sigma = SigmaMin;
        % Adapt sigma value based on fc
        if fc > cutoff 
            addVal = A * (fc - fcutoff).^(R);
            sigma = sigma + addVal;
        end
        % Calculate shape confidence value
        fs = 1 - exp((-(dx.^2)) / (sigma.^2));
        % Append to confidence vector
        confidences{window} = fs;
        sigmas{window} = sigma;
        si
    end
    ShapeConfidences = struct('Confidences', {confidences}, 'Sigmas', {sigmas});
end
