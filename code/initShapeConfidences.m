function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.

    confidences = {};
     
    cConf = ColorConfidences.Confidences;
    cDist = ColorConfidences.Distance;
    cCenter = ColorConfidences.LocalWindowCenter;
    
    for window = 1:(length(LocalWindows))
        % Find index wihtin color confidences of the current window
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
        
        fs = 1 - exp((-(dx.^2)) / (sigma.^2));
        
        confidences = {confidences fs}
        
    end

    ShapeConfidences = struct('Confidences', confidences);
    
end
