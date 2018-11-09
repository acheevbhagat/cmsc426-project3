function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    % Initalize setup values
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
        mask = zeros(WindowWidth, WindowWidth);
        for i = 1:WindowWidth
            for j = 1:WindowWidth
                dxVal = dx(i, j);
                if fc > fcutoff 
                    addVal = A * (fc - fcutoff).^(R);
                    sigma = sigma + addVal;
                end
                % Calculate shape confidence value
                fs = 1 - exp((-(dx.^2)) / (sigma.^2));
                % Append to mask
                mask(i, j) = fs;
            end
        end
        confidences = {confidences mask};
    end
    ShapeConfidences = struct('Confidences', confidences);
end
