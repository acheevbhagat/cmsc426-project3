function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.
    % Initalize setup values
    confidences = {};
    c_conf = ColorConfidences.Confidences;
    c_dist = ColorConfidences.Distances;
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
    ShapeConfidences = struct('Confidences', confidences);
end
