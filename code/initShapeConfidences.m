function ShapeConfidences = initShapeConfidences(LocalWindows, ColorConfidences, WindowWidth, SigmaMin, A, fcutoff, R)
% INITSHAPECONFIDENCES Initialize shape confidences.  ShapeConfidences is a struct you should define yourself.

    confidences = {}
     
    cConf = ColorConfidences.Confidences;
    cDist = ColorConfidences.Distance;
    cCenter = ColorConfidences.LocalWindowCenter;
    
    for window = 1:(length(LocalWindows))
    end

end
