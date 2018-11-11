function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    corners1 = cornermetric(IMG1, 'Harris');
    corners2 = cornermetric(IMG2, 'Harris');
    
    [features1, points1] = extractFeatures(corners1);
    [features2, points2] = extractFeatures(corners2);
    
    matchedFeatures = matchFeatures(features1, features2);
    
    matchedPoints1 = points1(matchedFeatures(:,1),:);
    matchedPoints2 = points2(matchedFeatures(:,2),:);
    
    tform = estimateGeometricTransform(matchedPoints1, matchedPoints2, 'affine');
    
    warpedWindows = [];
    
    for i = 1:length(Windows)
        warpedWindows = [warpedWindows warp(Windows(i), tform)];
    end
    
    WarpedFrame = imwarp(IMG1, tform);
    WarpedMask = warp(Mask, tform);
    WarpedLocalWindows = warpedWindows;
    
end

