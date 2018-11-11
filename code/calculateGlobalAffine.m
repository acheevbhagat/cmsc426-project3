function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    img1 = rgb2gray(IMG1);
    img2 = rgb2gray(IMG2);
    corners1 = detectHarrisFeatures(img1);
    corners2 = detectHarrisFeatures(img2);
    
    [features1, points1] = extractFeatures(img1, corners1);
    [features2, points2] = extractFeatures(img2, corners2);
    
    matchedFeatures = matchFeatures(features1, features2);
    
    matchedPoints1 = points1(matchedFeatures(:,1),:);
    matchedPoints2 = points2(matchedFeatures(:,2),:);
    
    tform = estimateGeometricTransform(matchedPoints1, matchedPoints2, 'affine');
    
    warpedWindows = {};
    
    for i = 1:length(Windows)
        warpedWindows = [warpedWindows imwarp(Windows(i), tform)];
    end
    
    WarpedFrame = imwarp(IMG1, tform);
    WarpedMask = imwarp(Mask, tform);
    WarpedMaskOutline = bwperim(WarpedMask, 4);
    WarpedLocalWindows = warpedWindows;
    
end

