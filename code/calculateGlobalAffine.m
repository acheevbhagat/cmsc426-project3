function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows,Width)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    
    img1 = rgb2gray(IMG1);
    img2 = rgb2gray(IMG2);
    corners1 = detectHarrisFeatures(img1);
    corners2 = detectHarrisFeatures(img2);
    [features1, points1] = extractFeatures(img1, corners1);
    [features2, points2] = extractFeatures(img2, corners2);
    
    matchedFeatures = matchFeatures(features1, features2);
    
    matched_pts1 = points1(matchedFeatures(:,1),:);
    matched_pts2 = points2(matchedFeatures(:,2),:);

%{
    IMG1_gray = single(rgb2gray(IMG1));
    IMG2_gray = single(rgb2gray(IMG2));
    
    [f1, d1] = vl_sift(IMG1_gray);
    [f2, d2] = vl_sift(IMG2_gray);
    
    matches = vl_ubcmatch(d1, d2);
    matched_pts1 = round(f1(1:2, matches(1, :))');
    matched_pts2 = round(f2(1:2, matches(2, :))');
%}
    showMatchedFeatures(IMG1, IMG2, matched_pts1, matched_pts2, 'montage');
    
    tform = estimateGeometricTransform(matched_pts1, matched_pts2, 'affine');
    
    WarpedFrame = imwarp(IMG1, tform, 'OutputView', imref2d(size(IMG1)));
    WarpedMask = imwarp(Mask, tform, 'OutputView', imref2d(size(IMG1)));
    
    WarpedMaskOutline = bwperim(WarpedMask, 4);
    WarpedLocalWindows = round(transformPointsForward(tform, Windows));
    %[~, WarpedLocalWindows] = initLocalWindows(WarpedFrame, WarpedMask, length(Windows), Width, true); 

end

