function [WarpedFrame, WarpedMask, WarpedMaskOutline, WarpedLocalWindows] = calculateGlobalAffine(IMG1,IMG2,Mask,Windows,Width)
% CALCULATEGLOBALAFFINE: finds affine transform between two frames, and applies it to frame1, the mask, and local windows.
    
%     img1 = rgb2gray(IMG1);
%     img2 = rgb2gray(IMG2);
%     corners1 = detectHarrisFeatures(img1);
%     corners2 = detectHarrisFeatures(img2);
%     [features1, points1] = extractFeatures(img1, corners1);
%     [features2, points2] = extractFeatures(img2, corners2);
%     
%     matchedFeatures = matchFeatures(features1, features2);
%     
%     matched_pts1 = points1(matchedFeatures(:,1),:);
%     matched_pts2 = points2(matchedFeatures(:,2),:);
    
    y_min = min(Windows(:, 2));
    y_max = max(Windows(:, 2));
    x_min = min(Windows(:, 1));
    x_max = max(Windows(:, 1));
    
    IMG1_gray = single(rgb2gray(IMG1));
    IMG1_gray = IMG1_gray(y_min - Width/2:y_max + Width/2, x_min - Width/2:x_max + Width/2);
    IMG2_gray = single(rgb2gray(IMG2));
    IMG2_gray = IMG2_gray(y_min - Width/2:y_max + Width/2, x_min - Width/2:x_max + Width/2);
    cropped_mask = Mask(y_min - Width/2:y_max + Width/2, x_min - Width/2:x_max + Width/2);
    
    [f1, d1] = vl_sift(IMG1_gray);
    [f2, d2] = vl_sift(IMG2_gray);
    
    matches = vl_ubcmatch(d1, d2);
    matched_pts1 = round(f1(1:2, matches(1, :))');
    matched_pts2 = round(f2(1:2, matches(2, :))');
    matched_pts1_in_bounds = [];
    matched_idxs = [];
    
    
    for i = 1:length(matched_pts1)
        coord = [matched_pts1(i, 2) matched_pts1(i, 1)];
        if cropped_mask(coord(1), coord(2)) == 1
            matched_pts1_in_bounds = [matched_pts1_in_bounds; matched_pts1(i, :)];
            matched_idxs = [matched_idxs; i];
        end
    end
    
    matched_pts2_in_bounds = matched_pts2(matched_idxs, :);
    
    figure(3)
    showMatchedFeatures(IMG1_gray, IMG2_gray, matched_pts1_in_bounds, matched_pts2_in_bounds, 'montage');
    
    adjustment = ones(size(matched_pts1_in_bounds));
    adjustment(:, 1) = adjustment(:, 1) .* (x_min - Width/2);
    adjustment(:, 2) = adjustment(:, 2) .* (y_min - Width/2);
    
    matched_pts1_in_bounds = matched_pts1_in_bounds + adjustment;
    matched_pts2_in_bounds = matched_pts2_in_bounds + adjustment;
    
    tform = estimateGeometricTransform(matched_pts1_in_bounds, matched_pts2_in_bounds, 'affine');
    
    WarpedFrame = imwarp(IMG1, tform, 'OutputView', imref2d(size(IMG1)));
    WarpedMask = imwarp(Mask, tform, 'OutputView', imref2d(size(IMG1)));
    WarpedMaskOutline = bwperim(WarpedMask, 4);
    WarpedLocalWindows = round(transformPointsForward(tform, Windows));
    WarpedLocalWindows = WarpedLocalWindows(~any(isnan(WarpedLocalWindows) | ...
        isinf(WarpedLocalWindows), 2), :);
    %[~, WarpedLocalWindows] = initLocalWindows(WarpedFrame, WarpedMask, length(Windows), Width, true); 
end
 