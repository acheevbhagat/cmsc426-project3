function [mask, LocalWindows, ColorModels, ShapeConfidences] = ...
    updateModels(...
        NewLocalWindows, ...
        LocalWindows, ...
        CurrentFrame, ...
        WarpedMask, ...
        WarpedMaskOutline, ...
        WindowWidth, ...
        ColorModels, ...
        ShapeConfidences, ...
        ProbMaskThreshold, ...
        fcutoff, ...
        SigmaMin, ...
        R, ...
        A ...
    )
% UPDATEMODELS: update shape and color models, and apply the result to generate a new mask.
% Feel free to redefine this as several different functions if you prefer.
    
    NewColorModels = update_color_models(CurrentFrame, WarpedMask, WarpedMaskOutline, ...
        NewLocalWindows, BoundaryWidth, WindowWidth, ColorModels);
    ShapeModels = update_shape_models(NewLocalWindows, NewColorModels, WindowWidth, ...
        SigmaMin, A, fcutoff, R);
    
    all_f_s = ShapeModels.Confidences;
    all_p_c = NewColorModels.ForegroundProbs;
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
    Mask = padarray(WarpedMask, [half_wwidth half_wwidth], 0, 'both');
    MaskOutline = padarray(WarpedMaskOutline, [half_wwidth half_wwidth], 0, 'both');
    pix_dists_to_boundary = bwdist(WarpedMaskOutline);
    
    row_disp = half_wwidth;
    col_disp = half_wwidth;
    
    window_count = 1;
    for window = NewLocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];

        % Create window
        window = IMG_Lab(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict mask to window's size and location
        window_mask = WarpedMask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        
        f_s = all_f_s{window_count};
        p_c = all_p_c{window_count};
        
        p_k_f = zeros(size(window));
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                f_s_x = f_s(row, col);
                p_c_x = p_c(row, col);
                
                p_k_f(row, col) = f_s_x * window_mask(row, col) + (1 - f_s_x) * p_c_x
            end
        end
    end
end

function ColorModels = update_color_models(CurrentFrame, WarpedMask, WarpedMaskOutline, ...
    NewLocalWindows, BoundaryWidth, WindowWidth, OldColorModels)
    confidences = cell(1, size(NewLocalWindows, 1));
    distances = cell(1, size(NewLocalWindows, 1));
    foreground_probs = cell(1, size(NewLocalWindows, 1));
    window_masks = cell(1, size(NewLocalWindows, 1));
    window_gmms = cell(1, size(NewLocalWindows, 1));
    
    half_wwidth = floor(WindowWidth / 2);
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
    Mask = padarray(WarpedMask, [half_wwidth half_wwidth], 0, 'both');
    MaskOutline = padarray(WarpedMaskOutline, [half_wwidth half_wwidth], 0, 'both');
    pix_dists_to_boundary = bwdist(WarpedMaskOutline);
    
    row_disp = half_wwidth;
    col_disp = half_wwidth;

    IMG_Lab = rgb2lab(IMG);
    
    % p_c matrix from ColorModels of previous frame
    p_c_prev_models = OldColorModels.GMMs;
    
    % Thresholds for sample data for updated GMMs
    high_thresh = 0.75;
    low_thresh = 0.25;
    
    % Loop through all windows along boundary
    window_count = 1;
    for window_center = NewLocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];
        
        % Create window
        window = IMG_Lab(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict mask to window's size and location
        window_mask = Mask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        
        % Retrieve background and foreground color data at least 5 pixels away from boundary
        pix_dists_to_boundary_window = pix_dists_to_boundary(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        
        window_Lab_vals = reshape(window, [size(window, 1)^2 3]);
        
        % get p_c matrix from applying previous GMMs to current frame's
        % window
        prev_models = p_c_prev_models{window_count};
        prev_F_gmm = prev_models{1};
        prev_B_gmm = prev_models{2};
        prev_model_p_c_matrix = get_fore_prob(prev_F_gmm, prev_B_gmm, window_Lab_vals, WindowWidth + 1);
        
        B_Lab_vals = [];
        F_Lab_vals = [];
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                if window_mask(row, col) == 0 && ...
                        pix_dists_to_boundary_window(row, col) > BoundaryWidth && ...
                        prev_model_p_c_matrix(row, col) < low_thresh
                    B_Lab_vals = [B_Lab_vals; window(row, col, 1) window(row, col, 2) window(row, col, 3)];
                end
                if window_mask(row, col) == 1 && ...
                        pix_dists_to_boundary_window(row, col) > BoundaryWidth && ...
                        prev_model_p_c_matrix(row, col) > high_thresh
                    F_Lab_vals = [F_Lab_vals; window(row, col, 1) window(row, col, 2) window(row, col, 3)];
                end
            end
        end
        
        % Fit GMM models
        options = statset('MaxIter', 700);
        new_F_gmm = fitgmdist(F_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
        new_B_gmm = fitgmdist(B_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
        
        % Combine foreground and background probabilities
        new_model_p_c_matrix = get_fore_prob(new_F_gmm, new_B_gmm, window_Lab_vals, WindowWidth + 1);
        
        num_F_pixels_new_GMM = length(find(new_model_p_c_matrix==1));
        num_F_pixels_prev_GMM = length(find(prev_model_p_c_matrix==1));
        
        % If number of foreground pixels doesn't increase, use new GMM,
        % else, use old GMM
        if num_F_pixels_new_GMM <= num_F_pixels_prev_GMM
            confidence_numer = 0;
            confidence_denom = 0;
            for row = 1:size(window, 1)
                for col = 1:size(window, 2)
                    weight = exp(-pix_dists_to_boundary_window(row, col)^2 / (WindowWidth / 2)^2);
                    p_c = prev_model_p_c_matrix(row, col);
                    confidence_numer = confidence_numer + abs(window_mask(row, col) - p_c) * weight;
                    confidence_denom = confidence_denom + weight;
                end
            end
            
            color_model_confidence = 1 - (confidence_numer / confidence_denom);
            p_c_matrix = new_model_p_c_matrix;
            F_gmm = new_F_gmm;
            B_gmm = new_B_gmm;
        else
            color_model_confidence = p_c_prev_models.Confidences{window_count};
            p_c_matrix = prev_model_p_c_matrix;
            F_gmm = prev_F_gmm;
            B_gmm = prev_B_gmm;
        end
        
        confidences{window_count} = color_model_confidence;
        distances{window_count} = pix_dists_to_boundary_window;
        foreground_probs{window_count} = p_c_matrix;
        window_masks{window_count} = window_mask;
        window_gmms{window_count} = {F_gmm B_gmm};
        window_count = window_count + 1;
    end
    imshow(foreground_probs{1});
    ColorModels = struct('Confidences', {confidences}, 'Distances', {distances}, ...
        'ForegroundProbs', {foreground_probs}, 'SegmentationMasks', {window_masks}, ...
        'GMMs', window_gmms);
end

function ShapeModels = update_shape_models(NewLocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R)
    ShapeModels = initShapeConfidences(NewLocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R);
end