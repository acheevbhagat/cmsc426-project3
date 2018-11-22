function [mask, LocalWindows, ColorModels, ShapeConfidences] = ...
    updateModels(...
        NewLocalWindows, ...
        LocalWindows, ...
        CurrentFrame, ...
        WarpedMask, ...
        WarpedMaskOutline, ...
        WindowWidth, ...
        NumWindows, ...
        BoundaryWidth, ...
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
    
    [NewColorModels, remasked, NewMask, NewMaskOutline, NewWindows] = update_color_models(CurrentFrame, WarpedMask, ...
        WarpedMaskOutline, NewLocalWindows, BoundaryWidth, WindowWidth, NumWindows, ColorModels, ProbMaskThreshold);
    if remasked == true
        WarpedMask = NewMask;
        WarpedMaskOutline = NewMaskOutline;
        NewLocalWindows = NewWindows;
    end
    
    ShapeModels = update_shape_models(NewLocalWindows, NewColorModels, WindowWidth, ...
        SigmaMin, A, fcutoff, R);
    
    all_f_s = ShapeModels.Confidences;
    all_p_c = NewColorModels.ForegroundProbs;
    
    half_wwidth = floor(WindowWidth / 2);
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
    WarpedMask = padarray(WarpedMask, [half_wwidth half_wwidth], 0, 'both');
    WarpedMaskOutline = padarray(WarpedMaskOutline, [half_wwidth half_wwidth], 0, 'both');
    pix_dists_to_boundary = bwdist(WarpedMaskOutline);
    
    row_disp = half_wwidth;
    col_disp = half_wwidth;
    
    p_k_f_windows = cell(1, length(NewLocalWindows));
    window_count = 1;
    for window_center = NewLocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];

        % Restrict mask to window's size and location
        window_mask = WarpedMask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        
        f_s = all_f_s{window_count};
        p_c = all_p_c{window_count};
        
        %p_k_f = zeros(size(window_mask));
        p_k_f = f_s .* window_mask + (ones(size(f_s)) - f_s) .* p_c;
        
        p_k_f_windows{window_count} = p_k_f;
        window_count = window_count + 1;
    end
    
    numers = zeros([size(IMG, 1) size(IMG, 2)]);
    denoms = zeros([size(IMG, 1) size(IMG, 2)]);
    epsilon = 0.1;
    window_count = 1;
    for window_center = NewLocalWindows'
        center = [window_center(2) + row_disp, window_center(1) + col_disp];
        
        for row = 1:WindowWidth+1
            for col = 1:WindowWidth+1
                curr_p_k_f_x = p_k_f_windows{window_count}(row, col);
                dist = sqrt((col - center(2))^2 + (row - center(1))^2);
                
                % Get numerators and denominators for all pixels of global
                % foreground mask
                numers(center(1) - half_wwidth + row, center(2) - half_wwidth + col) = ...
                    numers(center(1) - half_wwidth + row, center(2) - half_wwidth + col) + ...
                    curr_p_k_f_x * (dist + epsilon)^-1;
                denoms(center(1) - half_wwidth + row, center(2) - half_wwidth + col) = ...
                    denoms(center(1) - half_wwidth + row, center(2) - half_wwidth + col) + ...
                    (dist + epsilon)^-1;
            end
        end
        window_count = window_count + 1;
    end
    
    p_f_matrix = (numers ./ denoms);
    imshow(p_f_matrix);
    p_f_matrix(isnan(p_f_matrix)) = 0;
    mask = p_f_matrix;
%     mask = p_f_matrix(row_disp + 1:size(p_f_matrix, 1) - row_disp, ...
%         col_disp + 1:size(p_f_matrix, 2) - col_disp);
    
    y_min = min(NewLocalWindows(:, 2));
    y_max = max(NewLocalWindows(:, 2));
    x_min = min(NewLocalWindows(:, 1));
    x_max = max(NewLocalWindows(:, 1));
    
    cropped_frame = IMG(y_min - half_wwidth + row_disp:y_max + half_wwidth + row_disp, ...
        x_min - half_wwidth + col_disp:x_max + half_wwidth + col_disp);
    mask = mask(y_min - half_wwidth + row_disp:y_max + half_wwidth + row_disp, ...
        x_min - half_wwidth + col_disp:x_max + half_wwidth + col_disp);
    mask = mask > ProbMaskThreshold;
    mask = lazysnapping(cropped_frame, superpixels(cropped_frame, 5000), mask, ~mask);
    final_mask = zeros([size(IMG, 1) size(IMG, 2)]);
    
    for r = 1:size(mask, 1)
        for c = 1:size(mask, 2)
            final_mask(r + y_min - half_wwidth + row_disp - 1, c + x_min + col_disp - half_wwidth - 1) = ...
                mask(r, c);
        end
    end
    
    mask = final_mask(row_disp + 1:size(final_mask, 1) - row_disp, ...
        col_disp + 1:size(final_mask, 2) - col_disp);
    LocalWindows = NewLocalWindows;
end

function [ColorModels, remasked, NewMask, NewMaskOutline, NewWindows] = update_color_models(CurrentFrame, WarpedMask, ...
    WarpedMaskOutline, NewLocalWindows, BoundaryWidth, WindowWidth, NumWindows, OldColorModels, ProbMaskThreshold)
    confidences = cell(1, size(NewLocalWindows, 1));
    distances = cell(1, size(NewLocalWindows, 1));
    foreground_probs = cell(1, size(NewLocalWindows, 1));
    window_masks = cell(1, size(NewLocalWindows, 1));
    window_gmms = cell(1, size(NewLocalWindows, 1));
    
    half_wwidth = floor(WindowWidth / 2);
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
    WarpedMask = padarray(WarpedMask, [half_wwidth half_wwidth], 0, 'both');
    WarpedMaskOutline = padarray(WarpedMaskOutline, [half_wwidth half_wwidth], 0, 'both');
    pix_dists_to_boundary = bwdist(WarpedMaskOutline);
    
    row_disp = half_wwidth;
    col_disp = half_wwidth;

    IMG_Lab = rgb2lab(IMG);
    
    % p_c matrix from ColorModels of previous frame
    prev_models_GMMs = OldColorModels.GMMs;
    prev_confidences = OldColorModels.Confidences;
    
    % Thresholds for sample data for updated GMMs
    high_thresh = 0.85;
    low_thresh = 0.70;
    
    % Loop through all windows along boundary
    window_count = 1;
    total_num_new_pixels = 0;
    total_num_prev_pixels = 0;
    for window_center = NewLocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];
        
        % Create window
        window = IMG_Lab(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict mask to window's size and location
        window_mask = WarpedMask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        
        % Retrieve background and foreground color data at least 5 pixels away from boundary
        pix_dists_to_boundary_window = pix_dists_to_boundary(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        
        window_Lab_vals = reshape(window, [size(window, 1)^2 3]);
        
        % get p_c matrix from applying previous GMMs to current frame's
        % window
        prev_models = prev_models_GMMs{window_count};
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
        options = statset('MaxIter', 1000);
        if size(F_Lab_vals, 1) > size(F_Lab_vals, 2) && size(B_Lab_vals, 1) > size(B_Lab_vals, 2)
            new_F_gmm = fitgmdist(F_Lab_vals, 3, 'RegularizationValue', 0.0005, 'Options', options);
            new_B_gmm = fitgmdist(B_Lab_vals, 3, 'RegularizationValue', 0.0005, 'Options', options);
        else
            new_F_gmm = prev_F_gmm;
            new_B_gmm = prev_B_gmm;
        end
        
        % Combine foreground and background probabilities
        new_model_p_c_matrix = get_fore_prob(new_F_gmm, new_B_gmm, window_Lab_vals, WindowWidth + 1);
        
        num_F_pixels_new_GMM = length(find(new_model_p_c_matrix > ProbMaskThreshold));
        num_F_pixels_prev_GMM = length(find(prev_model_p_c_matrix > ProbMaskThreshold));
        
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
            color_model_confidence = prev_confidences{window_count};
            p_c_matrix = prev_model_p_c_matrix;
            F_gmm = prev_F_gmm;
            B_gmm = prev_B_gmm;
        end
        
        total_num_new_pixels = total_num_new_pixels + num_F_pixels_new_GMM;
        total_num_prev_pixels = total_num_prev_pixels + num_F_pixels_prev_GMM;
        
        confidences{window_count} = color_model_confidence;
        distances{window_count} = pix_dists_to_boundary_window;
        foreground_probs{window_count} = p_c_matrix;
        window_masks{window_count} = window_mask;
        window_gmms{window_count} = {F_gmm B_gmm};
        window_count = window_count + 1;
    end
    remasked = false;
    if (1 - total_num_prev_pixels / total_num_new_pixels) > 0.18 || ...
            (1 - total_num_new_pixels / total_num_prev_pixels) > 0.18
        disp('Remask');
        beep
        remasked = true;
        new_mask = roipoly();
        [new_mask_outline, new_LocalWindows] = initLocalWindows(CurrentFrame, new_mask, ...
            NumWindows,WindowWidth,true);
        new_model = initColorModels(CurrentFrame, new_mask, new_mask_outline, ...
            new_LocalWindows, BoundaryWidth, WindowWidth);
        
        confidences = new_model.Confidences;
        distances = new_model.Distances;
        foreground_probs = new_model.ForegroundProbs;
        window_masks = new_model.SegmentationMasks;
        window_gmms = new_model.GMMs;
        NewMask = new_mask;
        NewMaskOutline = new_mask_outline;
        NewWindows = new_LocalWindows;
    end
    
    figure(4);
    imshow(foreground_probs{1});
    ColorModels = struct('Confidences', {confidences}, 'Distances', {distances}, ...
        'ForegroundProbs', {foreground_probs}, 'SegmentationMasks', {window_masks}, ...
        'GMMs', {window_gmms});
    if remasked == false
        NewMask = [];
        NewMaskOutline = [];
        NewWindows = [];
    end
end

function ShapeModels = update_shape_models(NewLocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R)
    ShapeModels = initShapeConfidences(NewLocalWindows, ColorModels, WindowWidth, SigmaMin, A, fcutoff, R);
end

% Get probabilities of pixels being in the foreground
function p_c = get_fore_prob(F_gmm, B_gmm, data, window_width)
    F_likelihood = pdf(F_gmm, data);
    B_likelihood = pdf(B_gmm, data);
    p_c = F_likelihood ./ (F_likelihood + B_likelihood);
    
    % reshape into dimensions of original window (pdf() produces nx1 vector)
    p_c = reshape(p_c, [window_width window_width]);
end