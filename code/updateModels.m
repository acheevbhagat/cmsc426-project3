function [mask, LocalWindows, ColorModels, ShapeConfidences] = ...
    updateModels(...
        NewLocalWindows, ...
        LocalWindows, ...
        CurrentFrame, ...
        WarpedMask, ...
        WarpedMaskOutline, ...
        WindowWidth, ...
        BoundaryWidth, ...
        ColorModels, ...
        ShapeConfidences, ...
        ProbMaskThreshold, ...
        fcutoff, ...
        SigmaMin, ...
        R, ...
        A, ...
        origColorModel ...
    )
% UPDATEMODELS: update shape and color models, and apply the result to generate a new mask.
% Feel free to redefine this as several different functions if you prefer.
    NewColorModels = update_color_models(CurrentFrame, WarpedMask, WarpedMaskOutline, ...
        NewLocalWindows, BoundaryWidth, WindowWidth, ColorModels, origColorModel);
    ShapeModels = update_shape_models(NewLocalWindows, NewColorModels, WindowWidth, ...
        SigmaMin, A, fcutoff, R);
    
    all_f_s = ShapeModels.Confidences;
    all_p_c = NewColorModels.ForegroundProbs;
    
    half_wwidth = floor(WindowWidth / 2);
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
    Mask = padarray(WarpedMask, [half_wwidth half_wwidth], 0, 'both');
    MaskOutline = padarray(WarpedMaskOutline, [half_wwidth half_wwidth], 0, 'both');
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
    mask = p_f_matrix(half_wwidth + 1:size(p_f_matrix, 1) - half_wwidth, ...
        half_wwidth + 1:size(p_f_matrix, 2) - half_wwidth);
    mask = mask > ProbMaskThreshold;
    LocalWindows = NewLocalWindows;
end

function ColorModels = update_color_models(CurrentFrame, WarpedMask, WarpedMaskOutline, ...
    NewLocalWindows, BoundaryWidth, WindowWidth, OldColorModels, origColorModel)
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
    prev_models_GMMs = OldColorModels.GMMs;
    prev_confidences = OldColorModels.Confidences;
    
    % Thresholds for sample data for updated GMMs
    high_thresh = 0.75;
    low_thresh = 0.25;
    
    % Flag to keep track of whether ColorModels were totally recalculated.
    % If it does, then break the window loop and ask user to create new
    % mask, and create colormodels based on that mask.
    recalced = false;

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
        origGMM = origColorModel.GMMs{window_count};
        prev_F_gmm = origGMM{1};
        prev_B_gmm = origGMM{2};
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
        if size(F_Lab_vals, 1) > size(F_Lab_vals, 2) && size(B_Lab_vals, 1) > size(B_Lab_vals, 2)
            new_F_gmm = fitgmdist(F_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
            new_B_gmm = fitgmdist(B_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
        elseif size(F_Lab_vals, 1) < size(F_Lab_vals, 2)
            new_F_gmm = prev_F_gmm;
            new_B_gmm = fitgmdist(B_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
        elseif size(B_Lab_vals, 1) < size(B_Lab_vals, 2)
            new_F_gmm = fitgmdist(F_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
            new_B_gmm = prev_B_gmm;
        else
            new_F_gmm = prev_F_gmm;
            new_B_gmm = prev_B_gmm;
        end
        
        
        
        % Combine foreground and background probabilities
        new_model_p_c_matrix = get_fore_prob(new_F_gmm, new_B_gmm, window_Lab_vals, WindowWidth + 1);
        
        num_F_pixels_new_GMM = length(find(new_model_p_c_matrix==1));
        num_F_pixels_prev_GMM = length(find(prev_model_p_c_matrix==1));
        
        % If number of foreground pixels doesn't increase, use new GMM,
        % else, use old GMM
        if num_F_pixels_new_GMM >= num_F_pixels_prev_GMM
            % If the difference is greater than some percentage, detect an
            % error and ask the user to create a new mask for use by the
            % color model
            if (1 - (num_F_pixels_prev_GMM / num_F_pixels_new_GMM)) > 0.75
                recalced = true;
                break;
            else
                confidence_numer = 0;
                confidence_denom = 0;
                for row = 1:size(window, 1)
                    for col = 1:size(window, 2)
                        weight = exp(-pix_dists_to_boundary_window(row, col)^2 / (WindowWidth / 2)^2);
                        p_c = new_model_p_c_matrix(row, col);
                        confidence_numer = confidence_numer + abs(window_mask(row, col) - p_c) * weight;
                        confidence_denom = confidence_denom + weight;
                    end
                end
                color_model_confidence = 1 - (confidence_numer / confidence_denom);
                p_c_matrix = new_model_p_c_matrix;
                F_gmm = new_F_gmm;
                B_gmm = new_B_gmm;
            end
        else
            color_model_confidence = prev_confidences{window_count};
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
    
    if recalced
        disp("Margin of error is too large, please create new mask");
        mask = roipoly(CurrentFrame);
        mask_outline = bwperim(mask, 4);
        ColorModels = initColorModels(CurrentFrame, mask, ...
            mask_outline, NewLocalWindows, BoundaryWidth, WindowWidth);
        foreground_probs = {};
        i = 1;
        for window_center = NewLocalWindows'
            center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];
            window = IMG_Lab(center(1) - half_wwidth:center(1) + half_wwidth, ...
                center(2) - half_wwidth:center(2) + half_wwidth, :);
            window_Lab_vals = reshape(window, [size(window, 1)^2 3]);
            GMM = ColorsModels.GMMs{i};
            fGMM = GMM{1};
            bGMM = GMM{2};
            foreground_probs{i} = get_fore_prob(fGMM, bGMM, ...
                window_Lab_vals, WindowWidth + 1);
            i = i + 1;
        end
    else

        ColorModels = struct('Confidences', {confidences}, 'Distances', {distances}, ...
            'ForegroundProbs', {foreground_probs}, 'SegmentationMasks', {window_masks}, ...
            'GMMs', {window_gmms});
    end
    imshow(foreground_probs{1});
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