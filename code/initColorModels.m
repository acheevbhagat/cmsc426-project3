function ColorModels = initializeColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
    % INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
    %
    % Must define a field ColorModels.Confidences: a cell array of the color confidence map for each local window.
    
    confidences = cell(1, size(LocalWindows, 1));
    distances = cell(1, size(LocalWindows, 1));
    foreground_probs = cell(1, size(LocalWindows, 1));
    window_masks = cell(1, size(LocalWindows, 1));
    window_gmms = cell(1, size(LocalWindows, 1));
    
    half_wwidth = floor(WindowWidth / 2);
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(IMG, [half_wwidth half_wwidth], 0, 'both');
    IMG_Lab = rgb2lab(IMG);
    Mask = padarray(Mask, [half_wwidth half_wwidth], 0, 'both');
    MaskOutline = padarray(MaskOutline, [half_wwidth half_wwidth], 0, 'both');
    pix_dists_to_boundary = bwdist(MaskOutline);
    
    row_disp = half_wwidth;
    col_disp = half_wwidth;

    % Loop through all windows along boundary
    window_count = 1;
    for window_center = LocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];
        
        % Create window
        window = IMG_Lab(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict mask to window's size and location
        window_mask = Mask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        
        B_Lab_vals = [];
        F_Lab_vals = [];
        
        % Retrieve background and foreground color data at least 5 pixels away from boundary
        pix_dists_to_boundary_window = pix_dists_to_boundary(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        %imshow(mat2gray(pix_dists_to_boundary_window, [0 20]));
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                if window_mask(row, col) == 0 && ...
                        pix_dists_to_boundary_window(row, col) > BoundaryWidth
                    B_Lab_vals = [B_Lab_vals; window(row, col, 1) window(row, col, 2) window(row, col, 3)];
                end
                if window_mask(row, col) == 1 && ...
                        pix_dists_to_boundary_window(row, col) > BoundaryWidth
                    F_Lab_vals = [F_Lab_vals; window(row, col, 1) window(row, col, 2) window(row, col, 3)];
                end
            end
        end
        
        % Fit GMM models
        options = statset('MaxIter', 1000);
        F_gmm = fitgmdist(F_Lab_vals, 3, 'RegularizationValue', 0.0005, 'Options', options);
        B_gmm = fitgmdist(B_Lab_vals, 3, 'RegularizationValue', 0.0005, 'Options', options);
        
        % Part 3.3 of project notes
        window_Lab_vals = reshape(window, [size(window, 1)^2 3]);
        p_c_matrix = get_fore_prob(F_gmm, B_gmm, window_Lab_vals, WindowWidth + 1);
        %imshow(p_c_matrix)
        %imshow(reshape(pdf(B_gmm, window_Lab_vals), [size(window, 1) size(window, 1)]));
        pdf(F_gmm, window_Lab_vals);
        
        confidence_numer = 0;
        confidence_denom = 0;
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                weight = exp(-pix_dists_to_boundary_window(row, col)^2 / (WindowWidth / 2)^2);
                p_c = p_c_matrix(row, col);
                confidence_numer = confidence_numer + abs(window_mask(row, col) - p_c) * weight;
                confidence_denom = confidence_denom + weight;
            end
        end
        
        color_model_confidence = 1 - (confidence_numer / confidence_denom);
        confidences{window_count} = color_model_confidence;
        distances{window_count} = pix_dists_to_boundary_window;
        foreground_probs{window_count} = p_c_matrix;
        window_masks{window_count} = window_mask;
        window_gmms{window_count} = {F_gmm B_gmm};
        window_count = window_count + 1;
    end
%     figure(5)
%     imshow(foreground_probs{1});
    ColorModels = struct('Confidences', {confidences}, 'Distances', {distances}, ...
        'ForegroundProbs', {foreground_probs}, 'SegmentationMasks', {window_masks}, ...
        'GMMs', {window_gmms});
end

% Get probabilities of pixels being in the foreground
function p_c = get_fore_prob(F_gmm, B_gmm, data, window_width)
    F_likelihood = pdf(F_gmm, data);
    B_likelihood = pdf(B_gmm, data);
    p_c = F_likelihood ./ (F_likelihood + B_likelihood);
    
    % reshape into dimensions of original window (pdf() produces nx1 vector)
    p_c = reshape(p_c, [window_width window_width]);
end
