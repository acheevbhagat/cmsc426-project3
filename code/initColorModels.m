function ColorModels = initializeColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
    % INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
    %
    % Must define a field ColorModels.Confidences: a cell array of the color confidence map for each local window.
    
    confidences = {};
    if mod(WindowWidth, 2) == 0
         half_wwidth = WindowWidth / 2;
    else
         half_wwidth = (WindowWidth - 1) / 2;
    end
    
    % Pad everything to avoid index out of bounds
    IMG = padarray(IMG, [half_wwidth half_wwidth], 0, 'both');
    Mask = padarray(Mask, [half_wwidth half_wwidth], 0, 'both');
    MaskOutline = padarray(MaskOutline, [half_wwidth half_wwidth], 0, 'both');
    pix_dists_to_boundary = bwdist(MaskOutline);
    
    row_disp = half_wwidth;
    col_disp = half_wwidth;

    IMG_Lab = rgb2lab(IMG);
    
    % Loop through all windows along boundary
    for window_center = LocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [window_center(2) + row_disp window_center(1) + col_disp];
        
        % Create window
        window = IMG_Lab(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict mask to window's size and location
        window_mask = Mask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        
        % Mask for foreground and background
        % (Must adjust these to leave 5 pixel gap between samples and
        % boundaries)
        %B_mask = zeros(size(window));
        %F_mask = zeros(size(window));
        
        B_Lab_vals = [];
        F_Lab_vals = [];
        
        % Retrieve background and foreground color data at least 5 pixels away from boundary
        pix_dists_to_boundary_mask = pix_dists_to_boundary(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                if window_mask(row, col) == 0 && ...
                        pix_dists_to_boundary_mask(row, col) > 5
                    %B_mask(row, col) = 1;
                    B_Lab_vals = [B_Lab_vals; window(row, col, 1) window(row, col, 2) window(row, col, 3)];
                end
                if window_mask(row, col) == 1 && ...
                        pix_dists_to_boundary_mask(row, col) > 5
                    %F_mask(row, col) = 1;
                    F_Lab_vals = [F_Lab_vals; window(row, col, 1) window(row, col, 2) window(row, col, 3)];
                end
            end
        end
        F_Lab_vals
        
%         size(window)
%         window_L = double(window(:, :, 1));
%         window_a = double(window(:, :, 2));
%         window_b = double(window(:, :, 3));
        
        % Retrieve foreground and background samples from image using masks
%         F_Lab_vals = [(window_L .* F_mask) (window_a .* F_mask) (window_b .* F_mask)];
%         B_Lab_vals = [(window_L .* B_mask) (window_a .* B_mask) (window_b .* B_mask)];
        
        % Fit GMM models
        options = statset('MaxIter', 500);
        F_gmm = fitgmdist(F_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
        B_gmm = fitgmdist(B_Lab_vals, 3, 'RegularizationValue', 0.001, 'Options', options);
        
        
        % Part 3.3 of project notes
        window_Lab_vals = reshape(window, [size(window,1)^2 3]);
        F_probs = get_fore_prob(F_gmm, B_gmm, window_Lab_vals);
        confidence_numer = 0;
        confidence_denom = 0;
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                weight = exp(-pix_dists_to_boundary_mask(row, col)^2 / (WindowWidth / 2)^2);
                pix = [window(row, col, 1); window(row, col, 2); window(row, col, 3)];
                p_c = F_probs(row, col);
                confidence_numer = confidence_numer + abs(window_mask(row, col) - p_c) * weight;
                confidence_denom = confidence_denom + weight;
            end
        end
        
        color_model_confidence = 1 - (confidence_numer / confidence_denom)
        confidences = {confidences color_model_confidence};
    end
    ColorModels = struct('Confidences', confidences);
end

function likelihood = get_likelihood(pix, mu, covar)
    likelihood = exp(-.5*(pix-mu)'*(covar\(pix-mu))) / sqrt((2*pi)^3 * det(covar));
end

% Get probability of pixel being in foreground
function p_c = get_fore_prob(F_gmm, B_gmm, data)
    F_likelihood = pdf(F_gmm, data);
    B_likelihood = pdf(B_gmm, data);
    p_c = F_likelihood / (F_likelihood + B_likelihood);
end

