function ColorModels = initializeColorModels(IMG, Mask, MaskOutline, LocalWindows, BoundaryWidth, WindowWidth)
    % INITIALIZAECOLORMODELS Initialize color models.  ColorModels is a struct you should define yourself.
    %
    % Must define a field ColorModels.Confidences: a cell array of the color confidence map for each local window.
    IMG_Lab = rgb2lab(IMG);
    ColorModels.Confidences = {};
    for window_center = LocalWindows'
        if mod(WindowWidth, 2) == 0
             half_wwidth = WindowWidth / 2;
        else
             half_wwidth = (WindowWidth - 1) / 2;
        end
        
        IMG_Lab = padarray(IMG_Lab, [half_wwidth half_wwidth], 0, 'both');
        Mask = padarray(Mask, [half_wwidth half_wwidth], 0, 'both');
        MaskOutline = padarray(MaskOutline, [half_wwidth half_wwidth], 0, 'both');
        size(IMG_Lab)
        row_disp = half_wwidth;
        col_disp = half_wwidth;
        
        hold on;
        imshow(IMG_Lab(:, :, 1));
        hold off;
        
        window = IMG_Lab(window_center(1) + row_disp - half_wwidth:window_center(1) + row_disp + half_wwidth, ...
            window_center(2) + col_disp - half_wwidth:window_center(2) + col_disp + half_wwidth);
        size(window)
        % Mask for foreground and background
        % (Must adjust these to leave 5 pixel gap between samples and
        % boundaries)
        B_mask = zeros(size(window));
        F_mask = zeros(size(window));
        
        % Retrieve background data at least 5 pixels away from boundary
        pix_dists_to_boundary = bwdist(MaskOutline);
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                if Mask(window_center(1) + row_disp + row, window_center(2) + col_disp + col) == 0 ...
                    && pix_dists_to_boundary(row + row_disp, col + col_disp) > 5
                    B_mask(row, col) = 1;
                end
            end
        end
        
        % Retrieve foreground data at least 5 pixels away from boundary
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                if Mask(row, col) == 1 && pix_dists_to_boundary(row, col) > 5
                    F_mask(row, col) = 1;
                end
            end
        end
        
        F_mask
        B_mask
        
        % Retrieve foreground and background samples from image using masks
        img_L = double(IMG_Lab(:, :, 1));
        img_a = double(IMG_Lab(:, :, 2));
        img_b = double(IMG_Lab(:, :, 3));
        F_Lab_vals = [img_L(F_mask) img_a(F_mask) img_b(F_mask)];
        B_Lab_vals = [img_L(B_mask) img_a(B_mask) img_b(B_mask)];
        
        F_gmm = fitgmdist(F_Lab_vals, 3);
        B_gmm = fitgmdist(B_Lab_vals, 3);
        
%         F_probs = zeros(size(window));
%         for row = 1:size(window, 1)
%             for col = 1:size(window, 2)
%                 F_likelihood = get_likelihood(window(row, col), F_gmm.mu, F_gmm.sigma);
%                 B_likelihood = get_likelihood(window(row, col), B_gmm.mu, B_gmm.sigma);
%                 p_c = F_likelihood / (F_likelihood + B_likelihood);
%                 F_probs(row, col) = p_c;
%             end
%         end
        
        confidence_numer = 0;
        confidence_denom = 0;
        for row = 1:size(window, 1)
            for col = 1:size(window, 2)
                weight = exp(-pix_dists_to_boundary(row, col)^2 / (WindowWidth / 2)^2);
                p_c = get_color_prob(window(row, col), F_gmm, B_gmm);
                confidence_numer = confidence_numer + abs(Mask(row, col) - p_c) * weight;
                confidence_denom = confidence_denom + weight;
            end
        end
        
        color_model_confidence = 1 - (confidence_numer / confidence_denom);
    end
end

function likelihood = get_likelihood(pix, mu, covar)
    likelihood = exp(-.5*(pix-mu)'*(covar\(pix-mu))) / sqrt((2*pi)^3 * det(covar));
end

function p_c = get_color_prob(pix, F_gmm, B_gmm)
    F_likelihood = get_likelihood(window(row, col), F_gmm.mu, F_gmm.sigma);
    B_likelihood = get_likelihood(window(row, col), B_gmm.mu, B_gmm.sigma);
    p_c = F_likelihood / (F_likelihood + B_likelihood);
end

