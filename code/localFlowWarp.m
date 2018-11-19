function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.

% WarpedPrevFrame = out
% CurrentFrame = ingImg
% LocalWindows = windows
% Mask = 
% Width = wSize
% Calculate tform
    opticFlow = opticalFlowFarneback;
    estimateFlow(opticFlow, rgb2gray(WarpedPrevFrame));
    flow = estimateFlow(opticFlow, rgb2gray(CurrentFrame));
    plot(flow);
    half_wwidth = floor(Width / 2);
    row_disp = half_wwidth;
    col_disp = half_wwidth;
    
%     CurrentFrame = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
%     CurrentFrame_gray = rgb2gray(CurrentFrame);
    Mask = padarray(Mask, [half_wwidth half_wwidth], 0, 'both');
    
    NewLocalWindows = zeros(size(LocalWindows));
    window_count = 1;
    for window_center = LocalWindows'
        % Row and columns displacement required for center of window due to
        % padding
        center = [(window_center(2) + row_disp) (window_center(1) + col_disp)];
        
        % Create window
%         window = CurrentFrame_gray(center(1) - half_wwidth:center(1) + half_wwidth, ...
%             center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict mask to window's size and location
        window_mask = Mask(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth, :);
        % Restrict flow to window's size and location
        window_Vx = flow.Vx(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        window_Vy = flow.Vy(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        
        in_bounds_Vx = window_Vx .* window_mask;
        in_bounds_Vy = window_Vy .* window_mask;
        Vx_avg = sum(in_bounds_Vx(:)) / (Width * Width);
        Vy_avg = sum(in_bounds_Vy(:)) / (Width * Width);
        
        disp(Vx_avg)
        disp(Vy_avg)
        
        if isnan(Vx_avg) || isnan(Vy_avg)
            Vx_avg = 0;
            Vy_avg = 0;
        end
        new_window_center_pos = window_center' + [Vx_avg Vy_avg];
        NewLocalWindows(window_count, :) = round(new_window_center_pos);
        window_count = window_count + 1;
    end
end

