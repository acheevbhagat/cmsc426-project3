function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.

% WarpedPrevFrame = out
% CurrentFrame = ingImg
% LocalWindows = windows
% Mask = 
% Width = wSize
% Calculate tform
    opticFlow = opticalFlowFarneback('NeighborhoodSize', 8);
    estimateFlow(opticFlow, rgb2gray(WarpedPrevFrame));
    flow = estimateFlow(opticFlow, rgb2gray(CurrentFrame));
    half_wwidth = floor(Width / 2);
    row_disp = Width;
    col_disp = Width;
    
%     CurrentFrame = padarray(CurrentFrame, [half_wwidth half_wwidth], 0, 'both');
%     CurrentFrame_gray = rgb2gray(CurrentFrame);
    Mask = padarray(Mask, [Width Width], 0, 'both');
    flow_Vx = padarray(flow.Vx, [Width Width], 0, 'both');
    flow_Vy = padarray(flow.Vy, [Width Width], 0, 'both');
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
            center(2) - half_wwidth:center(2) + half_wwidth);
        % Restrict flow to window's size and location
        window_Vx = flow_Vx(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        window_Vy = flow_Vy(center(1) - half_wwidth:center(1) + half_wwidth, ...
            center(2) - half_wwidth:center(2) + half_wwidth);
        
        in_bounds_Vx = window_Vx .* window_mask;
        in_bounds_Vy = window_Vy .* window_mask;
        Vx_avg = sum(in_bounds_Vx(:)) / nnz(window_mask == 1);
        Vy_avg = sum(in_bounds_Vy(:)) / nnz(window_mask == 1);
        
        if isnan(Vx_avg) || isnan(Vy_avg) || isinf(Vx_avg) || isinf(Vy_avg)
            Vx_avg = 0;
            Vy_avg = 0;
        end
        
%         imshow(imfuse(WarpedPrevFrame, CurrentFrame));
%         hold on
%         plot(flow,'DecimationFactor',[5 5],'ScaleFactor',5);
%         hold off;
        
        new_window_center_pos = window_center' + [Vx_avg Vy_avg];
        NewLocalWindows(window_count, :) = round(new_window_center_pos');
        window_count = window_count + 1;
    end
end

