function [NewLocalWindows] = localFlowWarp(WarpedPrevFrame, CurrentFrame, LocalWindows, Mask, Width)
% LOCALFLOWWARP Calculate local window movement based on optical flow between frames.

% WarpedPrevFrame = out
% CurrentFrame = ingImg
% LocalWindows = windows
% Mask = 
% Width = wSize
% Calculate tform
opticFlow = opticalFlowFarneback();
flow = estimateFlow(opticFlow, rgb2gray(WarpedPrevFrame));
windows = LocalWindows;

for i=1:size(windows,2)
    window = windows{i};
    X = round(window(1) - (Width/2));
    Y = round(window(2) - (Width/2));
    XX = X + Width;
    YY = Y + Width;
    if(XX >= size(CurrentFrame, 2))
        XX = size(CurrentFrame, 2);
        X = XX - (Width);
    end
    
    if(YY >= size(CurrentFrame,1))
       YY = size(CurrentFrame,1);
       Y = YY - Width;
    end
    
    Vx = flow.Vx;
    Vy = flow.Vy;
    Vx = Vx(Y:YY, X:XX);
    Vx(window) = NaN;
    
    Vy = Vy(Y:YY, X:XX);
    Vy(window) = NaN;
    
    avg_Vx = (mean(mean(Vx,2,'omitnan'),1,'omitnan'));
    avg_Vy = (mean(mean(Vy,2,'omitnan'),1,'omitnan'));
    
    window = [(window(1) + avg_Vx) ...
        (window(2) + avg_Vy)];
    
    windows{i} = window;
end

t = Width/2;


for i=1:numel(windows)
    pos = windows{i};
    X = round(pos(1));
    Y = round(pos(2));
    if(X + (Width/2) >= size(CurrentFrame,2))
        X = size(CurrentFrame,1) - (Width/2);
    end
    
    if(Y + (Width/2) >= size(CurrentFrame,2))
        Y = size(CurrentFrame,1) - (Width/2);
    end
    
    new_img = CurrentFrame(Y-t:Y+t,X-t:X+t,:);
    
    windows{i} = new_img;
end

NewLocalWindows = windows;

end

