% MyRotobrush.m  - UMD CMSC426, Fall 2018
% This is the main script of your rotobrush project.
% We've included an outline of what you should be doing, and some helful visualizations.
% However many of the most important functions are left for you to implement.
% Feel free to modify this code as you see fit.

% Some parameters you need to tune:
WindowWidth = 74;  
ProbMaskThreshold = 0.55; 
NumWindows= 35; 
BoundaryWidth = 5;

% Load images:
fpath = '../input';
new_fpath = '../bs';
files = dir(fullfile(fpath, '*.jpg'));
imageNames = zeros(length(files),1);
images = cell(length(files),1);

for i=1:length(files)
    imageNames(i) = str2double(strtok(files(i).name,'.jpg'));
end

imageNames = sort(imageNames);
imageNames = num2str(imageNames);
imageNames = strcat(imageNames, '.jpg');

for i=1:length(files)
    images{i} = im2double(imread(fullfile(fpath, strip(imageNames(i,:)))));
end
% NOTE: to save time during development, you should save/load your mask rather than use ROIPoly every time.
%mask = roipoly(images{1});
%save('mask4.mat', 'mask');
mask = load('mask4.mat');
mask = mask.mask;
imshow(imoverlay(images{1}, boundarymask(mask,8), 'red'));
set(gca,'position',[0 0 1 1],'units','normalized')
F = getframe(gcf);
[I,~] = frame2im(F);
imwrite(I, fullfile(new_fpath, strip(imageNames(1,:))));
outputVideo = VideoWriter(fullfile(new_fpath,'video.mp4'),'MPEG-4');
open(outputVideo);
writeVideo(outputVideo,I);

% Sample local windows and initialize shape+color models:
[mask_outline, LocalWindows] = initLocalWindows(images{1},mask,NumWindows,WindowWidth,true);

ColorModels = ...
    initColorModels(images{1},mask,mask_outline,LocalWindows,BoundaryWidth,WindowWidth);

origColorModel = ColorModels;

% You should set these parameters yourself:
fcutoff = 0.1;
SigmaMin = 20;
SigmaMax = WindowWidth + 1;
R = 2;
A = (SigmaMax - SigmaMin) / (1 - fcutoff)^R;
ShapeConfidences = ...
    initShapeConfidences(LocalWindows,ColorModels,...
    WindowWidth, SigmaMin, A, fcutoff, R);

% Show initial local windows and output of the color model:
imshow(images{1})
hold on
showLocalWindows(LocalWindows,WindowWidth,'r.');
hold off
set(gca,'position',[0 0 1 1],'units','normalized')
F = getframe(gcf);
[I,~] = frame2im(F);

%showColorConfidences(images{1},mask_outline,ColorModels.Confidences,LocalWindows,WindowWidth);

%%% MAIN LOOP %%%
% Process each frame in the video.
for prev=1:(length(files)-1)
    curr = prev+1;
    fprintf('Current frame: %i\n', curr)
    
    %%% Global affine transform between previous and current frames:
    [warpedFrame, warpedMask, warpedMaskOutline, warpedLocalWindows] = ...
        calculateGlobalAffine(images{prev}, images{curr}, mask, LocalWindows, WindowWidth);
    
    %%% Calculate and apply local warping based on optical flow:
    NewLocalWindows = ...
        localFlowWarp(warpedFrame, images{curr}, warpedLocalWindows,warpedMask,WindowWidth);
    
%    NewLocalWindows = warpedLocalWindows;
 
    % Show windows before and after optical flow-based warp:
    imshow(images{curr});
    hold on
    showLocalWindows(warpedLocalWindows,WindowWidth,'r.');
    showLocalWindows(NewLocalWindows,WindowWidth,'b.');
    hold off
    
    %%% UPDATE SHAPE AND COLOR MODELS:
    % This is where most things happen.
    % Feel free to redefine this as several different functions if you prefer.
    [ ...
        mask, ...
        LocalWindows, ...
        ColorModels, ...
        ShapeConfidences, ...
    ] = ...
    updateModels(...
        NewLocalWindows, ...
        LocalWindows, ...
        images{curr}, ...
        warpedMask, ...
        warpedMaskOutline, ...
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
    );

    mask_outline = bwperim(mask,4);

    showColorConfidences(images{1},mask_outline,ColorModels.Confidences,LocalWindows,WindowWidth);
    
    
    
    % Write video frame:
    imshow(imoverlay(images{curr}, boundarymask(mask,8), 'red'));
    set(gca,'position',[0 0 1 1],'units','normalized')
    F = getframe(gcf);
    [I,~] = frame2im(F);
    imwrite(I, fullfile(new_fpath, strip(imageNames(curr,:))));
    writeVideo(outputVideo,I);

    imshow(images{curr})
    hold on
    %showColorConfidences(images{1},mask_outline,ColorModels.Confidences,LocalWindows,WindowWidth);
    %showLocalWindows(LocalWindows,WindowWidth,'r.');
    hold off
    set(gca,'position',[0 0 1 1],'units','normalized')
    F = getframe(gcf);
    [I,~] = frame2im(F);
    
end

close(outputVideo);
