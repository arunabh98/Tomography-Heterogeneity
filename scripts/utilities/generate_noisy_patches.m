%% Occlusion Generator
%  ------------------------------------------------------
%  Generates images of a person with added occlusions
%  ------------------------------------------------------
directory = '../data/outlier_1_refs_012/';
if(exist(directory) == 0)
    mkdir(directory);
end


% Get the image.
image_size = 100;
image_name = 'refs_012.png';
P = imread(strcat('../data/', image_name));
P = imresize(P, [image_size image_size]);
P = im2double(rgb2gray(P));

% Make sure the image is symmetric by taking an mirror image.
flipped_P2 = flip(P, 2);
P((image_size/2) + 1:end, (image_size/2) + 1:end) = flipped_P2((image_size/2) + 1:end, (image_size/2) + 1:end);
flipped_P1 = flip(P, 1);
P(1:(image_size/2), :) = flipped_P1(1:(image_size/2), :);

% Check if the image is symmetric.
left_half = P(:, 1:(image_size/2));
right_half = flip(P(:, (image_size/2) + 1:end), 2);
disp('**** Check if the image is symmetric ****');
disp(norm(left_half - right_half, 'fro'));
disp('');

% Rotate the image.
P = imrotate(P, 90);
% Pad the image with a fixed boundary of 3 pixels.
img = padarray(P, [3, 3], 0.0);
[m,n] = size(img);

numSamples = 200;
minPatchSize = 15;
maxPatchSize = 25;
minNumPatches = 4;
maxNumPatches = 8;

for i = 1:numSamples
    newImg = img(:,:);
    numPatches = randi([minNumPatches, maxNumPatches]);
    for j=1:numPatches
        patchSize = randi([minPatchSize, maxPatchSize]);
        startY = randi([1,m-patchSize+1]);
        startX = randi([1,n-patchSize+1]);
        newImg(startY:startY+patchSize-1,startX:startX+patchSize-1) = zeros(patchSize);
    end
    fileName = sprintf('%s/%d.png',directory,i);
    imwrite(newImg,fileName);
end