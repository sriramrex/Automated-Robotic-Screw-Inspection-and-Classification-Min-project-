clear; clc; close all;
%% Prepare the image files
filepath = cell(1,18);
for i = 1:18
folder = 'Calibration Checker board/';
filename = sprintf('img_%d',i);
extension = '.bmp'; 
filepath{i} = [folder,filename,extension];

end
I = imread(filepath{1});
imshow(I);
%% Estimate Camera Parameters

[imagePoints, boardSize]= detectCheckerboardPoints(filepath);
squareSize = 6;
worldPoints = generateCheckerboardPoints(boardSize,squareSize);

imageSize = [size(I,1),size(I,2)];
cameraParams = estimateCameraParameters(imagePoints,worldPoints,ImageSize = imageSize);

figure; showReprojectionErrors(cameraParams);
title("Reprojection Errors");
%% Read the image of Objects to be measured

folder1 = 'new head position images/';
filename1 = 'Binary_Img1';
extension1 = '.bmp';
img = [folder1,filename1,extension1];
ImOrig = imread(img);
%ImOrig = imresize(ImOrig,[ 1456 1088]);
figure;imshow(ImOrig)

 %% Undistort the image
magnification = 100;
[im, newOrigin] = undistortImage(ImOrig,cameraParams);
figure; imshow(im);
title("Undistorted Image");

%% Segment the Screw head
grayImg = rgb2gray(ImOrig);
T = adaptthresh(grayImg);
%figure,imshow(grayImg);
Img = imbinarize(grayImg, T);

Img = ~Img;
figure,imshow(Img)
 se = strel('disk',5);
 afterClosing = imclose(Img,se);
 %figure;imshow(afterClosing); 
SeperateHead = medfilt2(afterClosing,[1,1]);
figure;imshow(SeperateHead)

%% Compute Extrinsics

% Detect the checkerboard.
[imagePoints, boardSize] = detectCheckerboardPoints(im);

% Adjust the imagePoints so that they are expressed in the coordinate system
% used in the original image, before it was undistorted.  This adjustment
% makes it compatible with the cameraParameters object computed for the original image.
imagePoints = imagePoints + newOrigin; % adds newOrigin to every row of imagePoints

% Extract camera intrinsics.
camIntrinsics = cameraParams.Intrinsics;

% Compute extrinsic parameters of the camera.
camExtrinsics = estimateExtrinsics(imagePoints, worldPoints, camIntrinsics);





