clear; clc; close all;
%% read the images from the Ros publishers when it was taking images

% Run only once rosinit("host_name")
% To shutdown rosshutdown

%% For subscibing to a particular topic

img_sub = rossubscriber("/image_raw");
pause(2)
img = receive(img_sub,6);
image = readImage(img);
figure,imshow(image);
