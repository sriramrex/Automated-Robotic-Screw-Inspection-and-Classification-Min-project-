clear; clc; close all;

%% Read the image from the folder and find the threshold using Otsu thresholding 

folder = 'new head position images/';
file = 'third_head';
extension = '.bmp';
filepath = [folder,file,extension];
colorImg = imread(filepath);
grayImg = rgb2gray(colorImg);
[counts,binLocations] = imhist(grayImg);
T = otsuthresh(counts);
figure,imshow(grayImg);


%% Image padding to convert the image into 1456 x 1088

pad_Img = padarray(colorImg,[636 1004],0,'both');
figure,imshow(pad_Img);
imwrite(pad_Img,"Binary_Img1.bmp")

%% image croping

crop_Img= imcrop(im_test_rect,[526.5 358.5 420 416]);
figure,imshow(crop_Img);

%%

%% Detecting the diameter of the circular screws using Hough Circle transform
invImg = ~Img;
[centers,radii] = imfindcircles(invImg,[100 400], 'Sensitivity', 0.9,'EdgeThreshold', 0.03,'Method', 'TwoStage', 'ObjectPolarity', 'dark');
figure,imshow(invImg), hold on
viscircles(centers, radii);
len1 = size(centers,1);
count_Circles = 0;
for k = 1:len1
    plot(centers(k,1), centers(k,2), 'r+', 'MarkerSize', 1, 'LineWidth', 2);
    count_Circles = count_Circles+1;
end
Dia_Screw = 2*max(radii);
%% Detecting the width of the hexagon screws using Hough line tansform

% if count_Circles == 0
% else
%     sprintf("This is a Circular Screw head")
% end

% doubleImg = im2double(colorImg);
bwImg = imbinarize(grayImg);
[edgeImg, threshOut] = edge(bwImg, 'Canny', 0.28,1.8);
figure,imshow(edgeImg);

[Hmatrix, Theta, Rho] = hough(edgeImg, 'Theta', -90.00:7.00:83.00, 'RhoResolution', 2.32);
Peaks = houghpeaks(Hmatrix, 2, 'threshold', 0.50, 'NHoodSize', [9 13]);             
lines = houghlines(edgeImg, Theta, Rho, Peaks, 'FillGap', 20, 'MinLength', 40);
        
%%VISUALIZATION:
figure;
imshow(colorImg);
hold on;                            
for ii = 1:length(lines)                                                  
	xy = [lines(ii).point1; lines(ii).point2];                               
	plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
end 


