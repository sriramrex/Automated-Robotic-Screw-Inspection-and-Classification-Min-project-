clear; clc; close all;

%% Read the images and assign them to the respective variables

folder1 = 'screw images top and front/';
file1 = 'second screw thread';
extension1 = '.bmp';
filepath1 = [folder1,file1,extension1];
OriginalImg = imread(filepath1);
%figure,imshow(OriginalImg)

folder2 = 'screw images top and front/';
file2 = 'Blank_Img';
extension2 = '.bmp';
filepath2 = [folder2,file2,extension2];
BlankImg = imread(filepath2);
%figure,imshow(BlankImg)

ScrewImg = BlankImg-OriginalImg;
figure,imshow(ScrewImg)

RotatedImg = imrotate(ScrewImg,5,'bilinear','loose');

figure,imshow(RotatedImg)