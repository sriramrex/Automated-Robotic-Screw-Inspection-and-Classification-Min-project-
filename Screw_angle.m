clear; clc; close all;

%% Read the images and assign them to the respective variables

folder1 = 'screw images top and front/';
file1 = 'first screw image';
extension1 = '.bmp';
filepath1 = [folder1,file1,extension1];
OriginalImg = imread(filepath1);
figure,imshow(OriginalImg)

folder2 = 'screw images top and front/';
file2 = 'blank_Img';
extension2 = '.bmp';
filepath2 = [folder2,file2,extension2];
BlankImg = imread(filepath2);
figure,imshow(BlankImg)

ScrewImg = BlankImg-OriginalImg;
figure,imshow(ScrewImg)

%% binarizing the image:

[pixelCount,grayLevels]= imhist(ScrewImg);
Threshold = otsuthresh(pixelCount);
BrightImg = imbinarize(ScrewImg,Threshold);
figure, imshow(BrightImg);

%% Connected Components 

conn=8; 
labeledImg = bwlabel(BrightImg, conn); 
colLabels = label2rgb (labeledImg, 'hsv', 'k', 'shuffle');
figure,imshow(colLabels);

%% Regional properties:

Out_BrightImg = BrightImg;

% Remove portions of the image that touch an outside edge.
Out_BrightImg = imclearborder(Out_BrightImg);

% Fill holes in regions.
Out_BrightImg = imfill(Out_BrightImg, 'holes');

% Filter image based on image properties.
Out_BrightImg = bwpropfilt(Out_BrightImg,'Area',[5000 + eps(5000), Inf]);

% Get properties.
properties = regionprops(Out_BrightImg, {'Area', 'Eccentricity', 'EquivDiameter', 'EulerNumber', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter'});

% Show the output image
figure, imshow(Out_BrightImg);
%% Morphological Operations

r = 13;
SE = strel('diamond',r);

OpenedRegion = imopen(Out_BrightImg,SE);
figure,imshow(OpenedRegion);

%% Remove small areas using Regional properties

SeperateImg = OpenedRegion;

% Remove portions of the image that touch an outside edge.
SeperateImg = imclearborder(SeperateImg);

% Fill holes in regions.
SeperateImg = imfill(SeperateImg, 'holes');

% Filter image based on image properties.
SeperateImg = bwpropfilt(SeperateImg,'Area',[2000 + eps(2000), Inf]);

% Get properties.
properties = regionprops(SeperateImg, {'Area', 'Eccentricity', 'EquivDiameter', 'EulerNumber', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter'});
 
figure, imshow(SeperateImg);

% r = 20;
% SE = strel('square',r);
% 
% SeperateImg = imopen(SeperateImg,SE);
% figure,imshow(SeperateImg);

[Boundary,Location] = bwboundaries(SeperateImg);

%% Removing the top head of the screw image

ThreadImg = SeperateImg;
ThreadImg = imclearborder(ThreadImg);
ThreadImg = imfill(ThreadImg, 'holes');
ThreadImg = bwpropfilt(ThreadImg,'Area',[112295 + eps(112295), Inf]);
properties = regionprops(ThreadImg, {'Area', 'Eccentricity', 'EquivDiameter', 'EulerNumber', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter'});

%% Apply Median Filter to filer the noise in the image (Rightside)
FThreadImg = medfilt2(ThreadImg,[20,20]);
I = imcrop(FThreadImg,[1407.51 900 103.98 1110]);
figure;imshow(I)


%% Edge detection

Threshold = 0.4;
%BW = edge(SeperateImg,"canny", Threshold,"both","nothinning"); 
EdgeScrewImg = edge(I,'canny',Threshold);
[r, c] = find(EdgeScrewImg); % Edge point coordinates
figure, imshow(EdgeScrewImg);

% BW4 = bwperim(EdgeScrewImg,18);
% Threadedge = imdilate(BW4, strel('diamond',11));
% figure, imshow(Threadedge);

% %% Hough lines
% 
% 
% [H, T, R] = hough(EdgeScrewImg);
% P = houghpeaks(H, 2, 'threshold', 0.8, 'NHoodSize', [1 1]);             
% lines = houghlines(EdgeScrewImg, T, R, P, 'FillGap', 20, 'MinLength', 1);
%          
% figure;
% imshow(EdgeScrewImg);
% hold on;                            
% for ii = 1:length(lines)                                                  
% 	xy = [lines(ii).point1; lines(ii).point2];                               
% 	plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'green');
% end          
% 
% %%
% Point1 = transpose({ lines(:).point1});
% 
% Point1 =[Point1 transpose({ lines(:).point2})];

%% Hugh circle detection 

[centersIn,radiiIn] = imfindcircles(ThreadImg,[2 20], "Method","TwoStage","EdgeThreshold",0.4,"Sensitivity",0.8);
[centersOut,radiiOut] = imfindcircles(ThreadImg,[2 20], "Method","TwoStage","EdgeThreshold",0.4,"Sensitivity",0.8,"ObjectPolarity","dark");
figure,imshow(ThreadImg), hold on
viscircles(centersIn, radiiIn);
viscircles (centersOut, radiiOut);
len1 = size(centersIn,1);
len2 = size(centersOut,1);
for k = 1:len1
    plot(centersIn(k,1), centersIn(k,2), 'r+', 'MarkerSize', 1, 'LineWidth', 2);
end

for k = 1:len2
    plot(centersOut(k,1), centersOut(k,2), 'w+', 'MarkerSize', 1, 'LineWidth', 2);
end


%% Screw thread angle Right

len3 = size(centersIn,1);
len4 = size(centersOut,1);
t1 = 1;
t2 = 1;

for k = 1 : len1         % Loop through all blobs.
	% Place the blob label number at the centroid of the blob.
	text(centersIn(k,1), centersIn(k,2), num2str(k), 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'cap','Color','k');
end
for k = 1 : len2        % Loop through all blobs.
	% Place the blob label number at the centroid of the blob.
	text(centersOut(k,1), centersOut(k,2), num2str(k), 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'cap','Color','w');
end

q_min = min(centersIn(:,1));
q_max = max(centersIn(:,1));
q_mean = (q_max+q_min)/2;

u1 = 1;
for k = 1:len4  
     if centersOut(k,1)>q_mean 
        xaOut(u1,1) = centersOut(k,1);
        yaOut(u1,1) = centersOut(k,2);
        u1=u1+1;
     end
end
u2 = 1;
for k = 1:len3
     if centersIn(k,1)>q_mean 
        xaIn(u2,1) = centersIn(k,1);
        yaIn(u2,1) = centersIn(k,2);
        u2=u2+1;
     end
end


for k = 1:size(yaOut,1)
    for i1 = 1:size(yaIn,1)
    s_a(i1,k)= yaIn(i1,1)-yaOut(k,1);
    end
end

for m = 1:size(s_a,2)
e1 = 1;
d_a = [];
for k = 1:size(s_a,1)
    if s_a(k,m)>0
        d_a(e1,1) = (s_a(k,m));
         e1 = e1+1;   
    end
end
short_distright(m,1) =min(d_a); 
end

nearest_right = yaOut+short_distright;

for i  = 1:size(nearest_right,1)
nearest_right(i,2) = find(centersIn(:,2) == nearest_right(i,1));
sr = nearest_right(i,2); 
nearest_right(i,3) = centersIn(sr,1);
nearest_right(i,4) = radiiIn(sr,1);
end

for i  = 1:size(nearest_right,1)
nearest_right(i,5) = yaOut(i,1);
end

for i  = 1:size(nearest_right,1)
nearest_right(i,6) = find(centersOut(:,2) == nearest_right(i,5));
 sr = nearest_right(i,6); 
 nearest_right(i,7) = centersOut(sr,1);
 nearest_right(i,8) = radiiOut(sr,1);
end

dist_right(:,1) = nearest_right(:,3);
dist_right(:,2) = nearest_right(:,1)-nearest_right(:,4);
dist_right(:,3) = nearest_right(:,7);
dist_right(:,4) = nearest_right(:,5)+nearest_right(:,8);

dist_right(:,5) = (dist_right(:,2)-dist_right(:,4))./(dist_right(:,1)-dist_right(:,3));


