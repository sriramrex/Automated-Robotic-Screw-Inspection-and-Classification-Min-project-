folder = 'new head position images/';
file = 'seventh_head';
extension = '.bmp';
filepath = [folder,file,extension];
colorImg = imread(filepath);
grayImg = rgb2gray(colorImg);
T = adaptthresh(grayImg);
figure,imshow(grayImg);

Img = imbinarize(grayImg, 0.22);
figure,imshow(Img);