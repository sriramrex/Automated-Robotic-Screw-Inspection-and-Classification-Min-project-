len1 = size(centersIn,1);
len2 = size(centersOut,1);
t1 = 1;
t2 = 1;
p_min = min(centersIn(:,1));
p_max = max(centersIn(:,1));
p_mean = (p_max+p_min)/2;
for k = 1:len1
     if centersIn(k,1)<p_mean 
        x1(t1) = centersIn(k,1);
        y1(t1) = centersIn(k,2);
         t1= t1+1;
     end
end
for k = 1:len1
     if centersIn(k,1)>p_mean
        x2(t2) = centersIn(k,1);
        y2(t2) = centersIn(k,2);
         t2= t2+1;
     end
end
%%
figure,imshow(ThreadImg), hold on
%set(gca,'View',[-180 90])
p1 = polyfit(x1,y1,1); 
f1 = polyval(p1,x1); 
plot(x1,y1,'o',x1,f1,'-','LineWidth',1,'Color','b')  

p2 = polyfit(x2,y2,1); 
f2 = polyval(p2,x2); 
plot(x2,y2,'o',x2,f2,'-','LineWidth',1,'Color','b')  



%% Practice for finding depth
t1 = 1;
for k = 1:len1
     if dist_y(k,1)>0 
        s_y(t1) = dist_y(k,1);
         t1= t1+1;
     end
end
%%
A = [ 3.5 2 1.6 1.456 0 1.9 2.6 ; 3.8 2.6 3.9 0 6 1.564 0  ];
minvals = mink(A, 2, 2);  %first 2 miminums along dimension 2.
sb = minvals(:, 2);
%%

I1 = Out_BrightImg-OpenedRegion;
figure,imshow(I1);
[pixelCount,grayLevels]= imhist(I1);
Threshold = otsuthresh(pixelCount);
BI1 = imbinarize(I1,Threshold);
figure, imshow(BI1);

%%
BW_out = BI1;

BW_out = imclearborder(BW_out);


BW_out = imfill(BW_out, 'holes');


BW_out = bwpropfilt(BW_out,'Area',[-Inf, 250 - eps(250)]);

properties = regionprops(BW_out, {'Area', 'Eccentricity', 'EquivDiameter', 'EulerNumber', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter'});

figure,imshow(BW_out);

%%

figure,imshow(SeperateImg);
CloseImg = SeperateImg+BW_out;
figure,imshow(CloseImg);

[pixelCount,grayLevels]= imhist(CloseImg);
Threshold = otsuthresh(pixelCount);
CloseLogicImg = imbinarize(CloseImg,Threshold);

ClearImg = CloseLogicImg;
ClearImg = imclearborder(ClearImg);
ClearImg= imfill(ClearImg, 'holes');
ClearImg = bwpropfilt(ClearImg,'Area',[112295 + eps(112295), Inf]);
properties = regionprops(ClearImg, {'Area', 'Eccentricity', 'EquivDiameter', 'EulerNumber', 'MajorAxisLength', 'MinorAxisLength', 'Orientation', 'Perimeter'});
figure,imshow(ClearImg);


%%

figure,  imshow(SeperateImg);
  roi = drawline('Color','r');

  %%

  figure, imshow(SeperateImg)
 roi = drawline('Color','r');

addlistener(roi,'MovingROI',@allevents);
addlistener(roi,'ROIMoved',@allevents);

function allevents(src,evt)
    evname = evt.EventName;
    switch(evname)
        case{'MovingROI'}
            disp(['ROI moving previous position: ' mat2str(evt.PreviousPosition)]);
            disp(['ROI moving current position: ' mat2str(evt.CurrentPosition)]);
        case{'ROIMoved'}
            disp(['ROI moved previous position: ' mat2str(evt.PreviousPosition)]);
            disp(['ROI moved current position: ' mat2str(evt.CurrentPosition)]);
    end
end

%%

for k = 1:size(yOut,1)
    for i1 = 1:size(yIn,1)
    ans1(i1,k)= yOut(k,1)-yIn(i1,1);
    end
end




