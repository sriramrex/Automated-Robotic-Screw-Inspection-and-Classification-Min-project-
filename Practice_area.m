%% Detection of screw Pitch, Crest, Root, Depth parameters

%finding depth in the right side of the screw

len1 = size(centersIn,1);
len2 = size(centersOut,1);
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
p_min = min(centersIn(:,1));
p_max = max(centersIn(:,1));
p_mean = (p_max+p_min)/2;

% for right side of the screw thread
u1 = 1;
for k = 1:len2  
     if centersOut(k,1)>p_mean 
        xOut(u1,1) = centersOut(k,1);
        yOut(u1,1) = centersOut(k,2);
        u1=u1+1;
     end
end
u2 = 1;
for k = 1:len1
     if centersIn(k,1)>p_mean 
        xIn(u2,1) = centersIn(k,1);
        yIn(u2,1) = centersIn(k,2);
        u2=u2+1;
     end
end

for k = 1:size(yOut,1)
    for i1 = 1:size(yIn,1)
    s_yp(i1,k)= yOut(k,1)-yIn(i1,1);
    end
end

for m = 1:size(s_yp,2)
e1 = 1;
d_yp = [];
for k = 1:size(s_yp,1)
    if s_yp(k,m)>0
        d_yp(e1,1) = (s_yp(k,m));
         e1 = e1+1;   
    end
end
short_distp(m,1) =min(d_yp); 
end
nearest_ptp = yOut-short_distp;

for i  = 1:size(nearest_ptp,1)
nearest_ptp(i,2) = find(centersIn(:,2) == nearest_ptp(i,1));
sn = nearest_ptp(i,2); 
nearest_ptp(i,3) = centersIn(sn,1);
nearest_ptp(i,4) = radiiIn(sn,1);
end

for i  = 1:size(nearest_ptp,1)
nearest_ptp(i,5) = yOut(i,1);
end

for i  = 1:size(nearest_ptp,1)
nearest_ptp(i,6) = find(centersOut(:,2) == nearest_ptp(i,5));
 sn = nearest_ptp(i,6); 
 nearest_ptp(i,7) = centersOut(sn,1);
 nearest_ptp(i,8) = radiiOut(sn,1);
end

%% Thread angles for the screw 

dist_rightp(:,1) = nearest_ptp(:,3);
dist_rightp(:,2) = nearest_ptp(:,1)+nearest_ptp(:,4);
dist_rightp(:,3) = nearest_ptp(:,7);
dist_rightp(:,4) = nearest_ptp(:,5)-nearest_ptp(:,8);

dist_rightp(:,5) = abs((dist_rightp(:,4)-dist_rightp(:,2))./(dist_rightp(:,3)-dist_rightp(:,1)));

% for k1 = 1:size(dist_rightp(:,5))
%     if dist_rightp(k1,5)<0
%         dist_rightp(k1,:)=[];
%     end
% % %dist_rightp(:,5) = filloutliers(dist_rightp(:,5),"nearest");
%  end

Dist = dist_rightp(:,5)+dist_right(:,5);

detect = isoutlier(Dist);

Angle = rmoutliers(Dist,"OutlierLocations",detect);



