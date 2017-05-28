%% Corner distribution
X2 = 3; % binning interval for corners

image6 = figure;
corner_rot2=corner_rot2(~cellfun(@isempty, corner_rot2));
for ll=1:length(corner_rot2)
    if corner_rot2{ll}
        scatter(corner_rot2{ll}(:,1),corner_rot2{ll}(:,2));
        hold on
    end
end
hold off

print(image6, '-dtiff', '-r300', 'all_corners.tif');

%% Distance and angle for each corner
corner_angle = zeros(1,1);
distance_corner = zeros(1,1);
pp = 0;

for l=1:length(corner_rot2)
    for m=1:length(corner_rot2{l})
        pp=pp+1;
        distance_corner(pp) = sqrt(corner_rot2{l}(m,2)*corner_rot2{l}(m,2) +...
            corner_rot2{l}(m,1)*corner_rot2{l}(m,1));
        
        a = corner_rot2{l}(m,2)/distance_corner(pp);
        b = corner_rot2{l}(m,1)/distance_corner(pp);
        corner_angle(pp)= atan2d(a,b);
        if corner_angle(pp)<0
            corner_angle(pp)= corner_angle(pp)+360;
        end
    end
end

%% Number of corners by angle
number_corners = zeros(1,360/X2);
for j=X2:X2:360
    for q=1:length(corner_angle)
        if corner_angle(q)<=j && corner_angle(q)>j-X2
            number_corners(j/X2)=number_corners(j/X2)+1;
        end
    end
end
image7=figure;
plot(X2:X2:360,number_corners);
print(image7, '-dtiff', '-r300', 'number_corners.tif');

%% Distance of corners by angleN= zeros(360, 100);
curve_corners = struct([]);
gof_corners = struct([]);
distance_temp = struct([]);

radius_corners = zeros(1, 360/X2);
sigma1_corners = zeros(1, 360/X2);
pvalue_corners = zeros(1, 360/X2);
N_corners= zeros(X2, 100);

for i=X2:X2:360
    distance_temp{i/X2}= zeros(1,1);
    for l=1:length(distance_corner)
        if corner_angle(l)<=i && corner_angle(l)>i-X2
            distance_temp{i/X2} = [distance_temp{i/X2}; distance_corner(l)];
        end
    end
    distance_temp{i/X2} =  distance_temp{i/X2}(distance_temp{i/X2}~=0);
    [N_corners(i/X2,:), bins] = histc(distance_temp{i/X2},binrange);
    [curve_corners{i/X2},gof_corners{i/X2}] = fit(bincenter',N_corners(i/X2,:)','gauss1');
    radius_corners(i/X2) = curve_corners{i/X2}.b1;
    sigma1_corners(i/X2) = sqrt(curve_corners{i/X2}.c1*curve_corners{i/X2}.c1/2);
    pvalue_corners(i/X2) = gof_corners{i/X2}.rsquare;
end
image8=figure;
cd(sum_dir);
plot(X2:X2:360, radius_corners);
print(image8, '-dtiff', '-r300', 'corners_radius.tif');


