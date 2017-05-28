%% Making a distribution
X = 2; % binning interval

cd(sum_dir);
distribution = zeros(250,250);
distances_angle = struct([]);
binrange = 1:100;
bincenter=binrange(1:end) - 0.5;
N= zeros(360/X, 100);
curve = struct([]);
gof = struct([]);

radius = zeros(1, 360/X);
sigma1 = zeros(1, 360/X);
pvalue = zeros(1, 360/X);
intensity = zeros(1,360/X);

for l=1:counter
    
    for m=1:length(b_cells_rot2{l})
        distribution(round(b_cells_rot2{l}(m,2))+125,round(b_cells_rot2{l}(m,1))+125) = ...
            distribution(round(b_cells_rot2{l}(m,2))+125,round(b_cells_rot2{l}(m,1))+125)+1;
    end
end

image2 = figure;
imshow(distribution, [min(distribution(:)) max(distribution(:))]);
print(image2, '-dtiff', '-r300', 'distribution.tif');
cd(currdir);


%% Radius variation measuring angles and distances of each border point
angle = zeros(1,1);
distance = zeros(1,1);
p=0;
for l=1:counter
    for m=1:length(b_cells_rot2{l})
        p=p+1;
        distance(p) = sqrt(b_cells_rot2{l}(m,2)*b_cells_rot2{l}(m,2) +...
            b_cells_rot2{l}(m,1)*b_cells_rot2{l}(m,1));
        
        a = b_cells_rot2{l}(m,2)/distance(p);
        b = b_cells_rot2{l}(m,1)/distance(p);
        angle(p)= atan2d(a,b);
        if angle(p)<0
            angle(p)= angle(p)+360;
        end
    end
end


%% obtaining radial distribution by angles
for j=X:X:360
    distances_angle{j}= zeros(1,1);
    for q=1:length(distance)
        if angle(q)<=j && angle(q)>j-X
            distances_angle{j/X}=[distances_angle{j/X};distance(q)];
        end  
    end
    distances_angle{j/X} =  distances_angle{j/X}( distances_angle{j/X}~=0);
end

%% fitting radius
for i=X:X:360
    [N(i/X,:), bins] = histc(distances_angle{i/X},binrange);
    [curve{i/X},gof{i/X}] = fit(bincenter',N(i/X,:)','gauss1');
    radius(i/X) = curve{i/X}.b1;
    sigma1(i/X) = sqrt(curve{i/X}.c1*curve{i/X}.c1/2);
    pvalue(i/X) = gof{i/X}.rsquare;
end

image3=figure;
cd(sum_dir);
plot(X:X:360, radius);
print(image3, '-dtiff', '-r300', 'circle_radius.tif');

image4=figure;
plot(X:X:360, sigma1);
print(image4, '-dtiff', '-r300', 'circle_width.tif');



for i=X:X:360
    intensity(i/X) = N(i/X,round(radius(i/X))-1) + N(i/X,round(radius(i/X))) + N(i/X,round(radius(i/X))+1);
end

image5=figure;
plot(X:X:360,intensity);
print(image5, '-dtiff', '-r300', 'circle_intensity.tif');

%close all