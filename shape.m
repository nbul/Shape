%% Clear all and initial parameters
clc
clear variables
close all

%% Variables
se90V = strel('line', 8, 90);
se0V = strel('line', 8, 0);
%dilation of borders
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);
b_cells_final = struct([]);
corners_final = struct([]);
k=0;
kk=0;

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8_dir =[filedir, '/tifs_8bit'];
%Folder to save summarised information
mkdir(filedir,'/Summary_shape');
sum_dir = [filedir, '/Summary_shape'];

cd(tif8_dir);
files_tif = dir('*.tif');

%%Processing data
for g=1:numel(files_tif)
    %% Open images and modify
    bd_dir = [tif8_dir,'/', num2str(g)];
    cd(bd_dir);
    %Open image with vertices
    V = imread('vertices.png');
    V = im2bw(V,1/255);
    % Vdil - dilated image of all vertices
    Vdil = imdilate(V, [se90V se0V]);
    Vdil = imclearborder(Vdil);
    % Individual vertices as objects
    cc_Vdil = bwconncomp(Vdil);
    N_borders = labelmatrix(cc_Vdil);
    s_borders = regionprops(cc_Vdil, Vdil, 'Centroid');
    I=imread('tracked_bd.png');
    I2=im2bw(I,1/255);
    I3 = imdilate(I2, [se90I se0I]);
    % I_cells - inverted image of all cells that are completely in frame;
    % s_cells - individual cells as objects
    I_cells=imcomplement(I3);
    I_cells=imclearborder(I_cells);
    cc_cells=bwconncomp(I_cells);
    L_cells=labelmatrix(cc_cells);
    s_cells=regionprops(cc_cells, I_cells, 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Orientation','PixelList','Centroid', 'Area');
    
    Pixels = struct([]);
    b_cells = struct([]);
    b_cells_flip = struct([]);
    b_cells_rot = struct([]);
    b_cells_rot2 = struct([]);
    corners = struct([]);
    corners_good = struct([]);
    corners_good_rot = struct([]);
    corners_good_rot2 = struct([]);
    
    for i=1:numel(s_cells)
        if s_cells(i).Area >200
            clear Center_mat_pixels Center_mat_border
            Pixels{i} = s_cells(i).PixelList;
            
            minima = min(s_cells(i).PixelList);
            Center = round(s_cells(i).Centroid);
            b_cells{i} = bwtraceboundary(I_cells,[Pixels{i}(1,2),Pixels{i}(1,1)],'E',4);
            b_cells_flip{i} = [ b_cells{i}(:,2), b_cells{i}(:,1)];
            
            %% Getting cornrers coordinates and selectin >30%
            for c=1:length(b_cells{i})
                corners{i}(c) = N_borders(b_cells{i}(c,1),b_cells{i}(c,2));
            end
            
            
            corners{i} = corners{i}(corners{i}~=0);
            corners{i} = unique(corners{i});
            corners_good{i} = zeros(1,2);
            cc=0;
            if length(corners{i})>2
                corners{i} = [corners{i},corners{i}(1), corners{i}(2)];
                CosTheta = zeros(1,length(corners{i})-2);
                ThetaInDegrees = zeros(1,length(corners{i})-2);
                for c=2:length(corners{i})-1
                    vector1 = s_borders(corners{i}(c+1)).Centroid - s_borders(corners{i}(c)).Centroid;
                    vector2 = s_borders(corners{i}(c-1)).Centroid - s_borders(corners{i}(c)).Centroid;
                    CosTheta(c-1) = dot( vector1, vector2)/(norm( vector1)*norm( vector2));
                    ThetaInDegrees(c-1) = acosd(CosTheta(c-1));
                    
                    if ThetaInDegrees(c-1)>0 && ThetaInDegrees(c-1)<120
                        cc = cc+1;
                        corners_good{i}(cc,:) = s_borders(corners{i}(c)).Centroid;
                    end
                end
            end
            
            
            %% Get boundaries and move objects to center at 0-0 coordinates
            R = [cosd(s_cells(i).Orientation) -sind(s_cells(i).Orientation); sind(s_cells(i).Orientation) cosd(s_cells(i).Orientation)];
            Center_mat_pixels = repmat([Center(1); Center(2)], 1, length(Pixels{i}))';
            Center_mat_border = repmat([Center(1); Center(2)], 1, length(b_cells_flip{i}))';
            
            b_cells_rot{i} = b_cells_flip{i} - Center_mat_border;
            b_cells_rot{i} = (R*b_cells_rot{i}')';
            b_cells_rot2{i}(:,1) = b_cells_rot{i}(:,1)*100/s_cells(i).MajorAxisLength;
            b_cells_rot2{i}(:,2) = b_cells_rot{i}(:,2)*100/s_cells(i).MinorAxisLength;
            k=k+1;
            b_cells_final{k} = b_cells_rot2{i};
            
            if length(corners_good{i})>2
                Center_mat_border2 = repmat([Center(1); Center(2)], 1, length(corners_good{i}))';
                corners_good_rot{i} = corners_good{i} - Center_mat_border2;
                corners_good_rot{i} = (R*corners_good_rot{i}')';
                corners_good_rot2{i}(:,1) = corners_good_rot{i}(:,1)*100/s_cells(i).MajorAxisLength;
                corners_good_rot2{i}(:,2) = corners_good_rot{i}(:,2)*100/s_cells(i).MinorAxisLength;
                
                %% Stretching
                if length(corners_good_rot2{i})>2
                    kk=kk+1;
                    corners_final{kk} = corners_good_rot2{i};
                end
                
            end
            
        end
    end
end



%% Making a distribution
distribution = zeros (250,250);
image1 = figure;
for l=1:k
    
    for m=1:length(b_cells_final{l})
        distribution(round(b_cells_final{l}(m,2))+125,round(b_cells_final{l}(m,1))+125) = ...
            distribution(round(b_cells_final{l}(m,2))+125,round(b_cells_final{l}(m,1))+125)+1;
    end
    
    plot(b_cells_final{l}(:,1),b_cells_final{l}(:,2));
    hold on
end
hold off

cd(sum_dir);
print(image1, '-dtiff', '-r300', 'all_boundaries.tif');

image2 = figure;
imshow(distribution, [min(distribution(:)) max(distribution(:))]);
print(image2, '-dtiff', '-r300', 'distribution.tif');
cd(currdir);


%% Radius variation
angle = zeros(1,1);
distance = zeros(1,1);
p=0;

distances_angle = struct([]);
% measuring angles and distances of each border point
for l=1:k
    for m=1:length(b_cells_final{l})
        p=p+1;
        distance(p) = sqrt(b_cells_final{l}(m,2)*b_cells_final{l}(m,2) +...
            b_cells_final{l}(m,1)*b_cells_final{l}(m,1));
        
        a = b_cells_final{l}(m,2)/distance(p);
        b = b_cells_final{l}(m,1)/distance(p);
        angle(p)= atan2d(a,b);
        if angle(p)<0
            angle(p)= angle(p)+360;
        end
    end
end

% obtaining radial distribution by abgles
for j=1:360
    distances_angle{j}= zeros(1,1);
    for q=1:length(distance)
        if angle(q)<=j && angle(q)>j-1
            distances_angle{j}=[distances_angle{j};distance(q)];
        end
    end
end

% fitting radius
binrange = 1:100;
bincenter=binrange(1:(end-1)) + 0.5;
N= zeros(360, 100);
curve = struct([]);
gof = struct([]);

radius = zeros(1, 360);
sigma1 = zeros(1, 360);
pvalue = zeros(1, 360);
good_ring = zeros(1, 360);


for i=1:360
    [N(i,:), bins] = histc(distances_angle{i},binrange);
    binsnew = 1:length(N(i,:));
    [curve{i},gof{i}] = fit(binsnew',N(i,:)','gauss1');
    radius(i) = curve{i}.b1;
    sigma1(i) = sqrt(curve{i}.c1*curve{i}.c1/2);
    pvalue(i) = gof{i}.rsquare;
end
image3=figure;
cd(sum_dir);
plot(1:360, radius);
print(image3, '-dtiff', '-r300', 'circle_radius.tif');

image4=figure;
plot(1:360, sigma1);
print(image4, '-dtiff', '-r300', 'circle_width.tif');


intensity = zeros(1,360);
for i=1:360
    intensity(i) = N(i,round(radius(i))-1) + N(i,round(radius(i))) + N(i,round(radius(i))+1);
end

image5=figure;
plot(1:360,intensity);
print(image5, '-dtiff', '-r300', 'circle_intensity.tif');


%% Corner distribution
image6 = figure;
for ll=1:kk
scatter(corners_final{ll}(:,1),corners_final{ll}(:,2));
hold on
end
hold off

print(image6, '-dtiff', '-r300', 'all_corners.tif');
corner_angle = zeros(1,1);
distance_corner = zeros(1,1);
pp=0;
for l=1:kk
    for m=1:length(corners_final{l})
        pp=pp+1;
        distance_corner(pp) = sqrt(corners_final{l}(m,2)*corners_final{l}(m,2) +...
            corners_final{l}(m,1)*corners_final{l}(m,1));
        
        a = corners_final{l}(m,2)/distance_corner(pp);
        b = corners_final{l}(m,1)/distance_corner(pp);
        corner_angle(pp)= atan2d(a,b);
        if corner_angle(pp)<0
            corner_angle(pp)= corner_angle(pp)+360;
        end
    end
end

number_corners = zeros(1,360);
for j=1:360
    for q=1:length(corner_angle)
        if corner_angle(q)<=j && corner_angle(q)>j-1
            number_corners(j)=number_corners(j)+1;
        end
    end
end
image7=figure;
plot(1:360,number_corners);
print(image7, '-dtiff', '-r300', 'number_corners.tif');

cd(currdir)

close all;
clear variables;
clc

