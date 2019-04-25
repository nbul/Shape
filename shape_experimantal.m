clc;
clear variables;
close all;

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

%Folders with images
tif8_dir =[filedir, '/borders']; 
%Folder to save information about cells
mkdir(filedir,'/summary_shape');
sum_dir = [filedir, '/summary_shape'];
Area = zeros(1,1);
Perimeter = zeros(1,1);
MeanBD = zeros(1,1);
MaxBD = zeros(1,1);
Major = zeros(1,1);
Minor = zeros(1,1);
Circularity = zeros(1,1);
Solidity = zeros(1,1);
Ecc = zeros(1,1);
AR = zeros(1,1);
cd(tif8_dir);
files_tif = dir('*.tif');
counter = 0;
se90V = strel('line', 2, 90);
se0V = strel('line', 2, 0);
ThetaInDegrees = struct([]);

for g=1:numel(files_tif)
        
    bd_dir = [tif8_dir,'/', num2str(g)];
    cd(bd_dir);    
    I=imread('handCorrection.tif');
    I2=imbinarize(rgb2gray(I),0);
    I2(:,1) = 0;
    I2(:,end) = 0;
    I2(1,:) = 0;
    I2(end,:) = 0;
    [im_x, im_y] = size(I2);
    
    %% Cell shape from regionprops
    [B,L,N,A] = bwboundaries(I2,'holes');
    
    im_cells_data=regionprops(L,'Centroid', 'Area', 'Eccentricity','ConvexArea',...
        'MajorAxisLength', 'MinorAxisLength','Perimeter', 'Solidity','PixelList');
    for n=1:numel(im_cells_data)
        if im_cells_data(n).Area>700 && im_cells_data(n).Area<35000 && sum(A(:,n)) == 0
            counter = counter + 1;
            Area(counter) = im_cells_data(n).Area;
            Perimeter(counter) = im_cells_data(n).Perimeter;
            D = zeros(size(B{n},1),1);
            for i = 1:size(B{n},1)
                D(i) = sqrt((B{n}(i,1)-im_cells_data(n).Centroid(2))*(B{n}(i,1)-im_cells_data(n).Centroid(2)) + ...
                    (B{n}(i,2)-im_cells_data(n).Centroid(1))*(B{n}(i,2)-im_cells_data(n).Centroid(1)));
            end
            MeanBD(counter) = mean(D);
            MaxBD(counter) = max(D);
            Major(counter) = im_cells_data(n).MajorAxisLength;
            Minor(counter) = im_cells_data(n).MinorAxisLength;
            Circularity(counter) = im_cells_data(n).Perimeter /...
                (2*sqrt(pi * im_cells_data(n).Area));
            Solidity(counter) = im_cells_data(n).Solidity;
            Ecc(counter) = im_cells_data(n).Eccentricity;
            AR(counter) = im_cells_data(n).MajorAxisLength / im_cells_data(n).MinorAxisLength;
            
            I3 = imcomplement(I2) .* ...
                poly2mask(B{n}(:,2), B{n}(:,1), im_x, im_y);
            dilatedImage = imdilate(I3,strel('disk',10));
            I4 = bwmorph(dilatedImage,'thin',10);
            Cor = detectMinEigenFeatures(I4, 'MinQuality', 0.03,'FilterSize', 11);
            pgon = delaunayTriangulation(double(Cor.Location(:,1)), double(Cor.Location(:,2)));
            [C,v] = convexHull(pgon);
            
            for t=2:length(C)-1
                u = [pgon.Points(C(t-1),1)-pgon.Points(C(t),1),...
                    pgon.Points(C(t-1),2)-pgon.Points(C(t),2)];
                v = [pgon.Points(C(t+1),1)-pgon.Points(C(t),1),...
                    pgon.Points(C(t+1),2)-pgon.Points(C(t),2)];
                ThetaInDegrees{counter}(t-1) = acosd(dot(u, v) / norm(u) / norm(v)); %atan2d(u(1)*v(2)-u(2)*v(1),u(1)*u(2)+v(1)*v(2));
            end
            ThetaInDegrees{counter}(ThetaInDegrees{counter}>150) = [];
        end
    end
    
end

all = [Area', Perimeter', MeanBD', MaxBD', Major', Minor', Circularity', Solidity', Ecc', AR'];

[R, P] = corrcoef(all);


all2 = array2table(all);
all2.Properties.VariableNames = {'Area', 'Perimeter', 'MeanBD', 'MaxBD', 'MajorAxis',...
    'Minor', 'Circularity', 'Solidity', 'Ecc', 'AR'};
corrplot(all2);
cd(currdir);