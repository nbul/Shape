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
MeanTheta = zeros(1,1);
MinTheta = zeros(1,1);
MaxTheta = zeros(1,1);
NumTheta = zeros(1,1);



%Borders data
Length_S = zeros(1,1);
Length_F = zeros(1,1);
Wave = zeros(1,1);
B_angle = zeros(1,1);


cd(tif8_dir);
files_tif = dir('*.tif');
counter = 0;
counter2 = 0;
se90V = strel('line', 2, 90);
se0V = strel('line', 2, 0);
ThetaInDegrees = struct([]);

byembryo = zeros(numel(files_tif),36);

for g=1:numel(files_tif)
    
    bd_dir = [tif8_dir,'/', num2str(g)];
    cd(bd_dir);
    I=imread('handCorrection.tif');
    I=imbinarize(rgb2gray(I),0);
    I(:,1) = 0;
    I(:,end) = 0;
    I(1,:) = 0;
    I(end,:) = 0;
    [im_x, im_y] = size(I);
    
    %% Cell shape from regionprops
    [B,L,N,A] = bwboundaries(I,'holes');
    counter3=0;
    counter4=0;
    Area2 = zeros(1,1);
    Perimeter2 = zeros(1,1);
    MeanBD2 = zeros(1,1);
    MaxBD2 = zeros(1,1);
    Major2 = zeros(1,1);
    Minor2 = zeros(1,1);
    Circularity2 = zeros(1,1);
    Solidity2 = zeros(1,1);
    Ecc2 = zeros(1,1);
    AR2 = zeros(1,1);
    MeanTheta2 = zeros(1,1);
    MinTheta2 = zeros(1,1);
    MaxTheta2 = zeros(1,1);
    NumTheta2 = zeros(1,1);
    Length_S2 = zeros(1,1);
    Length_F2 = zeros(1,1);
    Wave2 = zeros(1,1);
    B_angle2 = zeros(1,1);
    
    for l=1:numel(B)
        if sum(A(:,l)) > 0
            Itemp = imdilate(poly2mask(B{l}(:,2),B{l}(:,1),im_x,im_y), strel('diamond', 1));
            I2 = I .* Itemp;
            I2 = bwareaopen(I2,100);
            [B2,L2, N2, A2] = bwboundaries(I2,'holes');
            im_cells_data=regionprops(L2,'Centroid', 'Area', 'Eccentricity','ConvexArea',...
                'MajorAxisLength', 'MinorAxisLength','Perimeter', 'Solidity','PixelList','Orientation');
            for n=1:numel(im_cells_data)
                if im_cells_data(n).Area>700 && im_cells_data(n).Area<35000 && sum(A2(:,n)) == 0
                    counter = counter + 1;
                    counter3 = counter3 + 1;
                    
                    Area(counter) = im_cells_data(n).Area;
                    Area2(counter3) = im_cells_data(n).Area;
                    Perimeter(counter) = im_cells_data(n).Perimeter;
                    Perimeter2(counter3) = im_cells_data(n).Perimeter;
                    D = zeros(size(B2{n},1),1);
                    for i = 1:size(B2{n},1)
                        D(i) = sqrt((B2{n}(i,1)-im_cells_data(n).Centroid(2))*(B2{n}(i,1)-im_cells_data(n).Centroid(2)) + ...
                            (B2{n}(i,2)-im_cells_data(n).Centroid(1))*(B2{n}(i,2)-im_cells_data(n).Centroid(1)));
                    end
                    MeanBD(counter) = mean(D);
                    MeanBD2(counter3) = mean(D);
                    MaxBD(counter) = max(D);
                    MaxBD2(counter3) = max(D);
                    Major(counter) = im_cells_data(n).MajorAxisLength;
                    Major2(counter3) = im_cells_data(n).MajorAxisLength;
                    Minor(counter) = im_cells_data(n).MinorAxisLength;
                    Minor2(counter3) = im_cells_data(n).MinorAxisLength;
                    Circularity(counter) = im_cells_data(n).Perimeter /...
                        (2*sqrt(pi * im_cells_data(n).Area));
                    Circularity2(counter3) = im_cells_data(n).Perimeter /...
                        (2*sqrt(pi * im_cells_data(n).Area));
                    Solidity(counter) = im_cells_data(n).Solidity;
                    Solidity2(counter3) = im_cells_data(n).Solidity;
                    Ecc(counter) = im_cells_data(n).Eccentricity;
                    Ecc2(counter3) = im_cells_data(n).Eccentricity;
                    AR(counter) = im_cells_data(n).MajorAxisLength / im_cells_data(n).MinorAxisLength;
                    AR2(counter3) = im_cells_data(n).MajorAxisLength / im_cells_data(n).MinorAxisLength;
                    
                    I3 = imcomplement(I2) .* ...
                        poly2mask(B2{n}(:,2), B2{n}(:,1), im_x, im_y);
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
                    MeanTheta(counter) = mean(ThetaInDegrees{counter});
                    MeanTheta2(counter3) = MeanTheta(counter);
                    MinTheta(counter) = min(ThetaInDegrees{counter});
                    MinTheta2(counter3) = MinTheta(counter);
                    MaxTheta(counter) = max(ThetaInDegrees{counter});
                    MaxTheta2(counter3) = MaxTheta(counter);
                    NumTheta(counter) = length(ThetaInDegrees{counter});
                    NumTheta2(counter3) = NumTheta(counter);
                end
            end
            %collecting borders
            V=imread('vertices.tif');
            V=imbinarize(rgb2gray(V),0);
            V2 = imdilate(V, strel('diamond', 5));
            
            for i=2:(im_y-1)
                k = false;
                j = 0;
                while k==false
                    j= j+1;
                    if I2(i,j) == 1
                        I2(i,j) = 0;
                        k = true;
                    elseif j== im_x
                        k = true;
                    end
                end
                j = im_x;
                k= false;
                while k==false
                    j= j-1;
                    if I2(i,j) == 1
                        I2(i,j) = 0;
                        k = true;
                    elseif j== 1
                        k = true;
                    end
                end
            end
            I2 = bwareaopen(I2,100);
            I3 = imsubtract(I2,V2);
            I3 = imdilate(I3, strel('diamond', 1));
            I3(I3<0) = 0;
            I3 = bwareaopen(I3,20);
            Orient = [im_cells_data.Orientation];
            Orient(Orient<0) = 180 + Orient(Orient<0);
            Orient_mean = mean(Orient);
            Borders = bwconncomp(I3);
            im_borders_data=regionprops(Borders,'PixelList', 'Perimeter');
            V3 = imdilate(V, strel('diamond', 7));
            
            Overlap = V3 .* I3;
            for i=1:numel(im_borders_data)
                temp = zeros(im_x,im_y);
                for j = 1:length(im_borders_data(i).PixelList)
                    temp(im_borders_data(i).PixelList(j,2),im_borders_data(i).PixelList(j,1)) = 1;
                end
                temp2 = Overlap .* temp;
                Ends2 = bwconncomp(temp2);
                Ends = regionprops(Ends2, 'Centroid');
                
                if numel(Ends) == 2
                    counter2 = counter2 + 1;
                    counter4 = counter4 + 1;
                    Length_S(counter2) = sqrt(((Ends(1).Centroid(1,1) - Ends(2).Centroid(1,1)) .*...
                        (Ends(1).Centroid(1,1) - Ends(2).Centroid(1,1))) + ...
                        ((Ends(1).Centroid(1,2) - Ends(2).Centroid(1,2)) .*...
                        (Ends(1).Centroid(1,2) - Ends(2).Centroid(1,2))));
                    Length_S2(counter4) = Length_S(counter2);
                    Length_F(counter2) = im_borders_data(i).Perimeter/2;
                    Length_F2(counter4) = Length_F(counter2);
                    Wave(counter2) = Length_F(counter2) / Length_S(counter2);
                    Wave2(counter4) = Wave(counter2);
                    B_angle(counter2) = atand((Ends(1).Centroid(1,1) - Ends(2).Centroid(1,1))/...
                        (Ends(1).Centroid(1,2) - Ends(2).Centroid(1,2)));
                    B_angle2(counter4) = B_angle(counter2);
                end
            end
        end
    end
    byembryo(g,1) = mean(Area2);
    byembryo(g,2) = std(Area2);
    byembryo(g,3) = mean(Perimeter2);
    byembryo(g,4) = std(Perimeter2);
    byembryo(g,5) = mean(MeanBD2);
    byembryo(g,6) = std(MeanBD2);
    byembryo(g,7) = mean(MaxBD2);
    byembryo(g,8) = std(MaxBD2);
    byembryo(g,9) = mean(Major2);
    byembryo(g,10) = std(Major2);
    byembryo(g,11) = mean(Minor2);
    byembryo(g,12) = std(Minor2);
    byembryo(g,13) = mean(Circularity2);
    byembryo(g,14) = std(Circularity2);
    byembryo(g,15) = mean(Solidity2);
    byembryo(g,16) = std(Solidity2);
    byembryo(g,17) = mean(Ecc2);
    byembryo(g,18) = std(Ecc2);
    byembryo(g,19) = mean(AR2);
    byembryo(g,20) = std(AR2);
    byembryo(g,21) = mean(MeanTheta2);
    byembryo(g,22) = std(MeanTheta2);
    byembryo(g,23) = mean(MinTheta2);
    byembryo(g,24) = std(MinTheta2);
    byembryo(g,25) = mean(MaxTheta2);
    byembryo(g,26) = std(MaxTheta2);
    byembryo(g,27) = mean(NumTheta2);
    byembryo(g,28) = std(NumTheta2);
    byembryo(g,29) = mean(Length_S2(B_angle2>40));
    byembryo(g,30) = std(Length_S2(B_angle2>40));
    byembryo(g,31) = mean(Length_F2(B_angle2>40));
    byembryo(g,32) = std(Length_F2(B_angle2>40));
    byembryo(g,33) = mean(Wave2(B_angle2>40));
    byembryo(g,34) = std(Wave2(B_angle2>40));
    byembryo(g,35) = counter3;
    byembryo(g,36) = sum(B_angle2>40);
end

all = [Area', Perimeter', MeanBD', MaxBD', Major', Minor', Circularity',...
    Solidity', Ecc', AR', MeanTheta', MinTheta',  MaxTheta', NumTheta'];

[R, P] = corrcoef(all);


all2 = array2table(all);
all2.Properties.VariableNames = {'Area', 'Perimeter', 'MeanBD', 'MaxBD', 'MajorAxis',...
    'Minor', 'Circularity', 'Solidity', 'Ecc', 'AR',...
    'MeanTheta', 'MinTheta', 'MaxTheta', 'NumTheta'};
image = figure;
corrplot(all2);

byembryo2 = array2table(byembryo);
byembryo2.Properties.VariableNames = {'Area', 'Areasd', 'Perimeter', 'Perimetersd',...
    'MeanBD', 'MeanBDsd', 'MaxBD', 'MaxBDsd', 'MajorAxis', 'MajorAxissd',...
    'Minor', 'Minorsd', 'Circularity', 'Circularitysd', ...
    'Solidity', 'Soliditysd', 'Ecc', 'Eccsd', 'AR', 'ARsd',...
    'MeanTheta', 'MeanThetasd', 'MinTheta', 'MinThetasd',...
    'MaxTheta', 'MaxThetasd', 'NumTheta', 'NumThetasd',...
    'Straight_Length_40_90', 'Straight_Length_40_90sd', ...
    'Full_Length_40_90', 'Full_Length_40_90sd', 'Waviness_40_90',...
    'Waviness_40_90sd', 'Number_cells', 'Number_borders_40_90'};


cd(sum_dir);
writetable(all2,'cells.csv');
saveas(image,'cells.tif');
writetable(byembryo2,'embryos.csv');

B_angle(B_angle<0) = -B_angle(B_angle<0);
allB = [Length_S', Length_F',Wave', B_angle'];
allB = sortrows(allB,4);
image2 = figure;
plot(allB(:,4),allB(:,3)-1,'o');
saveas(image2,'borders.tif');
allB2 = array2table(allB);
allB2.Properties.VariableNames = {'Straight_Length', 'Full_Length', 'Waviness', 'Angle'};
writetable(allB2,'Borders.csv');

cd(currdir);

% all3 = [Area', Perimeter',Minor', Circularity',...
%     Solidity', AR', MeanTheta', MinTheta',  MaxTheta', NumTheta'];
%
%
% Y = pdist(all,'cityblock');
% Z = linkage(Y,'average');
% c = cophenet(Z,Y);
% T = cluster(Z,'maxclust',3);
%
%
% figure;
% color = ['r','g','b','k','m'];
% for b = 1:5
%     scatter3(all3(T==b,1),all3(T==b,2),all3(T==b,6),color(b))
%     hold on
% end

close all;
clear variables;
clc;