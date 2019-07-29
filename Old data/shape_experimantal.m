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
Orient = zeros(1,1);
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

byembryo = zeros(numel(files_tif),47);
borderlong = zeros(1,4);
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
    Orient2 = zeros(1,1);
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
            
            
            celldata;

        end
     end
    
    Orient2(Orient2<0) = 180 + Orient2(Orient2<0);
    for l=1:numel(B)
        if sum(A(:,l)) > 0
            Itemp = imdilate(poly2mask(B{l}(:,2),B{l}(:,1),im_x,im_y), strel('diamond', 1));
            I2 = I .* Itemp;
            I2 = bwareaopen(I2,100);
            [B2,L2, N2, A2] = bwboundaries(I2,'holes');
            
            cleanborder;
            
            Borders = bwconncomp(I3);
            im_borders_data=regionprops(Borders,'PixelList', 'Perimeter');
            V3 = imdilate(V, strel('diamond', 7));            
            Overlap = V3 .* I3;
            
            borderdata;
            
            
        end
    end
    
    byembryodata;
    
end



Orient(1) = [];
all = [Area', Perimeter', MeanBD', MaxBD', Major', Minor', Circularity',...
    Solidity', Ecc', AR', Orient, MeanTheta', MinTheta',  MaxTheta', NumTheta'];

[R, P] = corrcoef(all);


all2 = array2table(all);
all2.Properties.VariableNames = {'Area', 'Perimeter', 'MeanBD', 'MaxBD', 'MajorAxis',...
    'Minor', 'Circularity', 'Solidity', 'Ecc', 'AR', 'DevfromDV',...
    'MeanTheta', 'MinTheta', 'MaxTheta', 'NumTheta'};
image = figure;
corrplot(all2);

byembryo2 = array2table(byembryo);
byembryo2.Properties.VariableNames = {'Area', 'Areasd', 'Perimeter', 'Perimetersd',...
    'MeanBD', 'MeanBDsd', 'MaxBD', 'MaxBDsd', 'MajorAxis', 'MajorAxissd',...
    'Minor', 'Minorsd', 'Circularity', 'Circularitysd', ...
    'Solidity', 'Soliditysd', 'Ecc', 'Eccsd', 'AR', 'ARsd',...
    'Cell_direction', 'Cell_directionsd', 'DevfromDV', 'DEVfromDVsd',...
    'MeanTheta', 'MeanThetasd', 'MinTheta', 'MinThetasd',...
    'MaxTheta', 'MaxThetasd', 'NumTheta', 'NumThetasd', 'Number_cells',...
    'Straight_Length_40_90', 'Straight_Length_40_90sd', ...
    'Full_Length_40_90', 'Full_Length_40_90sd', 'Waviness_40_90',...
    'Waviness_40_90sd', 'Number_borders_40_90',...
    'Straight_Length_0_10', 'Straight_Length_0_10sd', ...
    'Full_Length_0_10', 'Full_Length_0_10sd', 'Waviness_0_10',...
    'Waviness_0_10sd', 'Number_borders_0_10'};


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


borderlong(1,:) = [];
borderlong(borderlong(:,4)==0,:) = [];
borderlong = sortrows(borderlong,4);
Bshort = allB(allB(:,4)>40,:);
borderlong2 = array2table(borderlong);
borderlong2.Properties.VariableNames = {'Straight_Length', 'Full_Length', 'Waviness', 'Angle'};
writetable(borderlong2,'Borders0-10.csv');
Bshort2 = array2table(Bshort);
Bshort2.Properties.VariableNames = {'Straight_Length', 'Full_Length', 'Waviness', 'Angle'};
writetable(Bshort2,'Borders40-90.csv');
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