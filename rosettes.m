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
mkdir(filedir,'/summary_rosettes');
sum_dir = [filedir, '/summary_rosettes'];
mkdir(filedir,'/rosettes');
res_dir = [filedir, '/rosettes'];

cd(tif8_dir);
files_tif = dir('*.tif');
vert2 = zeros(1,2);
pulled_data = zeros(1,5);
pulled_vertices = zeros(1,3);
averagedata = zeros(1,8);
se90V = strel('line', 6, 90);
se0V = strel('line', 6, 0);


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
    
    V=imread('vertices.tif');
    V=imbinarize(rgb2gray(V),0);
    V(:,1) = 0;
    V(:,end) = 0;
    V(1,:) = 0;
    V(end,:) = 0;
    V2 = imdilate(V,strel('disk', 8));
    V2 = imclearborder(V2);
    
    %% Cell shape from regionprops
    [B,L,N,A] = bwboundaries(I,'holes');
    im_cells_data=regionprops(L,'Area', 'Eccentricity','PixelList',...
        'MajorAxisLength', 'MinorAxisLength');
    counter = 0;
    area = zeros(1,1);
    eccentricity = zeros(1,1);
    AR = zeros(1,1);
    L(L==0) = numel(im_cells_data) + 1;
    for i=1:numel(im_cells_data)
        if sum(A(:,i)) > 0 || sum(A(i,:)) == 0
           for k=1:length(im_cells_data(i).PixelList)
                L(im_cells_data(i).PixelList(k,2),im_cells_data(i).PixelList(k,1)) = 0;
            end 
        else            
            counter = counter + 1;
            area(counter) = im_cells_data(i).Area;   
            eccentricity(counter) = im_cells_data(i).Eccentricity;  
            AR(counter) = im_cells_data(i).MajorAxisLength/im_cells_data(i).MinorAxisLength;
        end
    end
    image1 = figure;
    imshow(imdilate(I, strel('disk',3)));
    Vertices = bwconncomp(V2);
    im_vert_data=regionprops(Vertices,L,'PixelValues','Centroid');
    vert = zeros(1,1);
    counter2=0;
    for i=1:numel(im_vert_data)
        if ismember(numel(im_cells_data) + 1, unique(im_vert_data(i).PixelValues)) == 0
            counter2 = counter2 + 1;
            vert(counter2) = length(unique(im_vert_data(i).PixelValues(im_vert_data(i).PixelValues>0)));
            c = im_vert_data(i).Centroid;
            c_labels = text(c(1), c(2), sprintf('%d', vert(counter2)),...
                'HorizontalAlignment', 'center',...
                'VerticalAlignment', 'middle', 'Fontsize', 20);
            set(c_labels,'Color',[1 0 0])
        end
    end
    
    %% Saving the image
    cd(res_dir);
    image_filename = [num2str(g),'_analysed_image.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    
    embryodata = [(1:counter)',area', eccentricity', AR'];
    dataT = array2table(embryodata);
    dataT.Properties.VariableNames = {'cell', 'area', 'eccentricity', 'AR'};
    writetable(dataT, [num2str(g),'_cell_data.csv']);
    pulled_data = [pulled_data; [ones(counter,1)*g, embryodata]];
    pulled_vertices = [pulled_vertices; [ones(counter2,1)*g,(1:counter2)',vert']];
    
    averagedata(g,1) = g;
    averagedata(g,2:2:end-1) = mean(dataT{:,2:end},1);
    averagedata(g,3:2:end-1) = std(dataT{:,2:end},1);
    averagedata(g,end) = counter;
end

cd(sum_dir);

pulled_data(1,:) = [];
T1 = array2table(pulled_data);
T1.Properties.VariableNames = {'embryo','cell', 'area', 'eccentricity', 'AR'};
writetable(T1,'pulled_cell_shape.csv');

pulled_vertices(1,:) = [];
T2 = array2table(pulled_vertices);
T2.Properties.VariableNames = {'embryo','cell', 'neighbours'};
writetable(T2,'pulled_vertices.csv');

T3 = array2table(averagedata);
T3.Properties.VariableNames = {'embryo','area','areaSD',...
    'eccentricity', 'eccentricitySD', 'AR', 'ARSD','cells'};
writetable(T3,'average_cell_shape.csv');

HistV = table({'3';'4';'rosettes'}, [sum(pulled_vertices(:,3)==3);...
    sum(pulled_vertices(:,3)==4); sum(pulled_vertices(:,3)>4)],...
    [sum(pulled_vertices(:,3)==3)*100/size(pulled_vertices,1);...
    sum(pulled_vertices(:,3)==4)*100/size(pulled_vertices,1);...
    sum(pulled_vertices(:,3)>4)*100/size(pulled_vertices,1)]);
HistV.Properties.VariableNames = {'vertex','number','percent'};
writetable(HistV,'vertices_distribution.csv');

cd(currdir);
close all;
clear variables;
clc;