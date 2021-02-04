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
if exist([filedir, '/summary_shape'],'dir') == 0
    mkdir(filedir,'/summary_shape');
end
sum_dir = [filedir, '/summary_shape'];

if exist([filedir, '/cell_images'],'dir') == 0
    mkdir(filedir,'/cell_images');
end
cell_dir = [filedir, '/cell_images'];


cd(tif8_dir);
files_tif = dir('*.tif');
se90V = strel('line', 2, 90);
se0V = strel('line', 2, 0);
ThetaInDegrees = struct([]);

byembryo = zeros(numel(files_tif),5);
all = zeros(1,2);
byembryoAR = zeros(numel(files_tif),5);
allAR = zeros(1,2);

byembryoout = zeros(numel(files_tif),5);
allout = zeros(1,2);
byembryooutAR = zeros(numel(files_tif),5);
alloutAR = zeros(1,2);

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
    counter = 0;
    Ecc = zeros(1,1);
    image1=figure;
    imshow(I);
    hold on;
    for l=1:numel(B)
        if sum(A(:,l)) > 0
            Itemp = imdilate(poly2mask(B{l}(:,2),B{l}(:,1),im_x,im_y), strel('diamond', 1));
            I2 = I .* Itemp;
            I2 = bwareaopen(I2,100);
            [B2,L2, N2, A2] = bwboundaries(I2,'holes');
            im_cells_data=regionprops(L2,'Eccentricity', 'Area','Centroid');
            
            
            for n=1:numel(im_cells_data)
                if im_cells_data(n).Area>700 && im_cells_data(n).Area<35000 && sum(A2(:,n)) == 0
                    counter = counter + 1;
                    Ecc(counter) = im_cells_data(n).Eccentricity;
                    
                    boundary = B2{n};
                    c = im_cells_data(n).Centroid;
                    c_labels = text(c(1), c(2), sprintf('%d', counter),'HorizontalAlignment', 'center',...
                        'VerticalAlignment', 'middle', 'Fontsize', 10);
                    set(c_labels,'Color',[1 1 0])
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                    
                end
            end
            
        end
    end
    
    [Ecc_clean, TF] = rmoutliers(Ecc');
    AR = 1 ./ sqrt(1 - Ecc_clean .* Ecc_clean);
    ARout = 1 ./ sqrt(1 - Ecc' .* Ecc');
    counter = 0;
    for l=1:numel(B)
        if sum(A(:,l)) > 0
            Itemp = imdilate(poly2mask(B{l}(:,2),B{l}(:,1),im_x,im_y), strel('diamond', 1));
            I2 = I .* Itemp;
            I2 = bwareaopen(I2,100);
            [B2,L2, N2, A2] = bwboundaries(I2,'holes');
            im_cells_data=regionprops(L2,'Eccentricity', 'Area','Centroid');
            for n=1:numel(im_cells_data)
                if im_cells_data(n).Area>700 && im_cells_data(n).Area<35000 && sum(A2(:,n)) == 0
                    counter = counter + 1;
                    if TF(counter) == 1
                        boundary = B2{n};
                        c = im_cells_data(n).Centroid;
                        c_labels = text(c(1), c(2), sprintf('%d', counter),'HorizontalAlignment', 'center',...
                            'VerticalAlignment', 'middle', 'Fontsize', 10);
                        set(c_labels,'Color',[1 1 1])
                        plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', 2);
                        
                    end
                end
             end
            
        end
    end
    
    byembryo(g,:) = [g, mean(Ecc_clean),std(Ecc_clean),length(Ecc_clean),sum(TF)];
    byembryoout(g,:) = [g, mean(Ecc),std(Ecc),length(Ecc),sum(TF)];
    byembryoAR(g,:) = [g, mean(AR),std(AR),length(Ecc_clean),sum(TF)];
    byembryooutAR(g,:) = [g, mean(ARout),std(ARout),length(Ecc),sum(TF)];
    
    all = [all; [ones(length(Ecc_clean),1) * g, Ecc_clean]]; %#ok<AGROW>
    allAR = [allAR; [ones(length(Ecc_clean),1) * g, AR]]; %#ok<AGROW>
    cd(cell_dir);
    image_filename = [num2str(g),'_analysed_image.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    
end

all(1,:) = [];
allAR(1,:) = [];

all2 = array2table(all);
all2.Properties.VariableNames = {'Embryo', 'Eccentricity'};

all3 = array2table(allAR);
all3.Properties.VariableNames = {'Embryo', 'Eccentricity'};

byembryo2 = array2table(byembryo);
byembryo2.Properties.VariableNames = {'embryo', 'Eccentricity', 'SD',...
    'Number_cells', 'Number_outliers'};

byembryo3 = array2table(byembryoAR);
byembryo3.Properties.VariableNames = {'embryo', 'Eccentricity', 'SD',...
    'Number_cells', 'Number_outliers'};

byembryo4 = array2table(byembryoout);
byembryo4.Properties.VariableNames = {'embryo', 'Eccentricity', 'SD',...
    'Number_cells', 'Number_outliers'};
byembryo5 = array2table(byembryooutAR);
byembryo5.Properties.VariableNames = {'embryo', 'Eccentricity', 'SD',...
    'Number_cells', 'Number_outliers'};


cd(sum_dir);
writetable(all2,'cells.csv');
writetable(byembryo2,'embryos.csv');
writetable(all3,'cellsAR.csv');
writetable(byembryo3,'embryosAR.csv');
writetable(byembryo4,'embryos_noout.csv');
writetable(byembryo5,'embryosAR_noout.csv');
cd(currdir);

close all;
clear variables;
clc;