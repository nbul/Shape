
cd(tif8_dir);

%% Collecting all borders and corners
se90V = strel('line', 2, 90);
se0V = strel('line', 2, 0);
%dilation of borders
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);
Orientation = zeros(1,1);
b_cells = struct([]);
b_cells_flip = struct([]);
Corner_center = struct([]);
Center = zeros(2,2);
MajorAxis = zeros(1,1);
MinorAxis = zeros(1,1);
Pixels = struct([]);
corners = struct([]);
counter = 0;
counter2 = 0;

for g=1:numel(files_tif)
    %% Open images and modify
    clear s_cells L_cells
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
    I_cells=imcomplement(I2);
    [B,L] = bwboundaries(I_cells,4);
    B{1}=[];
    for i=1:5
        B{length(B)-i+1}=[];
        L(L==length(B)-i+1) = 0;
    end
    B = B(~cellfun('isempty',B));
    L = L-1;
    
    s_cells=regionprops(L, 'Eccentricity', 'MajorAxisLength', 'MinorAxisLength', 'Orientation','PixelList','Centroid', 'Area');
  
    for i=1:numel(s_cells)
        if s_cells(i).Area >200
            counter = counter+1;
            Pixels{counter} = s_cells(i).PixelList;
            Orientation(counter) = s_cells(i).Orientation;
            MajorAxis(counter) = s_cells(i).MajorAxisLength;
            MinorAxis(counter) = s_cells(i).MinorAxisLength;
            Center(counter,:) = round(s_cells(i).Centroid);
            b_cells{counter} = B{i};
            b_cells_flip{counter} = [ b_cells{counter}(:,2), b_cells{counter}(:,1)];
            %% Getting cornrers coordinates and selectin <120%
            for c=1:length(b_cells{counter})
                corners{counter}(c) = N_borders(b_cells{counter}(c,1),b_cells{counter}(c,2));
            end      
            
            corners{counter} = corners{counter}(corners{counter}~=0);
            corners{counter} = unique(corners{counter}, 'stable');
            
            for c=1:length(corners{counter});
                Corner_temp = round(s_borders(corners{counter}(c)).Centroid);
                Corner_center{counter}(c,1) = Corner_temp(1);
                Corner_center{counter}(c,2) = Corner_temp(2);  
            end
            
            if length(corners{counter})>2
                corners{counter} = [corners{counter},corners{counter}(1), corners{counter}(2)];
                Corner_center{counter} = [Corner_center{counter};Corner_center{counter}(1,:); Corner_center{counter}(2,:)];
            end   
        end
    end
end
% plot(Corner_center{100}(:,1),Corner_center{100}(:,2))
% hold on
% plot(b_cells_flip{100}(:,1),b_cells_flip{100}(:,2))
% 
% imshow(L ==i)
% hold on
% plot(B{i}(:,2),B{i}(:,1));

cd(currdir);