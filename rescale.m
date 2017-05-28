
%% Get boundaries and move objects to center at 0-0 coordinates
b_cells_rot = struct([]);
b_cells_rot2 = struct([]);
corner_rot = struct([]);
corner_rot2 = struct([]);

image1 = figure;
for i=1:counter
    R = [cosd(Orientation(i)) -sind(Orientation(i)); sind(Orientation(i)) cosd(Orientation(i))];
    Center_mat_border = repmat([Center(i,1); Center(i,2)], 1, length(b_cells_flip{i}))';
    
    b_cells_rot{i} = b_cells_flip{i} - Center_mat_border;
    b_cells_rot{i} = (R*b_cells_rot{i}')';
    b_cells_rot2{i}(:,1) = b_cells_rot{i}(:,1)*100/MajorAxis(i);
    b_cells_rot2{i}(:,2) = b_cells_rot{i}(:,2)*100/MinorAxis(i);
    
    
    if length(corners{i})>2
        for j=1:length(corners{i})
            Center_mat_border2 = [Center(i,1), Center(i,2)];
            corner_rot{i}(j,:) = Corner_center{i}(j,:) - Center_mat_border2;
            corner_rot{i}(j,:) = (R*corner_rot{i}(j,:)')';
            corner_rot2{i}(j,1) = corner_rot{i}(j,1)*100/MajorAxis(i);
            corner_rot2{i}(j,2) = corner_rot{i}(j,2)*100/MinorAxis(i);
            corner_rot2{i}(j,1) = corner_rot2{i}(j,1)*corners_good{i}(j);
            corner_rot2{i}(j,2) = corner_rot2{i}(j,2)*corners_good{i}(j);
        end
    corner_rot2{i}(corner_rot2{i}(:,1)==0,:) = [];    
    end
    
    plot(b_cells_rot2{i}(:,1),b_cells_rot2{i}(:,2));
    hold on
end
hold off
cd(sum_dir);
print(image1, '-dtiff', '-r300', 'all_boundaries.tif');
cd(currdir);