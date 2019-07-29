

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

        B_angle(counter2) = -atand((Ends(1).Centroid(1,2) - Ends(2).Centroid(1,2))/...
            (Ends(1).Centroid(1,1) - Ends(2).Centroid(1,1)));
        
        
        if abs(B_angle(counter2) - mean(Orient2))>90
            if mean(Orient2) > 0
                B_angle(counter2) = abs(mean(Orient2) - 180 - B_angle(counter2));
            else
                B_angle(counter2) = abs(mean(Orient2) + 180 - B_angle(counter2));
            end
        else
            B_angle(counter2) = abs(mean(Orient2) - B_angle(counter2));
        end
        B_angle2(counter4) = B_angle(counter2);
    end
end