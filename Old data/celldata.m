
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
        Orient2(counter3) = im_cells_data(n).Orientation;

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
