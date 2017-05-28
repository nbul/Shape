%% selecting good corners (<150°)
corners_good = struct([]);
for i=1:counter
    if length(corners{i})>2 
        CosTheta = zeros(1,length(corners{i})-2);
        ThetaInDegrees = zeros(1,length(corners{i})-2);
        for c=2:length(corners{i})-1
            vector1 = Corner_center{i}(c+1,:) - Corner_center{i}(c,:);
            vector2 = Corner_center{i}(c-1,:) - Corner_center{i}(c,:);
            CosTheta(c-1) = dot( vector1, vector2)/(norm( vector1)*norm( vector2));
            ThetaInDegrees(c-1) = acosd(CosTheta(c-1));
            if ThetaInDegrees(c-1)>0 && ThetaInDegrees(c-1)<150
                corners_good{i}(c) = 1;
            else
                corners_good{i}(c) = 0;
            end
        end
        corners_good{i} = [corners_good{i},0];
    end
end