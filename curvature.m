

Coordinates = [(radius.*cosd(X:X:360))', (radius.*sind(X:X:360))'];

Coordinates = [Coordinates(end-2,:); Coordinates(end-1,:); Coordinates(end,:); Coordinates; Coordinates(1,:);Coordinates(2,:);Coordinates(3,:)];


 for i=4:(length(Coordinates)-3)
     x02 = Coordinates(i-3,1);
     x0 = Coordinates(i-2,1);
     x1 = Coordinates(i-1,1);
     x2 = Coordinates(i,1);
     x3 = Coordinates(i+1,1);
     x4 = Coordinates(i+2,1);
     x5 = Coordinates(i+3,1);
     y02 = Coordinates(i-3,2);
     y0 = Coordinates(i-2,2);
     y1 = Coordinates(i-1,2);
     y2 = Coordinates(i,2);
     y3 = Coordinates(i+1,2);
     y4 = Coordinates(i+2,2);
     y5 = Coordinates(i+3,2);
    [centerX(i-3),centerY(i-3),Rfit(i-3)] = circfit([x02,x0,x1,x2,x3,x4,x5],[y02,y0,y1, y2,y3,y4,y5]);
 end


cd(sum_dir);
image10=figure;
plot(Rfit,number_corners);
print(image10, '-dtiff', '-r300', 'curvature_vs_corners.tif');