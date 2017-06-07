

Coordinates = [(radius.*cosd(X:X:360))', (radius.*sind(X:X:360))'];

Coordinates = [Coordinates(end,:); Coordinates; Coordinates(1,:)];


 for i=2:(length(Coordinates)-1)
     x1 = Coordinates(i-1,1);
     x2 = Coordinates(i,1);
     x3 = Coordinates(i+1,1);
     y1 = Coordinates(i-1,2);
     y2 = Coordinates(i,2);
     y3 = Coordinates(i+1,2);
     K(i-1) = 2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)) ./ ...
     sqrt(((x2-x1).^2+(y2-y1).^2)*((x3-x1).^2+(y3-y1).^2)*((x3-x2).^2+(y3-y2).^2));
 end
