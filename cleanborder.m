%collecting borders
V=imread('vertices.tif');
V=imbinarize(rgb2gray(V),0);
V2 = imdilate(V, strel('diamond', 5));

I4 = imcomplement(imdilate(I2,strel('diamond', 1)));
I4 = imclearborder(I4);
I4 = I4 - bwareaopen(I4,700);
I4 = imdilate(I4, strel('diamond', 3));

I2 = I2 - I4;
I2(I2<0) = 0;

for i=2:(im_y-1)
    k = false;
    j = 0;
    while k==false
        j= j+1;
        if I2(i,j) == 1
            I2(i,j) = 0;
            k = true;
        elseif j== im_x
            k = true;
        end
    end
    j = im_x;
    k= false;
    while k==false
        j= j-1;
        if I2(i,j) == 1
            I2(i,j) = 0;
            k = true;
        elseif j== 1
            k = true;
        end
    end
end
I2 = bwareaopen(I2,100);
I3 = imsubtract(I2,V2);
I3 = imdilate(I3, strel('diamond', 1));
I3(I3<0) = 0;
I3 = bwareaopen(I3,20);