%% Making a distribution
X = 3; % binning interval

cd(sum_dir);
distribution = zeros(250,250);
distances_angle = struct([]);
binrange = 1:100;
bincenter=binrange(1:end) - 0.5;
N= zeros(360/X, 100);
curve = struct([]);
gof = struct([]);

radius = zeros(1, 360/X);
sigma1 = zeros(1, 360/X);
pvalue = zeros(1, 360/X);
intensity = zeros(1,360/X);

for l=1:counter
    
    for m=1:length(b_cells_rot2{l})
        distribution(round(b_cells_rot2{l}(m,2))+125,round(b_cells_rot2{l}(m,1))+125) = ...
            distribution(round(b_cells_rot2{l}(m,2))+125,round(b_cells_rot2{l}(m,1))+125)+1;
    end
end

%% Image of a circle
image2 = figure;
imshow(distribution, [min(distribution(:)) max(distribution(:))]);
print(image2, '-dtiff', '-r300', 'distribution.tif');
cd(currdir);


%% Radius variation measuring angles and distances of each border point
angle = zeros(1,1);
distance = zeros(1,1);
p=0;
for l=1:counter
    for m=1:length(b_cells_rot2{l})
        p=p+1;
        distance(p) = sqrt(b_cells_rot2{l}(m,2)*b_cells_rot2{l}(m,2) +...
            b_cells_rot2{l}(m,1)*b_cells_rot2{l}(m,1));
        
        a = b_cells_rot2{l}(m,2)/distance(p);
        b = b_cells_rot2{l}(m,1)/distance(p);
        angle(p)= atan2d(a,b);
        if angle(p)<0
            angle(p)= angle(p)+360;
        end
    end
end


%% obtaining radial distribution by angles
for j=X:X:360
    distances_angle{j}= zeros(1,1);
    for q=1:length(distance)
        if angle(q)<=j && angle(q)>j-X
            distances_angle{j/X}=[distances_angle{j/X};distance(q)];
        end  
    end
    distances_angle{j/X} =  distances_angle{j/X}( distances_angle{j/X}~=0);
end

%% fitting radius
for i=X:X:360
    [N(i/X,:), bins] = histc(distances_angle{i/X},binrange);
    [curve{i/X},gof{i/X}] = fit(bincenter',N(i/X,:)','gauss1');
    radius(i/X) = curve{i/X}.b1;
    sigma1(i/X) = sqrt(curve{i/X}.c1*curve{i/X}.c1/2);
    pvalue(i/X) = gof{i/X}.rsquare;
end


%% Radius of borders distribution - radius
image3=figure;
cd(sum_dir);
plot(X:X:360, radius, 'LineWidth', 2);
set(gca,'FontSize',18);
ylim([0 max(radius)*1.1]);
xlim([0 400]);
print(image3, '-dtiff', '-r300', 'circle_radius.tif');

%% Width of borders disribution - sigma1
image4=figure;
plot(X:X:360, sigma1, 'LineWidth', 2);
set(gca,'FontSize',18);
ylim([0 max(sigma1)*1.1]);
xlim([0 400]);
print(image4, '-dtiff', '-r300', 'circle_width.tif');


%% Determening and plotting intensity of border distribution within 3 bins at maximum - intensity
for i=X:X:360
    intensity(i/X) = N(i/X,round(radius(i/X))-1) + N(i/X,round(radius(i/X))) + N(i/X,round(radius(i/X))+1);
end

image5=figure;
plot(X:X:360,intensity, 'LineWidth', 2);
set(gca,'FontSize',18);
ylim([0 max(intensity)*1.1]);
xlim([0 400]);
print(image5, '-dtiff', '-r300', 'circle_intensity.tif');

%% Fitting and plotting an ellipse to the radius
Coordinates = [(radius.*cosd(X:X:360))', (radius.*sind(X:X:360))'];
[Z_radius, a_radius, b_radius, alpha_radius] = fitellipse(Coordinates);
image6 = figure;
while alpha_radius<0
    alpha_radius = alpha_radius+pi;
end
Q_el = [cos(alpha_radius), -sin(alpha_radius); sin(alpha_radius) cos(alpha_radius)];
% Ellipse points
Coordinates_el = (Q_el * [a_radius * cosd(0:0.1:360); b_radius * sind(0:0.1:360)] + repmat(Z_radius, 1, length(0:0.1:360)))';

IDX = knnsearch(Coordinates_el,Coordinates);
Coordinates_el2 = Coordinates_el(IDX,:); 

% Determening GOF of ellipse
%% Getting GOF (R^2) of circle fit
SStot = 0;
SSres = 0;
for k=1:length(Coordinates)
    SStot = SStot + (Coordinates(k,1) - Z_radius(1)) * (Coordinates(k,1) - Z_radius(1))...
        + (Coordinates(k,2) - Z_radius(2)) * (Coordinates(k,2) - Z_radius(2));
    SSres = SSres + (sqrt((Coordinates(k,1) - Z_radius(1)) * (Coordinates(k,1) - Z_radius(1))...
        + (Coordinates(k,2) - Z_radius(2)) * (Coordinates(k,2) - Z_radius(2)))-...
        sqrt((Coordinates_el2(k,1) - Z_radius(1)) * (Coordinates_el2(k,1) - Z_radius(1))...
        + (Coordinates_el2(k,2) - Z_radius(2)) * (Coordinates_el2(k,2) - Z_radius(2))))^2;
end


R_square = 1-SSres/SStot;

% The actual plotting one-liner
plot(Coordinates_el(:,1), Coordinates_el(:,2), 'r', 'LineWidth',2);
hold on;
scatter(Coordinates(:,1),Coordinates(:,2), 'b', 'LineWidth',2);
hold on
annotation('textbox',[.4 .5 .0 .0],'String',['R_s_q_u_a_r_e = ', num2str(R_square)],'FitBoxToText','on','FontSize', 14);
hold off;
%close all
print(image6, '-dtiff', '-r300', 'circle_fit.tif');
