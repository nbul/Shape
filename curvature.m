
%% Curvature based on fitting a circle to the neighbourhood of 5°
Coordinates = [Coordinates(end-2,:); Coordinates(end-1,:); Coordinates(end,:); Coordinates; Coordinates(1,:);Coordinates(2,:);Coordinates(3,:)];
centerX = zeros(length(radius),1);
centerY = zeros(length(radius),1);
Rfit = zeros(length(radius),1);

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
image11=figure;
Rfit_vs_corner = [Rfit, number_corners'];
Rfit_vs_corner = sortrows(Rfit_vs_corner,2);
scatter(Rfit_vs_corner(:,1),Rfit_vs_corner(:,2));
set(gca,'FontSize',18);
ylim([0 max(Rfit_vs_corner(:,2))*1.1]);
xlim([0 max(Rfit_vs_corner(:,1))*1.1]);
print(image11, '-dtiff', '-r300', 'curvature_vs_corners.tif');

%% Curvature based on Lyuba's method
image12=figure;
K = max(radius)/min(radius)*radius ./...
    ((sind(90+(X:X:360)).^2 + (max(radius)/min(radius))^2 * cosd(90+(X:X:360)).^2) .^(3/2));
K_vs_corner = [K', number_corners'];
K_vs_corner = sortrows(K_vs_corner,1);

[f_kc,gof_kc] = fit(K_vs_corner(:,1),K_vs_corner(:,2),'poly1');

K_vs_corner_fit = f_kc.p1*(min(K_vs_corner(:,1)):1:max(K_vs_corner(:,1))) + f_kc.p2;
plot((min(K_vs_corner(:,1)):1:max(K_vs_corner(:,1))),K_vs_corner_fit, 'r', 'LineWidth',2);
hold on;
scatter(K_vs_corner(:,1),K_vs_corner(:,2));
set(gca,'FontSize',18);
ylim([0 max(K_vs_corner(:,2))*1.1]);
xlim([0 max(K_vs_corner(:,1))*1.1]);
hold on
annotation('textbox',[.2 .8 .0 .0],'String',{['R_s_q_u_a_r_e = ', num2str(gof_kc.rsquare)],...
    ['f = ', num2str(f_kc.p1), ' * x + ', num2str(f_kc.p2)] },'FitBoxToText','on','FontSize', 14);
hold off;
print(image12, '-dtiff', '-r300', 'curvatureL_vs_corners.tif');

print(image12, '-dtiff', '-r300', 'curvature2_vs_corners.tif');


