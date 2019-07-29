

byembryo(g,1) = mean(Area2);
byembryo(g,2) = std(Area2);
byembryo(g,3) = mean(Perimeter2);
byembryo(g,4) = std(Perimeter2);
byembryo(g,5) = mean(MeanBD2);
byembryo(g,6) = std(MeanBD2);
byembryo(g,7) = mean(MaxBD2);
byembryo(g,8) = std(MaxBD2);
byembryo(g,9) = mean(Major2);
byembryo(g,10) = std(Major2);
byembryo(g,11) = mean(Minor2);
byembryo(g,12) = std(Minor2);
byembryo(g,13) = mean(Circularity2);
byembryo(g,14) = std(Circularity2);
byembryo(g,15) = mean(Solidity2);
byembryo(g,16) = std(Solidity2);
byembryo(g,17) = mean(Ecc2);
byembryo(g,18) = std(Ecc2);
byembryo(g,19) = mean(AR2);
byembryo(g,20) = std(AR2);

Orient2(Orient2<0) = 180 + Orient2(Orient2<0);
DEV = abs(Orient2 - mean(Orient2));
Orient = [Orient; DEV'];
byembryo(g,21) = mean(Orient2);
byembryo(g,22) = std(Orient2);

byembryo(g,23) = mean(DEV);
byembryo(g,24) = std(DEV);

byembryo(g,25) = mean(MeanTheta2);
byembryo(g,26) = std(MeanTheta2);
byembryo(g,27) = mean(MinTheta2);
byembryo(g,28) = std(MinTheta2);
byembryo(g,29) = mean(MaxTheta2);
byembryo(g,30) = std(MaxTheta2);
byembryo(g,31) = mean(NumTheta2);
byembryo(g,32) = std(NumTheta2);
byembryo(g,33) = counter3;

% borders 40-90°
byembryo(g,34) = mean(Length_S2(B_angle2>40));
byembryo(g,35) = std(Length_S2(B_angle2>40));
byembryo(g,36) = mean(Length_F2(B_angle2>40));
byembryo(g,37) = std(Length_F2(B_angle2>40));
byembryo(g,38) = mean(Wave2(B_angle2>40));
byembryo(g,39) = std(Wave2(B_angle2>40));
byembryo(g,40) = sum(B_angle2>40);

% borders 0-10°
B_angle2(Length_S2<byembryo(g,9)/3) = [];
Length_F2(Length_S2<byembryo(g,9)/3) = [];
Wave2(Length_S2<byembryo(g,9)/3) = [];
Length_S2(Length_S2<byembryo(g,9)/3) = [];
borderlong = [borderlong; [Length_S2',Length_F2',Wave2',B_angle2']];
byembryo(g,41) = mean(Length_S2(B_angle2<10));
byembryo(g,42) = std(Length_S2(B_angle2<10));
byembryo(g,43) = mean(Length_F2(B_angle2<10));
byembryo(g,44) = std(Length_F2(B_angle2<10));
byembryo(g,45) = mean(Wave2(B_angle2<10));
byembryo(g,46) = std(Wave2(B_angle2<10));
byembryo(g,47) = sum(B_angle2<10);
