% Dati: Temperatura (K), Kp e deltaG
T = [50 70 90 110 130 150 170 190 210 230 250 270 290 310 330 350 370 390 410 430 450 470 490 510 530 550 570 590 610 630 650 670 690 710 730 750 770 790 810 830];
Kp = [134.855 33.297 15.384 9.466 6.798 5.357 4.481 3.902 3.497 3.199 2.973 2.796 2.654 2.538 2.441 2.360 2.290 2.230 2.178 2.132 2.091 2.055 2.023 1.994 1.967 1.943 1.921 1.901 1.882 1.865 1.849 1.834 1.821 1.808 1.796 1.785 1.774 1.764 1.755 1.746];
deltaG = [-2037.695 -2039.131 -2044.270 -2054.598 -2070.604 -2092.127 -2118.682 -2149.674 -2184.504 -2222.611 -2263.506 -2306.763 -2352.029 -2399.004 -2447.444 -2497.133 -2547.911 -2599.624 -2652.157 -2705.402 -2759.275 -2813.705 -2868.625 -2923.989 -2979.742 -3035.841 -3092.258 -3148.961 -3205.917 -3263.106 -3320.507 -3378.105 -3435.878 -3493.803 -3551.883 -3610.099 -3668.428 -3726.878 -3785.430 -3844.074];

figure;  % Crea una nuova finestra
plot(T, Kp, 'bo-', 'MarkerFaceColor', 'b');  % Dati e grafico
xlabel('Temperatura (K)');
ylabel('Costante di equilibrio');
title('Costante di equilibrio vs Temperatura');


% Seconda figura: deltaG vs Temperatura
figure; % Crea una nuova finestra per il secondo grafico
plot(T, deltaG, 'ro-', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
xlabel('Temperatura (K)');
ylabel('\Delta G (J/mol)');
title('Energia Libera di Gibbs vs Temperatura');
grid on;
