% Dati
x = [40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420];
y = [0.45957E+03, 0.59597E+02, 0.21548E+02, 0.11768E+02, 0.79065E+01, 0.59803E+01, 0.48697E+01, 0.41638E+01, 0.36827E+01, 0.33371E+01, 0.30787E+01, 0.28792E+01, 0.27210E+01, 0.25930E+01, 0.24873E+01, 0.23988E+01, 0.23237E+01, 0.22592E+01, 0.22033E+01, 0.21544E+01];

% Fitting esponenziale
f = fit(x', y', 'exp1'); 
% Parametri del fitting
a = f.a; % Coefficiente a
b = f.b; % Coefficiente b
disp('Parametri del fitting esponenziale:');
disp(['a = ', num2str(a)]);
disp(['b = ', num2str(b)]);

% Creazione del grafico
figure;
plot(x, y, 'o', 'MarkerSize', 8, 'Color', 'red'); % Dati "sperimentali"
hold on;
plot(f, x, y); % Curva di fitting
title('\textit{Grafico dei dati con fitting Exp.}', 'Interpreter', 'latex');
xlabel('\textit{Temperatura (K)}', 'Interpreter', 'latex');
ylabel('\textit{Costante equilibrio}', 'Interpreter', 'latex');
legend('\textit{Dati}', '\textit{Fitting Esponenziale}', 'Interpreter', 'latex'); % Legenda in corsivo
grid on; % Griglia
set(gca, 'FontSize', 12); % Dimensione font