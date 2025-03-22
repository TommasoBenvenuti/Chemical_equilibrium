% Script molto basic per fitting lineare
clc; clear; close all;
clc; clear; close all;

% Lettura dati senza le prime 5 righe di intestazione
filename = 'risultati_kp.dat';
data = readmatrix(filename, 'NumHeaderLines', 5);

%  dati
T = data(:,1);      % Temperatura (K)
Kp = data(:,2);     % Costante di equilibrio
deltaG = data(:,3); % Energia libera di Gibbs (J/mol)

% Costante universale dei gas
R = 8.314; % J/(mol*K)

invT = 1 ./ T;       % Calcola 1/T
lnKp = log(Kp);      % Calcola ln(Kp)

% Fit lineare
coeffs = polyfit(invT, lnKp, 1);
fit_lnKp = polyval(coeffs, invT);

% Estrazione dei parametri
DeltaH = -coeffs(1) * 8.314; % J/mol
DeltaS = coeffs(2) * 8.314;  % J/mol·K

% Plot
figure;
plot(invT, lnKp, 'bo', 'MarkerFaceColor', 'b', 'DisplayName', 'Dati');
hold on;
plot(invT, fit_lnKp, 'black', 'LineWidth', 1.5, 'DisplayName', 'Fit lineare');
xlabel('1/T (K^{-1})');
ylabel('ln(K_p)');
title('Equazione di Van''t Hoff');% Definizione della legenda con i valori calcolati
legendatesto = sprintf('\\DeltaH = %.2f J/mol, \\DeltaS = %.2f J/(mol K)', DeltaH, DeltaS);
legend(legendatesto);
grid on;

fprintf('ΔH° = %.2f J/mol\n', DeltaH);
fprintf('ΔS° = %.2f J/mol·K\n', DeltaS);


% Plot di Kp vs Temperatura
figure;
semilogy(T, Kp, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b');
xlabel('Temperatura (K)');
ylabel('K_p');
title('Cost. vs T (semilog)');

% Plot di DeltaG vs Temperatura
figure;
plot(T, deltaG, 'r-', 'LineWidth', 1.5, 'Marker', 'o', 'MarkerFaceColor', 'b');
xlabel('Temperatura (K)');
ylabel('\DeltaG (J/mol)');
title('Delta G vs T');


