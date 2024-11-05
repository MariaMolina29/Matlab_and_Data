%% Parte 1: Graficar Pitch e Intensidad desde el archivo WAV

% Cargar el archivo de audio
[file_audio, path_audio] = uigetfile('*.wav', 'Selecciona un archivo WAV');
fullPath_audio = fullfile(path_audio, file_audio);
[snd, Fs] = audioread(fullPath_audio, 'native');

data_type = class(snd);
max_value = 32767;

if size(snd, 2) > 1
    snd = mean(double(snd), 2);  % Convertir a double para cálculos y promedio de canales
else
    snd = double(snd);  % Convertir a double para cálculos
end

if max_value > 0
    amplitud_normalizada = snd / max_value; % Normalizar
else
    amplitud_normalizada = snd; % Si max_value es 0, no se hace nada
end

tipo_snd = class(snd);


% Tiempo del audio
t = (0:length(snd)-1) / Fs;

%% Calcular Intensidad y Pitch desde el archivo WAV

window_length = round(0.05 * Fs);  % Ventana de 50 ms
intensity_wav = movmean(snd.^2, window_length);  % Calcular energía promedio en la ventana

% Convertir la intensidad a dB
intensity_dB_wav = 10 * log10(intensity_wav);  % eps es un valor pequeño para evitar log(0)
intensity_dB_wav(isinf(intensity_dB_wav)) = -300;


pitch_floor = 75;
pitch_ceiling = 600;


window_length = round(0.05 * Fs);  % Longitud de la ventana de 50 ms

% Calcular el pitch con el rango de [50 Hz, 500 Hz] usando el método 'NCF'
[pitch_values_wav, pitch_time_wav] = pitch(amplitud_normalizada, Fs, 'Method', 'NCF', 'Range', [pitch_floor, pitch_ceiling], 'WindowLength', window_length);

% Normalizar el eje de tiempo del pitch
pitch_time_wav = linspace(0, length(snd)/Fs, length(pitch_values_wav));

% Interpolar el pitch para que coincida en número de muestras con la intensidad
pitch_values_interp_wav = interp1(pitch_time_wav, pitch_values_wav, t, 'linear');

% Definir un umbral de intensidad para detectar silencio
intensity_threshold = 10;
pitch_values_interp_wav(intensity_dB_wav < intensity_threshold) = NaN;

%% Parte 2:  Pitch e Intensidad desde el archivo TXT

[file_txt, path_txt] = uigetfile('*.txt', 'Selecciona formants.txt');
fullPath_txt = fullfile(path_txt, file_txt);

% Leer el archivo de texto como matriz
data = readmatrix(fullPath_txt, 'Delimiter', '/');

% Acceder a las columnas usando índices
time_txtp = data(:, 1);       % Primera columna: Tiempo
pitch_txt = data(:, 2);      % Segunda columna: Frecuencia
intensity_txt = data(:, 3);  % Tercera columna: Intensidad
formant1_txt = data(:, 4);   % Cuarta columna: Formante 1
formant2_txt = data(:, 5);   % Quinta columna: Formante 2
formant3_txt = data(:, 6);   % Sexta columna: Formante 3

% Normalización de la intensidad del archivo TXT para que tenga un rango comparable
%intensity_txt = (intensity_txt - min(intensity_txt)) / (max(intensity_txt) - min(intensity_txt)) * 50 + 20;  % Escalar a un rango similar [20, 70]
% Convertir la intensidad a dB
%intensity_txt = 20 * log10(intensity_txt + eps);  % Agregar eps para evitar log(0)


%% Parte 3: Formantes del archivo WAV
% Preénfasis
pre_emphasis = [1, -0.97];
audio_preemphasized = filter(pre_emphasis, 1, snd);

% Reducir la frecuencia de muestreo
Fs_new = 10000;
audio_resampled = resample(audio_preemphasized, Fs_new, Fs);

% Parámetros de frames
frame_duration = 0.025; % 25 ms
frame_shift = 0.005;    % 5 ms

frame_length = round(frame_duration * Fs_new);
frame_step = round(frame_shift * Fs_new);

% Número total de frames
num_frames = floor((length(audio_resampled) - frame_length) / frame_step) + 1;

% Preparar matrices para almacenar los resultados
N = 4; % Número de formantes
formants = zeros(num_frames, N);
time_vector = zeros(num_frames, 1);

% Orden del modelo AR
p = round(2 + Fs_new / 1000);

for i = 1:num_frames
    % Índices del frame actual
    start_index = (i-1)*frame_step + 1;
    end_index = start_index + frame_length - 1;

    % Extraer el frame
    frame = audio_resampled(start_index:end_index);

    % Comprobar si el frame contiene valores no válidos o es silencioso
    if any(isnan(frame)) || any(isinf(frame)) || all(frame == 0) || sum(frame.^2) < 1e-6
        % Asignar ceros a los formantes y continuar
        formants(i, :) = 0;
        time_vector(i) = ((start_index + end_index) / 2) / Fs;
        continue;
    end

    % Aplicar ventana de Hamming
    frame = frame .* hamming(length(frame));

    % Coeficientes LPC usando el método autocorrelativo
    a = lpc(frame, p);


    % Comprobar si 'a' contiene valores no válidos
    if any(isnan(a)) || any(isinf(a))
        % Asignar ceros a los formantes y continuar
        formants(i, :) = 0;
        time_vector(i) = ((start_index + end_index) / 2) / Fs_new;
        continue;
    end

    % Raíces del polinomio AR
    roots_ar = roots(a);

    % Solo raíces en el semicírculo superior
    roots_pos = roots_ar(imag(roots_ar) >= 0);

    % Frecuencias de los formantes
    angulos = atan2(imag(roots_pos), real(roots_pos));
    frecuencias_formantes = angulos * (Fs_new / (2 * pi));

    % Filtrar frecuencias válidas
    indices_validos = (frecuencias_formantes > 200) & (frecuencias_formantes <5000);
    frecuencias_formantes = frecuencias_formantes(indices_validos);

    % Ordenar y seleccionar los primeros N formantes
    frecuencias_formantes = sort(frecuencias_formantes);
    num_formantes = min(N, length(frecuencias_formantes));
    formants(i, 1:num_formantes) = frecuencias_formantes(1:num_formantes);

    % Tiempo medio del frame
    time_vector(i) = ((start_index + end_index)/2)/Fs_new;
end

%% Espectrogramas del archivo TXT

% Cargar el archivo TXT
[file_txt, path_txt] = uigetfile('*.txt', 'Selecciona spectrogram.txt');
fullPath_txt = fullfile(path_txt, file_txt);
data = readmatrix(fullPath_txt, 'Delimiter', '/');

time_txt = data(:, 1);       % Primera columna: Tiempo
frequency_txt = data(:, 2);  % Segunda columna: Frecuencia
power_txt = data(:, 3);      % Tercera columna: Potencia

% Obtener los valores únicos de tiempo y frecuencia
unique_times = unique(time_txt);
unique_frequencies = unique(frequency_txt);

% Reorganizar la matriz de potencia para coincidir con las dimensiones del espectrograma
Z_txt = reshape(power_txt, length(unique_frequencies), length(unique_times));

%% 1. Espectrograma 2D del archivo WAV
window_length = round(0.1 * Fs);  % Ventana de 100 ms
time_step = 0.002;  % Paso de tiempo de 2 ms en segundos
overlap = round(window_length - (time_step * Fs));  % Solapamientonfft = 1024;  % Tamaño de la FFT
max_freq = 8000;  % Frecuencia máxima de visualización
nfft = 0.1*Fs;  % Resolución de frecuencia de 20 Hz

gaussian_window = gausswin(window_length);
[S, F, T, P] = spectrogram(snd, window_length, overlap, nfft, Fs, 'yaxis');
F = F(F <= max_freq);
P = P(1:length(F), :);
x= P*32767;
P_dB = 10 * log10(x);  % Normalizar la potencia a dBFS



%% oscilogram
[file_txt, path_txt] = uigetfile('*.txt', 'Selecciona oscilogram.txt');
fullPath_txt = fullfile(path_txt, file_txt);
data = readmatrix(fullPath_txt, 'Delimiter', '/');

time_oscilogram = data(:, 1);       % Primera columna: Tiempo
amplitud_oscilogram_txt = data(:, 2);  % Segunda columna: Frecuencia

%% Comparar las dos gráficas de Pitch e Intensidad lado a lado

figure;

% Subgráfica 1: Datos del archivo WAV
subplot(5, 2, 1);  % Crear la primera subgráfica
yyaxis left;  % Eje izquierdo para el pitch
plot(t, pitch_values_interp_wav, 'b', 'LineWidth', 2);
ylabel('Pitch [Hz]');
ylim([0 700]);
xlabel('Tiempo [s]');
yyaxis right;  % Eje derecho para la intensidad
plot(t, intensity_dB_wav, 'r', 'LineWidth', 2);
ylabel('Intensidad [dB]');
ylim([20 120]);
title('Pitch e Intensidad (Archivo WAV)');

% Subgráfica 2: Datos del archivo TXT
subplot(5, 2, 2);  % Crear la segunda subgráfica
yyaxis left;  % Eje izquierdo para el pitch
plot(time_txtp, pitch_txt, 'b', 'LineWidth', 2);
ylabel('Pitch [Hz]');
ylim([0 700]);
xlabel('Tiempo [s]');
yyaxis right;  % Eje derecho para la intensidad
plot(time_txtp, intensity_txt, 'r', 'LineWidth', 2);
ylabel('Intensidad [dB]');
ylim([20 120]);  % Ajustar el rango de la intensidad
title('Pitch e Intensidad (Archivo TXT)');


% Subgráfica 3: Formantes del archivo WAV
subplot(5, 2, 3);
plot(time_vector, formants(:, 1), 'r', 'LineWidth', 1.5); hold on;
plot(time_vector, formants(:, 2), 'g', 'LineWidth', 1.5);
plot(time_vector, formants(:, 3), 'b', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('Frecuencia [Hz]');
title('Formantes de la Señal de Voz (Archivo WAV)');
legend('Formante 1', 'Formante 2', 'Formante 3');
grid on;

%% Parte 4: Graficar los Formantes del archivo TXT



% Subgráfica 4: Formantes del archivo TXT
subplot(5, 2, 4);
plot(time_txtp, formant1_txt, 'r', 'LineWidth', 1.5); hold on;
plot(time_txtp, formant2_txt, 'g', 'LineWidth', 1.5);
plot(time_txtp, formant3_txt, 'b', 'LineWidth', 1.5);
xlabel('Tiempo [s]');
ylabel('Frecuencia [Hz]');
title('Formantes de la Señal de Voz (Archivo TXT)');
legend('Formante 1', 'Formante 2', 'Formante 3');
grid on;

% Subgráfica 1: Espectrograma 2D del archivo WAV
subplot(5, 2, 5);
imagesc(T, F, P_dB);
axis xy;
colormap('hot');  % Cambiar colormap
colorbar;
title('Espectrograma 2D (WAV)');
xlabel('Tiempo [s]');
ylabel('Frecuencia [Hz]');
ylim([0 8000]);


%% 3. Espectrograma 2D del archivo TXT
subplot(5, 2, 6);
imagesc(unique_times, unique_frequencies, Z_txt);
set(gca, 'YDir', 'normal');  % Invertir el eje Y
colormap('hot');  % Cambiar colormap
colorbar;
title('Espectrograma 2D (TXT)');
xlabel('Tiempo [s]');
ylabel('Frecuencia [Hz]');
ylim([0 8000]);


%% 2. Espectrograma 3D del archivo WAV
subplot(5, 2, 7);
[X, Y] = meshgrid(T, F);  % Crear un grid con tiempo y frecuencia
Z = P_dB;  % Convertir la potencia a dB
surf(X, Y, Z, 'EdgeColor', 'none');  % Crear superficie en 3D
axis tight;
view(3);  % Vista en 3D
colormap('hot');  % Cambiar colormap
colorbar;
title('Espectrograma 3D (WAV)');
xlabel('Tiempo [s]');
ylabel('Frecuencia [Hz]');
zlabel('Intensidad [dB]');
ylim([0 8000]);




%% 4. Espectrograma 3D del archivo TXT
subplot(5, 2, 8);
[X_txt, Y_txt] = meshgrid(unique_times, unique_frequencies);  % Crear el grid de tiempo y frecuencia
surf(X_txt, Y_txt, Z_txt, 'EdgeColor', 'none');  % Crear superficie en 3D
view(3);  % Vista en 3D
colormap('hot');  % Cambiar colormap
colorbar;
title('Espectrograma 3D (TXT)');
xlabel('Tiempo [s]');
ylabel('Frecuencia [Hz]');
zlabel('Potencia [dB]');
ylim([0 8000]);


%% 5. Oscilogram del archivo wav
subplot(5, 2, 9);  % Crear la primera subgráfica
yyaxis left;  % Eje izquierdo para el pitch
plot(t, amplitud_normalizada, 'b', 'LineWidth', 2);
ylabel('Ampllitud');
ylim([-1 1]);
xlabel('Tiempo [s]');
title('Oscilograma (Archivo wav)');

%% 5. Oscilogram del archivo txt
subplot(5, 2, 10);  % Crear la primera subgráfica
yyaxis left;  % Eje izquierdo para el pitch
plot(time_oscilogram, amplitud_oscilogram_txt, 'b', 'LineWidth', 2);
ylabel('Ampllitud');
ylim([-1 1]);
xlabel('Tiempo [s]');
title('Oscilograma (Archivo txt)');

amplitud_interpolada = interp1(t, amplitud_normalizada, time_oscilogram, 'linear', 'extrap');

pitch_txt(pitch_txt == 0) = nan;
valid_indices = ~isnan(pitch_txt);
time_txtp_valid = time_txtp(valid_indices);
pitch_txt_valid = pitch_txt(valid_indices);



pitch_interpolado = interp1(t, pitch_values_interp_wav, time_txtp_valid, 'linear','extrap');
intensidad_interpolada = interp1(t, intensity_dB_wav, time_txtp, 'linear', 'extrap');



differences = (pitch_interpolado - pitch_txt_valid).^2;

RMSE_adj_F = sqrt(mean(differences));


% MedAE_F = mean(abs((pitch_interpolado - pitch_txt_valid)));
disp(['RMSE (frecuencia): ', num2str(RMSE_adj_F)]);



formant1_txt(formant1_txt == 0) = nan;
formant2_txt(formant2_txt == 0) = nan;
formant3_txt(formant3_txt == 0) = nan;

valid_indices1 = ~isnan(formant1_txt);
valid_indices2 = ~isnan(formant2_txt);
valid_indices3 = ~isnan(formant3_txt);

time1_txtp_valid = time_txtp(valid_indices1);
formant1_txt_valid = formant1_txt(valid_indices1);

time2_txtp_valid = time_txtp(valid_indices2);
formant2_txt_valid = formant2_txt(valid_indices2);

time3_txtp_valid = time_txtp(valid_indices3);
formant3_txt_valid = formant3_txt(valid_indices3);

formant1_interpolado = interp1(time_vector, formants(:, 1), time1_txtp_valid, 'linear','extrap');
formant2_interpolado = interp1(time_vector, formants(:, 2), time2_txtp_valid, 'linear','extrap');
formant3_interpolado = interp1(time_vector, formants(:, 3), time3_txtp_valid, 'linear','extrap');

differences = (formant1_interpolado - formant1_txt_valid).^2;
RMSE_adj_F1 = sqrt(mean(differences));
disp(['RMSE (formante 1): ', num2str(RMSE_adj_F1)]);

differences = (formant2_interpolado - formant2_txt_valid).^2;
RMSE_adj_F2 = sqrt(mean(differences));
disp(['RMSE (formante 2): ', num2str(RMSE_adj_F2)]);

differences = (formant3_interpolado - formant3_txt_valid).^2;
RMSE_adj_F3 = sqrt(mean(differences));
disp(['RMSE (formante 3): ', num2str(RMSE_adj_F3)]);


% Suponiendo que tu vector se llama "data"
n = length(formant1_interpolado);             % Longitud del vector
remove_count = round(0.05 * n); % Número de elementos a eliminar (5% al inicio y al final)
% Crear un nuevo vector sin el 5% de elementos al inicio y al final
formant1_interpolado = formant1_interpolado(remove_count+1 : n-remove_count);

n = length(formant2_interpolado);             % Longitud del vector
remove_count = round(0.05 * n); % Número de elementos a eliminar (5% al inicio y al final)
% Crear un nuevo vector sin el 5% de elementos al inicio y al final
formant2_interpolado = formant2_interpolado(remove_count+1 : n-remove_count);

n = length(formant3_interpolado);             % Longitud del vector
remove_count = round(0.05 * n); % Número de elementos a eliminar (5% al inicio y al final)
% Crear un nuevo vector sin el 5% de elementos al inicio y al final
formant3_interpolado = formant3_interpolado(remove_count+1 : n-remove_count);

time_vector = time_vector(remove_count+1 : n-remove_count);




% Suponiendo que tu vector se llama "data"
n = length(formant1_txt_valid);             % Longitud del vector
remove_count = round(0.05 * n); % Número de elementos a eliminar (5% al inicio y al final)
% Crear un nuevo vector sin el 5% de elementos al inicio y al final
formant1_txt_valid = formant1_txt_valid(remove_count+1 : n-remove_count);
time1_txtp_valid = time1_txtp_valid(remove_count+1 : n-remove_count);


n = length(formant2_txt_valid);             % Longitud del vector
remove_count = round(0.05 * n); % Número de elementos a eliminar (5% al inicio y al final)
% Crear un nuevo vector sin el 5% de elementos al inicio y al final
formant2_txt_valid = formant2_txt_valid(remove_count+1 : n-remove_count);
time2_txtp_valid = time2_txtp_valid(remove_count+1 : n-remove_count);


n = length(formant3_txt_valid);             % Longitud del vector
remove_count = round(0.05 * n); % Número de elementos a eliminar (5% al inicio y al final)
% Crear un nuevo vector sin el 5% de elementos al inicio y al final
formant3_txt_valid = formant3_txt_valid(remove_count+1 : n-remove_count);
time3_txtp_valid = time3_txtp_valid(remove_count+1 : n-remove_count);







% disp(['MedAE (frecuencia): ', num2str(MedAE_F)]);
% 
differences = (intensidad_interpolada - intensity_txt).^2;
RMSE_adj_I = sqrt(mean(differences));
% MedAE_I = mean(abs((intensidad_interpolada - intensity_txt)));
disp(['RMSE (intensidad): ', num2str(RMSE_adj_I)]);
% disp(['MedAE (intensidad): ', num2str(MedAE_I)]);
% 
differences = (amplitud_interpolada - amplitud_oscilogram_txt).^2;
RMSE_adj_O = sqrt(mean(differences));
% MedAE_O = mean(abs((amplitud_interpolada - amplitud_oscilogram_txt)));
% 
disp(['RMSE (Oscilograma): ', num2str(RMSE_adj_O)]);
% disp(['MedAE (Oscilograma): ', num2str(MedAE_O)]);
% 
% error_cuadratico_pitch = mean((pitch_interpolado - pitch_txt_valid).^2);error_cuadratico_intensidad = mean((intensidad_interpolada - intensity_txt).^2);
% error_cuadratico = mean((amplitud_interpolada - amplitud_oscilogram_txt).^2);
% disp(['Error Cuadrático Medio (Oscilogram): ', num2str(error_cuadratico)]);

% disp(['Error Cuadrático Medio (Pitch): ', num2str(error_cuadratico_pitch)]);
% disp(['Error Cuadrático Medio (Intensidad): ', num2str(error_cuadratico_intensidad)]);

%  %
figure;

%Subgráfica 1: Datos del archivo WAV
subplot(1, 2, 1);  % Crear la primera subgráfica
yyaxis left;  % Eje izquierdo para el pitch
plot(time1_txtp_valid, formant1_interpolado, 'b', 'LineWidth', 2);hold on;
plot(time2_txtp_valid, formant2_interpolado, 'r', 'LineWidth', 2);
plot(time3_txtp_valid, formant3_interpolado, 'g', 'LineWidth', 2);
ylabel('Pitch [Hz]');
ylim([0 4000]);
xlabel('Tiempo [s]');
title('Pitch e Intensidad (Archivo WAV)');

%Subgráfica 2: Datos del archivo TXT
subplot(1, 2, 2);  % Crear la segunda subgráfica
yyaxis left;  % Eje izquierdo para el pitch
plot(time1_txtp_valid, formant1_txt_valid, 'b', 'LineWidth', 2);hold on;
plot(time2_txtp_valid, formant2_txt_valid, 'r', 'LineWidth', 2);
plot(time3_txtp_valid, formant3_txt_valid, 'g', 'LineWidth', 2);
ylabel('Pitch [Hz]');
ylim([0 4000]);
xlabel('Tiempo [s]');
title('Pitch e Intensidad (Archivo TXT)');

