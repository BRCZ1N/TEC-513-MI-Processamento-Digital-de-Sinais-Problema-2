%% TEC513 - PDS - UEFS | Problema 02 | Análise do ADC do ATmega2560
clc;
clear;
close all;
pkg load instrument-control
if ~isempty(mfilename), cd(fileparts(mfilename('fullpath'))); end   % roda na pasta do script

function v = triangular(t, f, fase, vmin, vmax)       % triangular ideal
  fr = mod(f*t + fase, 1);
  v  = vmin + (vmax - vmin)*(2*fr.*(fr < 0.5) + (2 - 2*fr).*(fr >= 0.5));
end

function [fb, pb] = ajustar(v, t, fs_, ps_, vmin, vmax) % triangular que melhor se ajusta
  e = inf;
  for f = fs_
    for p = ps_
      ei = sum((v - triangular(t, f, p, vmin, vmax)).^2);
      if ei < e, e = ei; fb = f; pb = p; end
    end
  end
end

%% Parâmetros
MODO      = 'serial';             % 'serial' (Arduino) ou 'arquivo' (csv salvo)
%PORTA     = 'COM5';              % WINDOWS - Abre a porta serial que esta no microcontrolador
PORTA     = "/dev/ttyUSB0";       % LINUX - Abre a porta serial que esta no microcontrolador
ARQUIVO   = 'captura.csv';
FS        = 16e6/(8*(38+1));      % fs nominal: OCR1A = 38 -> 51 282 SPS
F_GERADOR = 100;                  % frequência do gerador (Hz)
VREF      = 5.0;                  % AVcc medido
VIN_MIN   = 0.0;                  % extremos da triangular (osciloscópio)
VIN_MAX = 5.0;
Q         = VREF/256;             % degrau de quantização (8 bits)

%% 1) Recepção dos dados (requisitos 3 e 4)
if strcmp(MODO, 'serial')
  s1 = serial(PORTA);                 % abre a porta serial do microcontrolador
  set(s1,'baudrate',2000000);         % velocidade de transmissão 2 Mbps (igual ao .ino)
  set(s1,'bytesize',8);               % 8 bits de dados
  set(s1,'parity','n');               % sem paridade
  set(s1,'stopbits',1);               % 1 bit de parada
  set(s1,'timeout',10);               % 1 s sem dados = fim da transmissão (décimos de s)
  pause(2);                           % Arduino reinicia, captura e envia
  raw = [];                           % dados crus recebidos pela USB
  while true
    data = fread(s1, 4096);           % lê em blocos até parar de chegar dado
    if isempty(data), break; end
    raw = [raw; data(:)];
  end
  fclose(s1);                         % fecha a porta de comunicação
  c = str2double(strsplit(strtrim(char(raw')), "\n"))';   % texto -> números
  c = c(~isnan(c));                   % remove a linha vazia inicial
  printf('Amostras recebidas: %d\n', numel(c));
  csvwrite(ARQUIVO, c);               % guarda a captura
else
  c = csvread(ARQUIVO);
end
c = c(6:end);                                     % descarta o início (1ª conversão = 25 ciclos)
N = numel(c);
t = (0:N-1)'/FS;
v = (c + 0.5)*Q;                                  % reconstrução (centro do degrau)

%% 2) Triangular ideal simulada, alinhada ao sinal real (requisito 5)
[f, fase] = ajustar(v, t, F_GERADOR*(0.97:0.0005:1.03), 0:0.005:1, VIN_MIN, VIN_MAX);   % busca grossa
[f, fase] = ajustar(v, t, f + (-0.05:0.005:0.05), fase + (-0.005:0.0002:0.005), VIN_MIN, VIN_MAX); % fina
v_id = triangular(t, f, fase, VIN_MIN, VIN_MAX);  % entrada analógica ideal
c_id = min(max(floor(v_id/Q), 0), 255);           % saída de um ADC ideal de 8 bits

%% 3) Osciloscópio: real x ideal no tempo (requisito 4)
figure(1);
subplot(2,1,1);
stairs(t*1e3, c*Q, 'b');
hold on;
plot(t*1e3, v_id, 'r--');
grid on;
xlabel('Tempo (ms)');
ylabel('Tensão (V)');
xlim([0 30]);
title('Sinal capturado x triangular ideal');
legend('Real', 'Ideal');
subplot(2,1,2);                                   % zoom: entrada, amostras e retenção
j = find(c_id >= 128, 1) + (-8:8);
plot(t(j)*1e3, v_id(j), 'k');
hold on;
plot(t(j)*1e3, c(j)*Q, 'bo', 'markerfacecolor', 'b');   % amostras (antes era o stem)
stairs(t(j)*1e3, c(j)*Q, 'r');
stairs(t(j)*1e3, c_id(j)*Q, 'g--');
ylim([min(c(j))-2, max(c(j))+2]*Q);                     % nova: aproxima nos degraus
grid on; xlabel('Tempo (ms)');
ylabel('Tensão (V)');
title('Zoom: entrada, amostras e retenção');
legend('Entrada', 'Amostras', 'Retenção (real)', 'Ideal quantizado', 'location', 'northwest');
print(1, 'fig1_osciloscopio.png', '-dpng', '-r150');

%% 4) Comparação na frequência (requisito 4)
w  = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));
Xr = abs(fft((v    - mean(v)).*w));
Xi = abs(fft((v_id - mean(v_id)).*w));
fq = (0:N-1)'*FS/N;  k = 1:floor(N/2);
figure(2);
plot(fq(k), 20*log10(Xr(k)/max(Xi)), 'b');
hold on;
plot(fq(k), 20*log10(Xi(k)/max(Xi)), 'r--');
grid on;
xlim([0 4000]);
ylim([-100 5]);
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
title('Espectro: real x ideal');
legend('Real', 'Ideal');
print(2, 'fig2_fft.png', '-dpng', '-r150');

%% 5) Erros estáticos pelo histograma (requisito 5)
r = floor(2*(f*t + fase));                        % número da rampa de cada amostra
completa = r > r(1) & r < r(end);                 % só rampas completas
h = accumarray(c(completa) + 1, 1, [256 1]);      % h(k+1) = amostras do código k

T = VIN_MIN + (VIN_MAX - VIN_MIN)*cumsum(h(1:255))/sum(h);   % transição k-1 -> k (V)
Qr       = (T(255) - T(1))/254;                   % degrau médio real
offset   = (T(1) - Q)/Q;                          % LSB (transição ideal 0->1 em Q)
ganho    = (Qr/Q - 1)*100;                        % %
dnl      = diff(T)/Qr - 1;                        % códigos 1..254 (LSB)
inl      = (T - (T(1) + (0:254)'*Qr))/Qr;         % transições 1..255 (LSB)
ausentes = find(h(2:255) == 0)';                  % códigos 1..254 sem amostra
sobe     = mod(r, 2) == 0 & completa;             % amostras nas subidas
quedas   = sum(diff(c(sobe)) < -1 & diff(r(sobe)) == 0);  % código caiu > 1 LSB subindo
monotonico = isempty(ausentes) && quedas == 0;

figure(3);
bar(0:255, h, 1);
grid on;
xlim([-1 256]);
xlabel('Código');
ylabel('Ocorrências');
title(sprintf('Histograma de códigos (%d ausentes)', numel(ausentes)));
print(3, 'fig3_histograma.png', '-dpng', '-r150');

figure(4);
subplot(2,1,1);
stem(1:254, dnl, 'marker', 'none');
grid on;
xlim([0 255]);
ylabel('DNL (LSB)');
title(sprintf('DNL: %.2f a %.2f LSB', min(dnl), max(dnl)));
subplot(2,1,2); plot(1:255, inl, 'r');
grid on;
xlim([0 255]);
xlabel('Código');
ylabel('INL (LSB)');
title(sprintf('INL: %.2f a %.2f LSB', min(inl), max(inl)));
print(4, 'fig4_dnl_inl.png', '-dpng', '-r150');

figure(5);                                        % curva de transferência: início e fim
subplot(1,2,1);
stairs([0; T(1:8)], (0:8)', 'b');
hold on;
stairs((0:8)'*Q, (0:8)', 'r--');
grid on;
xlabel('Entrada (V)');
ylabel('Código');
title(sprintf('Offset: %.2f LSB', offset));
legend('Real', 'Ideal', 'location', 'northwest');
subplot(1,2,2);
stairs(T(247:255), (247:255)', 'b');
hold on;
stairs((247:255)'*Q, (247:255)', 'r--');
grid on;
xlabel('Entrada (V)');
title(sprintf('Ganho: %.2f%%', ganho));
print(5, 'fig5_offset_ganho.png', '-dpng', '-r150');

%% 6) Resultados (requisitos 5 e 6)
printf('\n--- Erros estáticos (Q = %.2f mV) ---\n', Q*1e3);
printf('Offset: %.2f LSB | Ganho: %.2f%%\n', offset, ganho);
printf('DNL: %.2f a %.2f LSB | INL: %.2f a %.2f LSB\n', min(dnl), max(dnl), min(inl), max(inl));
printf('Códigos ausentes: %d | Monotônico: %d (1 = sim)\n', numel(ausentes), monotonico);
printf('\n--- Características dinâmicas ---\n');
printf('Faixa dinâmica (8 bits): %.1f dB\n', 20*log10(256));
printf('Tempo de conversão (datasheet, auto-trigger): 13,5 us -> taxa máxima 74,1 kSPS\n');
printf('Taxa de amostragem medida: %.0f SPS (nominal %.0f)\n', FS*F_GERADOR/f, FS);
