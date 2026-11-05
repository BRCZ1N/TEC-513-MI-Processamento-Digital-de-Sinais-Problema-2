%% =========================================================================
%  TEC513 - MI de Processamento de Sinais - UEFS
%  Problema 02 - Análise do conversor A/D do ATmega2560 (8 bits)
%
%  Etapas:
%    1) Recepção das amostras (serial ou arquivo .csv)
%    2) Triangular ideal simulada, alinhada ao sinal capturado
%    3) Comparação no tempo (estilo osciloscópio)
%    4) Comparação na frequência (FFT)
%    5) Erros estáticos: offset, ganho, DNL, INL, códigos ausentes, monotonicidade
%    6) Resumo dos resultados e características dinâmicas
%  =========================================================================
clc;
clear;
close all;
pkg load instrument-control

% Faz o Octave trabalhar na pasta do script (onde ficam o .csv e as figuras)
if ~isempty(mfilename)
  cd(fileparts(mfilename('fullpath')));
end

% Gera uma onda triangular ideal entre vmin e vmax
function v = triangular(t, f, fase, vmin, vmax)
  posicao = mod(f*t + fase, 1);              % posição dentro do período (0 a 1)
  subida  = posicao < 0.5;                   % primeira metade: subida
  forma   = 2*posicao.*subida + (2 - 2*posicao).*(~subida);   % vai de 0 a 1 e volta
  v = vmin + (vmax - vmin)*forma;
end

% Procura a frequência e a fase da triangular que mais se aproximam do sinal
% medido (menor erro quadrático entre as listas de valores testados)
function [fMelhor, faseMelhor] = ajustarTriangular(v, t, listaF, listaFase, vmin, vmax)
  menorErro = inf;
  for f = listaF
    for fase = listaFase
      erro = sum((v - triangular(t, f, fase, vmin, vmax)).^2);
      if erro < menorErro
        menorErro  = erro;
        fMelhor    = f;
        faseMelhor = fase;
      end
    end
  end
end


%% ===================== PARÂMETROS (ajustar antes de rodar) ===============
MODO     = 'arquivo';           % 'serial' = captura com o Arduino
                                % 'arquivo' = usa um .csv já salvo
%PORTA   = 'COM5';              % Windows
PORTA    = "/dev/ttyUSB0";      % Linux
ARQUIVO  = 'captura_1_2.5V_50Ohms.csv';% onde a captura é salva (ou lida)

% Medidas feitas no laboratório
VREF     = 5.0;                 % AVcc medido com multímetro (V)
VIN_MIN  = 0.0;                 % mínimo da triangular no osciloscópio (V)
VIN_MAX  = 5.0;                 % máximo da triangular no osciloscópio (V)
F_GERADOR = 100;                % frequência ajustada no gerador (Hz)

% Valores do projeto
OCR1A    = 38;
FS       = 16e6/(8*(OCR1A + 1));% taxa de amostragem: 51 282 SPS
NBITS    = 8;
NCODIGOS = 2^NBITS;             % 256 códigos (0 a 255)
Q        = VREF/NCODIGOS;       % degrau de quantização (V)


%% ===================== 1) RECEPÇÃO DOS DADOS (requisitos 3 e 4) ==========
if strcmp(MODO, 'serial')
  % Abrir a porta reinicia o Arduino: ele captura as amostras e as envia
  s1 = serial(PORTA);
  set(s1, 'baudrate', 2000000);     % igual ao Serial.begin() do .ino
  set(s1, 'bytesize', 8);           % 8 bits de dados
  set(s1, 'parity',   'n');         % sem paridade
  set(s1, 'stopbits', 1);           % 1 bit de parada
  set(s1, 'timeout',  10);          % 1 s sem dados = fim da transmissão (décimos de s)
  pause(2);                         % tempo para o Arduino reiniciar

  % Lê em blocos até parar de chegar dado
  bytesRecebidos = [];
  while true
    bloco = fread(s1, 4096);
    if isempty(bloco)
      break;
    end
    bytesRecebidos = [bytesRecebidos; bloco(:)];
  end
  fclose(s1);

  % Converte o texto recebido ("123\r\n124\r\n...") em números
  texto   = char(bytesRecebidos');
  linhas  = strsplit(strtrim(texto), "\n");
  codigos = str2double(linhas)';
  codigos = codigos(~isnan(codigos));       % remove linhas vazias
  printf('Amostras recebidas: %d\n', numel(codigos));

  csvwrite(ARQUIVO, codigos);               % guarda a captura para uso posterior
else
  codigos = csvread(ARQUIVO);
end

N       = numel(codigos);
t       = (0:N-1)'/FS;                      % instante de cada amostra (s)
tMs     = t*1e3;                            % o mesmo, em ms (para os gráficos)
vReal   = (codigos + 0.5)*Q;                % tensão reconstruída (centro do degrau)


%% ===================== 2) TRIANGULAR IDEAL SIMULADA (requisito 5) =========
% O gerador nunca está exatamente em 100 Hz e a captura começa numa fase
% qualquer: por isso a triangular ideal é ajustada ao sinal real.

% Busca grossa: frequência entre 97% e 103% do gerador, fase de 0 a 1
listaF    = F_GERADOR*(0.97:0.0005:1.03);
listaFase = 0:0.005:1;
[fEst, fase] = ajustarTriangular(vReal, t, listaF, listaFase, VIN_MIN, VIN_MAX);

% Busca fina em torno do resultado anterior
listaF    = fEst + (-0.05:0.005:0.05);
listaFase = fase  + (-0.005:0.0002:0.005);
[fEst, fase] = ajustarTriangular(vReal, t, listaF, listaFase, VIN_MIN, VIN_MAX);

vIdeal      = triangular(t, fEst, fase, VIN_MIN, VIN_MAX);    % entrada analógica ideal
codigoIdeal = min(max(floor(vIdeal/Q), 0), NCODIGOS - 1);    % saída de um ADC ideal


%% ===================== 3) COMPARAÇÃO NO TEMPO (requisito 4) ===============
figure(1);

% 3a) Visão geral: sinal capturado x triangular ideal
subplot(2,1,1);
stairs(tMs, codigos*Q, 'b');
hold on;
plot(tMs, vIdeal, 'r--');
grid on;
xlim([0 30]);
xlabel('Tempo (ms)');
ylabel('Tensão (V)');
title('Sinal capturado x triangular ideal');
legend('Real', 'Ideal');

% 3b) Zoom no meio de uma subida: entrada, amostras e retenção
subplot(2,1,2);
meio   = round((min(codigos) + max(codigos))/2);       % meio da faixa capturada
% primeira subida que passa pelo meio, longe das bordas do vetor
antes  = codigoIdeal(9:end-9);
depois = codigoIdeal(10:end-8);
inicio = find(antes < meio & depois >= meio, 1) + 9;
janela = inicio + (-8:8);                   % 17 amostras em torno desse ponto
plot(tMs(janela), vIdeal(janela), 'k');
hold on;
plot(tMs(janela), codigos(janela)*Q, 'bo', 'markerfacecolor', 'b');
stairs(tMs(janela), codigos(janela)*Q, 'r');
stairs(tMs(janela), codigoIdeal(janela)*Q, 'g--');
ylim([min(codigos(janela)) - 2, max(codigos(janela)) + 2]*Q);
grid on;
xlabel('Tempo (ms)');
ylabel('Tensão (V)');
title('Zoom: entrada, amostras e retenção');
legend('Entrada', 'Amostras', 'Retenção (real)', 'Ideal quantizado', 'location', 'northwest');

print(1, 'fig1_osciloscopio.png', '-dpng', '-r150');


%% ===================== 4) COMPARAÇÃO NA FREQUÊNCIA (requisito 4) ==========
janelaHann = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));   % reduz o vazamento espectral

xReal  = abs(fft((vReal  - mean(vReal )).*janelaHann));
xIdeal = abs(fft((vIdeal - mean(vIdeal)).*janelaHann));

freq    = (0:N-1)'*FS/N;                    % eixo de frequência (Hz)
metade  = 1:floor(N/2);                     % só frequências positivas
ref     = max(xIdeal);                      % 0 dB = fundamental ideal

figure(2);
stem(freq(metade), 20*log10(xReal(metade)/ref), 'b', 'marker', 'none', 'basevalue', -100);
hold on;
stem(freq(metade), 20*log10(xIdeal(metade)/ref), 'r', 'marker', 'none', 'basevalue', -100);
grid on;
xlim([0 4000]);
ylim([-100 5]);
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
title('Espectro: real x ideal');
legend('Real', 'Ideal');

print(2, 'fig2_fft.png', '-dpng', '-r150');


%% ===================== 5) ERROS ESTÁTICOS (requisito 5) ===================
% Método do histograma: com entrada triangular, todos os códigos deveriam
% aparecer o mesmo número de vezes. Código mais largo = mais ocorrências.

% Numera as rampas (subidas e descidas) e usa só as completas
rampa    = floor(2*(fEst*t + fase));
completa = rampa > rampa(1) & rampa < rampa(end);

% Histograma: histograma(k+1) = número de amostras do código k
histograma = accumarray(codigos(completa) + 1, 1, [NCODIGOS 1]);

% Tensão de cada transição (código k-1 -> k), para k = 1..255
vTransicao = VIN_MIN + (VIN_MAX - VIN_MIN)*cumsum(histograma(1:NCODIGOS-1))/sum(histograma);

% Offset: primeira transição real x ideal (ideal em 1 Q)
offsetLsb = (vTransicao(1) - Q)/Q;

% Ganho: degrau médio real x degrau ideal
qReal    = (vTransicao(end) - vTransicao(1))/(NCODIGOS - 2);
ganhoPct = (qReal/Q - 1)*100;

% DNL: largura de cada degrau comparada ao degrau médio (códigos 1..254)
dnl = diff(vTransicao)/qReal - 1;

% INL: distância de cada transição até a reta ideal (transições 1..255)
reta = vTransicao(1) + (0:NCODIGOS-2)'*qReal;
inl  = (vTransicao - reta)/qReal;

% Códigos ausentes: códigos 1..254 que nunca apareceram
ausentes = find(histograma(2:NCODIGOS-1) == 0)';

% Monotonicidade: nas subidas, o código não pode cair mais que 1 LSB
% (quedas de 1 LSB são ruído normal)
naSubida   = mod(rampa, 2) == 0 & completa;
codSubida  = codigos(naSubida);
mesmaRampa = diff(rampa(naSubida)) == 0;
quedas     = sum(diff(codSubida) < -1 & mesmaRampa);
monotonico = isempty(ausentes) && quedas == 0;

% 5a) Histograma
figure(3);
bar(0:NCODIGOS-1, histograma, 1);
grid on;
xlim([-1 NCODIGOS]);
xlabel('Código');
ylabel('Ocorrências');
title(sprintf('Histograma de códigos (%d ausentes)', numel(ausentes)));
print(3, 'fig3_histograma.png', '-dpng', '-r150');

% 5b) DNL e INL
figure(4);
subplot(2,1,1);
stem(1:NCODIGOS-2, dnl, 'marker', 'none');
grid on;
xlim([0 NCODIGOS-1]);
ylabel('DNL (LSB)');
title(sprintf('DNL: %.2f a %.2f LSB', min(dnl), max(dnl)));

subplot(2,1,2);
plot(1:NCODIGOS-1, inl, 'r');
grid on;
xlim([0 NCODIGOS-1]);
xlabel('Código');
ylabel('INL (LSB)');
title(sprintf('INL: %.2f a %.2f LSB', min(inl), max(inl)));
print(4, 'fig4_dnl_inl.png', '-dpng', '-r150');

% 5c) Curva de transferência: início (offset) e fim (ganho)
figure(5);
subplot(1,2,1);
stairs([0; vTransicao(1:8)], (0:8)', 'b');
hold on;
stairs((0:8)'*Q, (0:8)', 'r--');
grid on;
xlabel('Entrada (V)');
ylabel('Código');
title(sprintf('Offset: %.2f LSB', offsetLsb));
legend('Real', 'Ideal', 'location', 'northwest');

subplot(1,2,2);
stairs(vTransicao(247:255), (247:255)', 'b');
hold on;
stairs((247:255)'*Q, (247:255)', 'r--');
grid on;
xlabel('Entrada (V)');
title(sprintf('Ganho: %.2f%%', ganhoPct));
print(5, 'fig5_offset_ganho.png', '-dpng', '-r150');


%% ===================== 6) RESUMO (requisitos 5 e 6) =====================
fsMedida = FS*F_GERADOR/fEst;   % se o gerador está certo, mede a fs real

printf('\n--- Erros estáticos (Q = %.2f mV) ---\n', Q*1e3);
printf('Offset: %.2f LSB\n', offsetLsb);
printf('Ganho: %.2f%%\n', ganhoPct);
printf('DNL: %.2f a %.2f LSB\n', min(dnl), max(dnl));
printf('INL: %.2f a %.2f LSB\n', min(inl), max(inl));
printf('Códigos ausentes: %d\n', numel(ausentes));
printf('Monotônico: %d (1 = sim)\n', monotonico);

printf('\n--- Características dinâmicas ---\n');
printf('Faixa dinâmica (8 bits): %.1f dB\n', 20*log10(NCODIGOS));
printf('Tempo de conversão (datasheet, auto-trigger): 13,5 us\n');
printf('Taxa máxima (datasheet, auto-trigger): 74,1 kSPS\n');
printf('Taxa de amostragem medida: %.0f SPS (nominal %.0f)\n', fsMedida, FS);
