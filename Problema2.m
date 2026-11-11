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
%
%  Convenções (datasheet ATmega2560, seção 26.6.4, p. 277-278):
%    - Transição ideal k-1 -> k ocorre em (k - 0,5)·Q  (1ª transição em 0,5 LSB)
%    - Offset: 1ª transição real x ideal, em LSB
%    - Ganho: após remover o offset, última transição real x ideal
%      (ideal = 1,5 LSB abaixo do fundo de escala), em LSB
%    - INL e DNL pelo método dos extremos (reta entre 1ª e última transição)
%  =========================================================================
clc;
clear;
close all;
pkg load instrument-control      % comunicação serial
pkg load signal                  % função sawtooth (triangular)

% Faz o Octave trabalhar na pasta do script (onde ficam o .csv e as figuras)
if ~isempty(mfilename)
  cd(fileparts(mfilename('fullpath')));
end


%% ===================== FUNÇÃO AUXILIAR ===================================
% Gera uma onda triangular entre vmin e vmax.
%   t    : vetor de tempo (s)
%   f    : frequência (Hz)
%   fase : deslocamento em fração do período (0 a 1)
% sawtooth(x, 0.5) devolve uma triangular de -1 a +1 que SOBE na primeira
% metade do período e DESCE na segunda; aqui ela é levada para a faixa 0..1
% e depois para vmin..vmax.
function v = triangular(t, f, fase, vmin, vmax)
  forma = (sawtooth(2*pi*(f*t + fase), 0.5) + 1)/2;
  v = vmin + (vmax - vmin)*forma;
end


%% ===================== PARÂMETROS (ajustar antes de rodar) ===============
% --- Origem dos dados ---
MODO    = 'arquivo';                     % 'serial' = captura com o Arduino
                                         % 'arquivo' = usa um .csv já salvo
%PORTA  = 'COM5';                        % Windows
PORTA   = "/dev/ttyUSB0";                % Linux
ARQUIVO = 'captura_1_2.5V_50Ohms.csv';   % onde a captura é salva (ou lida)

% --- Medidas feitas no laboratório ---
VREF      = 5.0;     % AVcc medido com multímetro (V)
VIN_MIN   = 0.0;     % mínimo da triangular medido no osciloscópio (V)
VIN_MAX   = 5.0;     % máximo da triangular medido no osciloscópio (V)
F_ENTRADA = 100;     % frequência ajustada no gerador (Hz)

% --- Configuração do Arduino (deve bater com o .ino) ---
F_CPU          = 16e6;   % clock do ATmega2560 (Hz)
PRESCALER_TMR1 = 8;      % prescaler do Timer1 (CS11)
OCR1A          = 38;     % valor de comparação do Timer1 (modo CTC)
PRESCALER_ADC  = 16;     % prescaler do ADC (ADPS2) -> clock do ADC = 1 MHz
CICLOS_CONV    = 13.5;   % ciclos por conversão em auto-trigger (Tabela 26-1)
NBITS          = 8;      % resolução usada (ADLAR = 1, leitura só de ADCH)
BAUD           = 2e6;    % igual ao Serial.begin() do .ino

% --- Ajustes da análise ---
AMOSTRAS_DESCARTE = 10;   % amostras iniciais descartadas (1ª conversão leva
                          % 25 ciclos e o ADC ainda está estabilizando)
TOL_QUEDA_LSB     = 1;    % nas subidas, quedas de até 1 LSB são tratadas como
                          % ruído (critério de monotonicidade)
PASSO_FASE        = 0.01; % resolução da busca inicial de fase (fração do período)

% --- Ajustes dos gráficos ---
JANELA_TEMPO_MS = 30;     % trecho mostrado na visão geral (ms)
MEIA_JANELA     = 8;      % zoom: amostras antes/depois do ponto central
F_MAX_GRAFICO   = 1000;   % limite do eixo de frequência (Hz)
PISO_DB         = -100;   % piso dos gráficos em dB
N_HARMONICAS    = 10;     % harmônicas mostradas (1f a 10f)
N_BORDA         = 8;      % códigos mostrados no início/fim da curva de transferência
RESOLUCAO_FIG   = '-r150';% resolução das figuras salvas (dpi)

% --- Grandezas derivadas (calculadas, não digitadas) ---
FS        = F_CPU/(PRESCALER_TMR1*(1 + OCR1A));  % taxa de amostragem (SPS)
F_ADC     = F_CPU/PRESCALER_ADC;                 % clock do ADC (Hz)
T_CONV    = CICLOS_CONV/F_ADC;                   % tempo de conversão (s)
FS_MAX    = 1/T_CONV;                            % taxa máxima teórica (SPS)
NCODIGOS  = 2^NBITS;                             % 256 códigos (0 a 255)
COD_MAX   = NCODIGOS - 1;                        % código de fundo de escala
Q         = VREF/NCODIGOS;                       % degrau de quantização (V)

% Tensões ideais das transições k-1 -> k (k = 1..255), convenção do datasheet
vTransIdeal = ((1:COD_MAX)' - 0.5)*Q;


%% ===================== 1) RECEPÇÃO DOS DADOS (requisitos 3 e 4) ==========
% No modo 'serial', abrir a porta reinicia o Arduino; ele captura o buffer
% inteiro e só depois transmite (um número por linha, em texto). A leitura
% para quando a porta fica 1 s sem receber dados (timeout).
% No modo 'arquivo', reaproveita uma captura salva anteriormente.
if strcmp(MODO, 'serial')
  TIMEOUT_DS    = 10;     % timeout em décimos de segundo (10 = 1 s)
  PAUSA_RESET_S = 2;      % tempo para o Arduino reiniciar e capturar
  TAM_BLOCO     = 4096;   % bytes lidos por vez

  s1 = serial(PORTA);
  set(s1, 'baudrate', BAUD);
  set(s1, 'bytesize', 8);           % 8 bits de dados
  set(s1, 'parity',   'n');         % sem paridade
  set(s1, 'stopbits', 1);           % 1 bit de parada
  set(s1, 'timeout',  TIMEOUT_DS);
  srl_flush(s1);
  pause(PAUSA_RESET_S);

  % Lê em blocos até parar de chegar dado
  bytesRecebidos = [];
  while true
    bloco = fread(s1, TAM_BLOCO);
    if isempty(bloco)
      break;
    end
    bytesRecebidos = [bytesRecebidos; bloco(:)];
  end
  fclose(s1);

  % Converte o texto recebido ("123\r\n124\r\n...") em números
  texto   = char(bytesRecebidos');
  linhas  = strsplit(strtrim(texto), "\n");
  codigos = str2double(linhas);
  codigos = codigos(~isnan(codigos));       % remove linhas vazias
  printf('Amostras recebidas: %d\n', numel(codigos));

  csvwrite(ARQUIVO, codigos(:));            % guarda a captura para uso posterior
else
  codigos = csvread(ARQUIVO);
end

% Garante VETOR COLUNA: se o .csv vier como linha, as operações com t
% (coluna) gerariam uma matriz N x N sem aviso e o ajuste sairia errado.
codigos = codigos(:);
codigos = codigos(AMOSTRAS_DESCARTE+1:end);

N     = numel(codigos);
t     = (0:N-1)'/FS;          % instante de cada amostra (s)
tMs   = t*1e3;                % o mesmo, em ms (para os gráficos)
vReal = codigos*Q;            % tensão equivalente de cada código


%% ===================== 2) TRIANGULAR IDEAL SIMULADA (requisito 5) =========
% O gerador nunca está exatamente na frequência nominal e a captura começa
% numa fase qualquer. Por isso a frequência e a fase da triangular ideal são
% AJUSTADAS ao sinal real, minimizando o erro quadrático entre os dois.
%   Passo 1: busca grosseira da fase, com a frequência nominal.
%   Passo 2: fminsearch (Nelder-Mead) refina frequência e fase juntas.

erroQuad = @(p) sum((vReal - triangular(t, p(1), p(2), VIN_MIN, VIN_MAX)).^2);

listaFase = 0:PASSO_FASE:1;
erros     = arrayfun(@(fase) erroQuad([F_ENTRADA, fase]), listaFase);
[~, iMin] = min(erros);

opcoes = optimset('TolX', 1e-8, 'TolFun', 1e-10);
p      = fminsearch(erroQuad, [F_ENTRADA, listaFase(iMin)], opcoes);
fEst   = p(1);
fase   = p(2);

% Entrada analógica ideal e a saída de um ADC IDEAL para ela.
% round() coloca as transições em (k - 0,5)·Q, como no datasheet.
vIdeal      = triangular(t, fEst, fase, VIN_MIN, VIN_MAX);
codigoIdeal = min(max(round(vIdeal/Q), 0), COD_MAX);

printf('Frequência estimada: %.3f Hz\n', fEst);
printf('Fase: %.3f do período (%.1f graus)\n', mod(fase,1), mod(fase,1)*360);


%% ===================== 3) COMPARAÇÃO NO TEMPO (requisito 4) ===============
% 3a) Visão geral: sinal capturado (degraus) sobre a triangular ideal.
% 3b) Zoom no meio de uma subida: mostra a entrada contínua, as amostras,
%     a retenção (sample-and-hold) real e a saída do ADC ideal.
figure(1);

subplot(2,1,1);
stairs(tMs, vReal, 'b');
hold on;
plot(tMs, vIdeal, 'r--');
grid on;
xlim([0 JANELA_TEMPO_MS]);
xlabel('Tempo (ms)');
ylabel('Tensão (V)');
title('Sinal capturado x triangular ideal');
legend('Real', 'Ideal');

% Ponto central do zoom: primeira amostra em que a saída ideal SOBE e cruza
% o meio da faixa, com espaço para MEIA_JANELA amostras de cada lado.
meio    = round((min(codigos) + max(codigos))/2);
cruza   = find(codigoIdeal(1:end-1) < meio & codigoIdeal(2:end) >= meio) + 1;
cruza   = cruza(cruza > MEIA_JANELA & cruza <= N - MEIA_JANELA);
centro  = cruza(1);
janela  = centro + (-MEIA_JANELA:MEIA_JANELA);

subplot(2,1,2);
plot(tMs(janela), vIdeal(janela), 'k');
hold on;
plot(tMs(janela), vReal(janela), 'bo', 'markerfacecolor', 'b');
stairs(tMs(janela), vReal(janela), 'r');
stairs(tMs(janela), codigoIdeal(janela)*Q, 'g--');
ylim([min(codigos(janela)) - 2, max(codigos(janela)) + 2]*Q);
grid on;
xlabel('Tempo (ms)');
ylabel('Tensão (V)');
title('Zoom: entrada, amostras e retenção');
legend('Entrada', 'Amostras', 'Retenção (real)', 'Ideal quantizado', ...
       'location', 'northwest');

print(1, 'fig1_osciloscopio.png', '-dpng', RESOLUCAO_FIG);


%% ===================== 4) COMPARAÇÃO NA FREQUÊNCIA (requisito 4) ==========
% A FFT supõe que o trecho analisado se repete. Cortando o sinal em um número
% INTEIRO de períodos, cada harmônica cai exatamente em um bin e não há
% vazamento espectral. Com nPeriodos períodos no trecho, a harmônica h fica
% no bin h·nPeriodos (índice h·nPeriodos + 1, pois o Octave começa em 1).
% O nível médio (DC) é removido para não dominar o gráfico.
amostrasPorPeriodo = FS/fEst;
nPeriodos = floor(N/amostrasPorPeriodo);
M         = round(nPeriodos*amostrasPorPeriodo);

xReal  = abs(fft(vReal(1:M)  - mean(vReal(1:M))));
xIdeal = abs(fft(vIdeal(1:M) - mean(vIdeal(1:M))));

freq   = (0:M-1)'*FS/M;             % eixo de frequência (Hz)
metade = 1:floor(M/2);              % só frequências positivas
ref    = max(xIdeal);               % 0 dB = fundamental ideal
emDb   = @(x) 20*log10(x/ref + eps);% eps evita log(0)

% 4a) Espectro completo até F_MAX_GRAFICO
figure(2);
stem(freq(metade), emDb(xReal(metade)),  'b', 'marker', 'none', 'basevalue', PISO_DB);
hold on;
stem(freq(metade), emDb(xIdeal(metade)), 'r', 'marker', 'none', 'basevalue', PISO_DB);
grid on;
xlim([0 F_MAX_GRAFICO]);
ylim([PISO_DB 5]);
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
title(sprintf('Espectro: real x ideal (%d períodos inteiros)', nPeriodos));
legend('Real', 'Ideal');
print(2, 'fig2_fft.png', '-dpng', RESOLUCAO_FIG);

% 4b) Só as harmônicas (a triangular ideal tem apenas as ímpares, caindo
%     com 1/h^2; as pares ideais ficam no piso)
h    = 1:N_HARMONICAS;
bins = h*nPeriodos + 1;
figure(3);
stem(h*fEst, emDb(xReal(bins)), 'b', 'filled', 'basevalue', PISO_DB);
hold on;
plot(h*fEst, emDb(xIdeal(bins)), 'rx', 'markersize', 10, 'linewidth', 2);
grid on;
set(gca, 'xtick', h*fEst, ...
         'xticklabel', arrayfun(@(x) sprintf('%df', x), h, 'UniformOutput', false));
xlim([0 (N_HARMONICAS + 1)*fEst]);
ylim([PISO_DB 5]);
xlabel('Frequência (Hz)');
ylabel('Magnitude (dB)');
title('Harmônicas: real x ideal');
legend('Real', 'Ideal');
print(3, 'fig3_harmonicas.png', '-dpng', RESOLUCAO_FIG);


%% ===================== 5) ERROS ESTÁTICOS (requisito 5) ===================
% MÉTODO DO HISTOGRAMA: com entrada triangular, a tensão passa o mesmo tempo
% em cada faixa de tensão. Logo, um ADC ideal produz o mesmo número de
% ocorrências para todos os códigos; um código "mais largo" aparece mais.

% --- 5.1 Separação das rampas ---
% rampa = número de meio-períodos decorridos. Pares = subidas, ímpares =
% descidas (é assim que sawtooth(..., 0.5) foi definida). A primeira e a
% última rampa são incompletas e ficam de fora.
rampa    = floor(2*(fEst*t + fase));
completa = rampa > rampa(1) & rampa < rampa(end);

% --- 5.2 Histograma: histograma(k+1) = número de amostras do código k ---
histograma = accumarray(codigos(completa) + 1, 1, [NCODIGOS 1]);
amostrasPorCodigo = mean(histograma(2:COD_MAX));   % média, sem as bordas

% --- 5.3 Tensões de transição reais ---
% A fração acumulada de amostras até o código k-1 é a fração da faixa
% VIN_MIN..VIN_MAX abaixo da transição k-1 -> k. Por isso offset e ganho
% dependem de VIN_MIN e VIN_MAX terem sido bem medidos.
vTransicao = VIN_MIN + (VIN_MAX - VIN_MIN)*cumsum(histograma(1:COD_MAX))/sum(histograma);

% --- 5.4 Offset (LSB): 1ª transição real x ideal (0,5 LSB) ---
offsetV   = vTransicao(1) - vTransIdeal(1);
offsetLsb = offsetV/Q;

% --- 5.5 Ganho: última transição, já sem o offset, x ideal ---
% Em LSB (definição do datasheet) e em % (inclinação real x ideal).
ganhoLsb = ((vTransicao(end) - offsetV) - vTransIdeal(end))/Q;
qReal    = (vTransicao(end) - vTransicao(1))/(COD_MAX - 1);   % degrau médio real
ganhoPct = (qReal/Q - 1)*100;

% --- 5.6 DNL (códigos 1..254): largura de cada degrau / degrau médio - 1 ---
% Usar o degrau médio real (e não Q) remove offset e ganho da DNL.
dnl = diff(vTransicao)/qReal - 1;

% --- 5.7 INL (transições 1..255): distância até a reta entre os extremos ---
reta = vTransicao(1) + (0:COD_MAX-1)'*qReal;
inl  = (vTransicao - reta)/qReal;

% --- 5.8 Códigos ausentes (1..254; 0 e 255 são bordas e não contam) ---
ausentes = find(histograma(2:COD_MAX) == 0)';    % já em número de código

% --- 5.9 Monotonicidade ---
% Nas subidas, o código não pode cair. Quedas de até TOL_QUEDA_LSB são
% ruído normal. Só compara amostras vizinhas da MESMA rampa.
% (Código ausente NÃO torna o ADC não monotônico: é outro erro.)
naSubida   = mod(rampa, 2) == 0 & completa;
codSubida  = codigos(naSubida);
mesmaRampa = diff(rampa(naSubida)) == 0;
quedas     = sum(diff(codSubida) < -TOL_QUEDA_LSB & mesmaRampa);
monotonico = (quedas == 0);

% --- Figura: histograma ---
figure(4);
bar(0:COD_MAX, histograma, 1);
grid on;
xlim([-1 NCODIGOS]);
xlabel('Código');
ylabel('Ocorrências');
title(sprintf('Histograma de códigos (%d ausentes, média %.1f amostras/código)', ...
              numel(ausentes), amostrasPorCodigo));
print(4, 'fig4_histograma.png', '-dpng', RESOLUCAO_FIG);

% --- Figura: DNL e INL ---
figure(5);
subplot(2,1,1);
stem(1:COD_MAX-1, dnl, 'marker', 'none');
grid on;
xlim([0 COD_MAX]);
ylabel('DNL (LSB)');
title(sprintf('DNL: %.2f a %.2f LSB', min(dnl), max(dnl)));

subplot(2,1,2);
plot(1:COD_MAX, inl, 'r');
grid on;
xlim([0 COD_MAX]);
xlabel('Código');
ylabel('INL (LSB)');
title(sprintf('INL: %.2f a %.2f LSB', min(inl), max(inl)));
print(5, 'fig5_dnl_inl.png', '-dpng', RESOLUCAO_FIG);

% --- Figura: curva de transferência no início (offset) e no fim (ganho) ---
% Em stairs(x, y), cada degrau começa na tensão de transição x(k) e vai
% para o código y(k).
inicioCod = (0:N_BORDA)';
fimCod    = (COD_MAX-N_BORDA+1:COD_MAX)';

figure(6);
subplot(1,2,1);
stairs([0; vTransicao(inicioCod(2:end))],  inicioCod, 'b');
hold on;
stairs([0; vTransIdeal(inicioCod(2:end))], inicioCod, 'r--');
grid on;
xlabel('Entrada (V)');
ylabel('Código');
title(sprintf('Offset: %.2f LSB', offsetLsb));
legend('Real', 'Ideal', 'location', 'northwest');

subplot(1,2,2);
stairs(vTransicao(fimCod),  fimCod, 'b');
hold on;
stairs(vTransIdeal(fimCod), fimCod, 'r--');
grid on;
xlabel('Entrada (V)');
title(sprintf('Ganho: %.2f LSB (%.2f%%)', ganhoLsb, ganhoPct));
print(6, 'fig6_offset_ganho.png', '-dpng', RESOLUCAO_FIG);


%% ===================== 6) RESUMO (requisitos 5 e 6) =====================
% A taxa de amostragem "medida" compara a frequência nominal do gerador com
% a frequência que o ajuste encontrou nos dados: se o gerador está certo,
% qualquer diferença vem da taxa real do Arduino.
% O tempo de conversão NÃO é medido aqui: isso é feito no Arduino
% (ex.: alternar um pino na ISR do ADC e medir no osciloscópio).
fsMedida = FS*F_ENTRADA/fEst;

printf('\n--- Erros estáticos (Q = %.2f mV) ---\n', Q*1e3);
printf('Offset: %.2f LSB\n', offsetLsb);
printf('Ganho: %.2f LSB (%.2f%%)\n', ganhoLsb, ganhoPct);
printf('DNL: %.2f a %.2f LSB\n', min(dnl), max(dnl));
printf('INL: %.2f a %.2f LSB\n', min(inl), max(inl));
printf('Amostras por código (média): %.1f -> resolução da DNL ~ %.2f LSB\n', ...
       amostrasPorCodigo, 1/amostrasPorCodigo);
printf('Códigos ausentes: %d\n', numel(ausentes));
if ~isempty(ausentes)
  printf('  Códigos: %s\n', mat2str(ausentes));
end
printf('Quedas > %d LSB nas subidas: %d\n', TOL_QUEDA_LSB, quedas);
printf('Monotônico: %d (1 = sim)\n', monotonico);

printf('\n--- Características dinâmicas ---\n');
printf('Clock do ADC: %.0f kHz\n', F_ADC/1e3);
printf('Tempo de conversão (teórico): %.1f us\n', T_CONV*1e6);
printf('Taxa máxima teórica: %.0f SPS\n', FS_MAX);
printf('Taxa de amostragem nominal: %.0f SPS\n', FS);
printf('Taxa de amostragem medida: %.0f SPS\n', fsMedida);
printf('Faixa dinâmica (%d bits): %.1f dB\n', NBITS, 20*log10(NCODIGOS));
