% Demonstração simples de aliasing no ESPECTRO de uma onda triangular
%
% A triangular só tem harmônicas ÍMPARES: 100, 300, 500, 700, 900 Hz...
% Se fs/2 (Nyquist) fica abaixo de alguma delas, essa harmônica "dobra"
% e aparece numa frequência errada.
clc; clear; close all;
pkg load signal                          % função sawtooth

f = 100;                                 % frequência da triangular (Hz)

% sawtooth(x, 0.5) gera uma triangular de -1 a 1
tri = @(t) sawtooth(2*pi*f*t, 0.5);

% Caso 1: SEM aliasing (Nyquist = 5000 Hz, bem acima das harmônicas visíveis)
fs1 = 10000;
% Caso 2: COM aliasing (Nyquist = 375 Hz: 500, 700, 900 Hz... dobram)
fs2 = 750;

fsList  = [fs1 fs2];
titulos = {'Sem aliasing', 'Com aliasing'};

figure;
for i = 1:2
  fs = fsList(i);
  n  = 0:1/fs:1-1/fs;                    % 1 s de sinal -> resolução de 1 Hz
  x  = tri(n);

  % --- Tempo: amostras (primeiros 30 ms) ---
  subplot(2,2,2*i-1);
  mostra = n <= 0.03;
  stem(n(mostra)*1e3, x(mostra), 'filled', 'markersize', 3);
  grid on;
  ylim([-1.2 1.2]);
  xlabel('Tempo (ms)'); ylabel('Amplitude');
  title(sprintf('%s: fs = %d Hz', titulos{i}, fs));

  % --- Espectro ---
  X    = abs(fft(x - mean(x)));          % espectro (sem o nível DC)
  X    = 20*log10(X/max(X) + eps);       % em dB (0 dB = fundamental)
  freq = (0:numel(x)-1)*fs/numel(x);     % eixo de frequência (Hz)

  subplot(2,2,2*i);
  stem(freq, X, 'filled', 'markersize', 3, 'basevalue', -60);
  grid on;
  ylim([-60 5]);
  xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)');
  if i == 1
    xlim([0 1000]);
    title('Só harmônicas ímpares (100, 300, 500...)');
  else
    xlim([0 fs/2]);
    title('Picos falsos em 50, 150, 200, 250, 350 Hz');
  end
end

% Onde cada harmônica aparece depois de amostrada a fs2
printf('Harmônica (Hz)  ->  aparece em (Hz) com fs = %d Hz\n', fs2);
for k = 1:2:13
  fk = k*f;
  printf('%8d        ->  %6.0f\n', fk, abs(fk - fs2*round(fk/fs2)));
end
