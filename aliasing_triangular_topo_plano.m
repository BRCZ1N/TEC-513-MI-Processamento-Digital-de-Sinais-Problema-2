% Demonstração de amostragem de TOPO PLANO (PAM) com uma onda triangular
%
% Cada amostra é mantida constante por um tempo tau (como o sample-and-hold
% do ADC). No espectro isso aparece como:
%   - réplicas do espectro em torno de fs, 2fs, 3fs...
%   - todas multiplicadas por uma envoltória sinc(f*tau) (efeito de abertura),
%     que zera em f = 1/tau
% O aliasing continua dependendo só de fs: com fs baixa, as réplicas invadem
% a faixa 0..fs/2 e criam picos em frequências que a triangular não tem.
clc; clear; close all;
pkg load signal                          % função sawtooth

f     = 100;                             % frequência da triangular (Hz)
FSIM  = 300000;                          % taxa da "simulação analógica" (Hz)
CICLO = 0.5;                             % tau = CICLO*Ts
PISO  = -60;                             % piso do espectro (dB)

% Casos: fs (Hz), janela de tempo mostrada (ms), xlim do espectro (Hz), nome
casos = { 10000, 2.5, 1000, 'Sem aliasing'; ...
            750, 30,  1500, 'Com aliasing'};

t = (0:FSIM-1)/FSIM;                     % 1 s de "sinal analógico"
x = sawtooth(2*pi*f*t, 0.5);             % triangular de -1 a 1

figure('position', [80 80 1300 750]);
for i = 1:2
  [fs, janela, fMax, nome] = casos{i,:};
  L   = FSIM/fs;                         % pontos da simulação por período Ts
  M   = round(CICLO*L);                  % pontos em que o pulso fica "ligado"
  tau = M/FSIM;                          % largura do pulso (s)

  % --- Topo plano: repete o valor da amostra por tau, depois zero ---
  pos  = mod(0:FSIM-1, L);               % posição dentro de cada Ts
  base = (0:FSIM-1) - pos;               % índice do instante de amostragem
  xp   = x(base + 1) .* (pos < M);       % sinal PAM de topo plano

  %% ---------- Tempo ----------
  subplot(2,2,2*i-1);
  mostra = t <= janela/1e3;
  amostr = mostra & pos == 0;
  plot(t(mostra)*1e3, x(mostra), 'k', 'linewidth', 1.2); hold on;
  plot(t(mostra)*1e3, xp(mostra), 'b', 'linewidth', 1.2);
  stem(t(amostr)*1e3, x(amostr), 'r', 'filled', 'markersize', 4, 'linestyle', 'none');
  grid on; box on;
  xlim([0 janela]);
  if i == 1
    ylim([-1.1 -0.4]);                   % zoom no começo da subida
  else
    ylim([-1.2 1.2]);
  end
  xlabel('Tempo (ms)'); ylabel('Amplitude');
  title(sprintf('%s: fs = %d Hz, \\tau = Ts/2 = %.2f ms', nome, fs, tau*1e3));
  if i == 1
    legend('Analógico', 'Topo plano', 'Amostras', 'location', 'northwest');
    text(0.05, -0.43, sprintf('(zoom: %.1f ms, %d amostras)', janela, sum(amostr)), ...
         'verticalalignment', 'top');
  end

  %% ---------- Espectro ----------
  X  = abs(fft(xp - mean(xp)));
  X  = 20*log10(X/max(X) + eps);         % 0 dB = fundamental
  fr = 0:FSIM-1;                         % Hz (1 s de sinal -> resolução de 1 Hz)

  % Separa os picos em três grupos
  pico    = fr <= fMax & X > PISO;
  impar   = mod(fr, 2*f) == f;           % 100, 300, 500... (o que a triangular tem)
  abaixo  = fr <= fs/2;
  certo   = pico & abaixo & impar;       % harmônica no lugar certo
  falso   = pico & abaixo & ~impar;      % alias dentro de 0..fs/2
  replica = pico & ~abaixo;              % acima de Nyquist (réplicas)

  subplot(2,2,2*i);
  h = []; rot = {};
  if any(replica)
    h(end+1) = stem(fr(replica), X(replica), 'color', [0.6 0.6 0.6], ...
                    'filled', 'markersize', 3, 'basevalue', PISO); hold on;
    rot{end+1} = 'Réplicas (acima de fs/2)';
  end
  h(end+1) = stem(fr(certo), X(certo), 'b', 'filled', 'markersize', 4, 'basevalue', PISO); hold on;
  rot{end+1} = 'Harmônicas da triangular';
  if any(falso)
    h(end+1) = stem(fr(falso), X(falso), 'r', 'filled', 'markersize', 4, 'basevalue', PISO);
    rot{end+1} = 'Picos falsos (aliasing)';
  end

  if fs/2 < fMax                         % Nyquist e envoltória só no caso 2
    fEnv = linspace(1, fMax, 1000);
    env  = 20*log10(abs(sinc(fEnv*tau))/abs(sinc(f*tau)) + eps);
    h(end+1) = plot(fEnv, env, 'm--', 'linewidth', 1.2);
    rot{end+1} = 'Envoltória sinc(f\tau)';
    h(end+1) = plot([fs/2 fs/2], [PISO 5], 'k:', 'linewidth', 1.5);
    rot{end+1} = sprintf('Nyquist (%d Hz)', fs/2);
  end

  grid on; box on;
  xlim([0 fMax]); ylim([PISO 5]);
  xlabel('Frequência (Hz)'); ylabel('Magnitude (dB)');
  if i == 1
    title('Espectro: só as harmônicas ímpares (réplicas em torno de 10 kHz)');
  else
    title('Espectro: réplicas em fs \pm 100 Hz e picos falsos abaixo de fs/2');
  end
  legend(h, rot, 'location', 'northeast');
end
