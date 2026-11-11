"""
TEC513 - MI de Processamento de Sinais - UEFS
Problema 02 - Análise do conversor A/D do ATmega2560 (8 bits)

Versão em Python do script do Octave. Mesmas etapas e mesmas convenções:
  1) Recepção das amostras (serial ou arquivo .csv)
  2) Triangular ideal simulada, alinhada ao sinal capturado
  3) Comparação no tempo (estilo osciloscópio)
  4) Comparação na frequência (FFT)
  5) Erros estáticos: offset, ganho, DNL, INL, códigos ausentes, monotonicidade
  6) Resumo dos resultados e características dinâmicas

Convenções (datasheet ATmega2560, seção 26.6.4, p. 277-278):
  - Transição ideal k-1 -> k em (k - 0,5)·Q  (1ª transição em 0,5 LSB)
  - Offset: 1ª transição real x ideal, em LSB
  - Ganho: após remover o offset, última transição real x ideal
    (ideal = 1,5 LSB abaixo do fundo de escala), em LSB
  - INL e DNL pelo método dos extremos

Instalação:  pip install -r requirements.txt
Execução:    python analise_adc.py
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import sawtooth
from scipy.optimize import minimize


# ===================== PARÂMETROS (ajustar antes de rodar) ==================
# --- Origem dos dados ---
MODO = "arquivo"                       # "serial" = captura com o Arduino
                                       # "arquivo" = usa um .csv já salvo
PORTA = "/dev/ttyUSB0"                 # Linux  (Windows: "COM5")
ARQUIVO = "captura_1_2.5V_50Ohms.csv"  # onde a captura é salva (ou lida)

# --- Medidas feitas no laboratório ---
VREF = 5.0         # AVcc medido com multímetro (V)
VIN_MIN = 0.0      # mínimo da triangular medido no osciloscópio (V)
VIN_MAX = 5.0      # máximo da triangular medido no osciloscópio (V)
F_ENTRADA = 100    # frequência ajustada no gerador (Hz)

# --- Configuração do Arduino (deve bater com o .ino) ---
F_CPU = 16e6           # clock do ATmega2560 (Hz)
PRESCALER_TMR1 = 8     # prescaler do Timer1 (CS11)
OCR1A = 38             # valor de comparação do Timer1 (modo CTC)
PRESCALER_ADC = 16     # prescaler do ADC (ADPS2) -> clock do ADC = 1 MHz
CICLOS_CONV = 13.5     # ciclos por conversão em auto-trigger (Tabela 26-1)
NBITS = 8              # resolução usada (ADLAR = 1, leitura só de ADCH)
BAUD = 2_000_000       # igual ao Serial.begin() do .ino

# --- Ajustes da análise ---
AMOSTRAS_DESCARTE = 10  # 1ª conversão leva 25 ciclos; ADC ainda estabilizando
TOL_QUEDA_LSB = 1       # nas subidas, quedas de até 1 LSB = ruído
PASSO_FASE = 0.01       # resolução da busca inicial de fase (fração do período)

# --- Ajustes dos gráficos ---
JANELA_TEMPO_MS = 30    # trecho mostrado na visão geral (ms)
MEIA_JANELA = 8         # zoom: amostras antes/depois do ponto central
F_MAX_GRAFICO = 1000    # limite do eixo de frequência (Hz)
PISO_DB = -100          # piso dos gráficos em dB
N_HARMONICAS = 10       # harmônicas mostradas (1f a 10f)
N_BORDA = 8             # códigos mostrados no início/fim da curva de transferência
DPI = 150               # resolução das figuras salvas

# --- Grandezas derivadas (calculadas, não digitadas) ---
FS = F_CPU / (PRESCALER_TMR1 * (1 + OCR1A))   # taxa de amostragem (SPS)
F_ADC = F_CPU / PRESCALER_ADC                  # clock do ADC (Hz)
T_CONV = CICLOS_CONV / F_ADC                   # tempo de conversão (s)
FS_MAX = 1 / T_CONV                            # taxa máxima teórica (SPS)
NCODIGOS = 2 ** NBITS                          # 256 códigos (0 a 255)
COD_MAX = NCODIGOS - 1                         # código de fundo de escala
Q = VREF / NCODIGOS                            # degrau de quantização (V)

# Tensões ideais das transições k-1 -> k (k = 1..255), convenção do datasheet
V_TRANS_IDEAL = (np.arange(1, NCODIGOS) - 0.5) * Q


# ===================== FUNÇÕES AUXILIARES ===================================
def triangular(t, f, fase, vmin, vmax):
    """Triangular entre vmin e vmax. 'fase' em fração do período (0 a 1).

    sawtooth(x, 0.5) vai de -1 a +1: SOBE na 1ª metade do período e DESCE na
    2ª. Aqui ela é levada para 0..1 e depois para vmin..vmax.
    """
    forma = (sawtooth(2 * np.pi * (f * t + fase), 0.5) + 1) / 2
    return vmin + (vmax - vmin) * forma


def em_db(x, ref):
    """Magnitude em dB relativa a 'ref' (o eps evita log de zero)."""
    return 20 * np.log10(x / ref + np.finfo(float).eps)


def salvar(fig, nome):
    fig.tight_layout()
    fig.savefig(nome, dpi=DPI)


def receber_serial():
    """Abre a porta (o que reinicia o Arduino), espera a captura e lê o texto.

    O Arduino captura o buffer inteiro e só depois transmite, um número por
    linha. A leitura termina quando a porta fica 1 s sem dados (timeout).
    """
    import time
    import serial  # pyserial: só é necessário no modo "serial"

    TIMEOUT_S = 1.0       # sem dados por 1 s = fim da transmissão
    PAUSA_RESET_S = 2.0   # tempo para o Arduino reiniciar e capturar
    TAM_BLOCO = 4096      # bytes lidos por vez

    with serial.Serial(PORTA, BAUD, bytesize=8, parity="N",
                       stopbits=1, timeout=TIMEOUT_S) as s:
        s.reset_input_buffer()
        time.sleep(PAUSA_RESET_S)
        dados = bytearray()
        while bloco := s.read(TAM_BLOCO):
            dados += bloco

    # "123\r\n124\r\n..." -> [123, 124, ...]  (linhas vazias são ignoradas)
    codigos = np.array([int(x) for x in dados.decode().split() if x.isdigit()])
    print(f"Amostras recebidas: {codigos.size}")
    np.savetxt(ARQUIVO, codigos, fmt="%d")   # guarda a captura
    return codigos


# Faz o script trabalhar na própria pasta (onde ficam o .csv e as figuras)
PASTA = Path(__file__).resolve().parent
ARQUIVO = PASTA / ARQUIVO


# ===================== 1) RECEPÇÃO DOS DADOS (requisitos 3 e 4) =============
if MODO == "serial":
    codigos = receber_serial()
else:
    # ravel() deixa o vetor 1D, venha o .csv em linha ou em coluna
    codigos = np.loadtxt(ARQUIVO, delimiter=",", ndmin=1).ravel()

codigos = codigos.astype(int)[AMOSTRAS_DESCARTE:]
N = codigos.size
t = np.arange(N) / FS      # instante de cada amostra (s)
t_ms = t * 1e3             # o mesmo, em ms (para os gráficos)
v_real = codigos * Q       # tensão equivalente de cada código


# ===================== 2) TRIANGULAR IDEAL SIMULADA (requisito 5) ===========
# O gerador nunca está exatamente na frequência nominal e a captura começa
# numa fase qualquer. A frequência e a fase da triangular ideal são AJUSTADAS
# ao sinal real, minimizando o erro quadrático:
#   Passo 1: busca grosseira da fase, com a frequência nominal.
#   Passo 2: Nelder-Mead (o mesmo método do fminsearch) refina as duas.
def erro_quad(p):
    return np.sum((v_real - triangular(t, p[0], p[1], VIN_MIN, VIN_MAX)) ** 2)

lista_fase = np.arange(0, 1 + PASSO_FASE, PASSO_FASE)
fase0 = lista_fase[np.argmin([erro_quad([F_ENTRADA, f]) for f in lista_fase])]

res = minimize(erro_quad, [F_ENTRADA, fase0], method="Nelder-Mead",
               options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 2000})
f_est, fase = res.x

# Entrada analógica ideal e a saída de um ADC IDEAL para ela.
# round() coloca as transições em (k - 0,5)·Q, como no datasheet.
v_ideal = triangular(t, f_est, fase, VIN_MIN, VIN_MAX)
codigo_ideal = np.clip(np.round(v_ideal / Q), 0, COD_MAX).astype(int)

print(f"Frequência estimada: {f_est:.3f} Hz")
print(f"Fase: {fase % 1:.3f} do período ({(fase % 1) * 360:.1f} graus)")


# ===================== 3) COMPARAÇÃO NO TEMPO (requisito 4) =================
# 3a) Visão geral: sinal capturado (degraus) sobre a triangular ideal.
# 3b) Zoom no meio de uma subida: entrada contínua, amostras, retenção
#     (sample-and-hold) real e saída do ADC ideal.
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))

ax1.step(t_ms, v_real, "b", where="post", label="Real")
ax1.plot(t_ms, v_ideal, "r--", label="Ideal")
ax1.set(xlim=(0, JANELA_TEMPO_MS), xlabel="Tempo (ms)", ylabel="Tensão (V)",
        title="Sinal capturado x triangular ideal")
ax1.grid(True)
ax1.legend()

# Ponto central do zoom: 1ª amostra em que a saída ideal SOBE e cruza o meio
# da faixa, com espaço para MEIA_JANELA amostras de cada lado.
meio = round((codigos.min() + codigos.max()) / 2)
cruza = np.flatnonzero((codigo_ideal[:-1] < meio) & (codigo_ideal[1:] >= meio)) + 1
cruza = cruza[(cruza >= MEIA_JANELA) & (cruza < N - MEIA_JANELA)]
janela = cruza[0] + np.arange(-MEIA_JANELA, MEIA_JANELA + 1)

ax2.plot(t_ms[janela], v_ideal[janela], "k", label="Entrada")
ax2.plot(t_ms[janela], v_real[janela], "bo", label="Amostras")
ax2.step(t_ms[janela], v_real[janela], "r", where="post", label="Retenção (real)")
ax2.step(t_ms[janela], codigo_ideal[janela] * Q, "g--", where="post",
         label="Ideal quantizado")
ax2.set(ylim=((codigos[janela].min() - 2) * Q, (codigos[janela].max() + 2) * Q),
        xlabel="Tempo (ms)", ylabel="Tensão (V)",
        title="Zoom: entrada, amostras e retenção")
ax2.grid(True)
ax2.legend(loc="upper left")
salvar(fig, PASTA / "fig1_osciloscopio.png")


# ===================== 4) COMPARAÇÃO NA FREQUÊNCIA (requisito 4) ============
# A FFT supõe que o trecho se repete. Cortando em um número INTEIRO de
# períodos, cada harmônica cai exatamente num bin (sem vazamento espectral):
# com n_periodos no trecho, a harmônica h fica no bin h·n_periodos.
# O nível médio (DC) é removido para não dominar o gráfico.
amostras_por_periodo = FS / f_est
n_periodos = int(N // amostras_por_periodo)
M = round(n_periodos * amostras_por_periodo)

x_real = np.abs(np.fft.rfft(v_real[:M] - v_real[:M].mean()))
x_ideal = np.abs(np.fft.rfft(v_ideal[:M] - v_ideal[:M].mean()))
freq = np.fft.rfftfreq(M, 1 / FS)   # eixo de frequência (só positivas)
ref = x_ideal.max()                 # 0 dB = fundamental ideal

# 4a) Espectro completo até F_MAX_GRAFICO
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.vlines(freq, PISO_DB, em_db(x_real, ref), colors="b", label="Real")
ax.vlines(freq, PISO_DB, em_db(x_ideal, ref), colors="r", label="Ideal")
ax.set(xlim=(0, F_MAX_GRAFICO), ylim=(PISO_DB, 5), xlabel="Frequência (Hz)",
       ylabel="Magnitude (dB)",
       title=f"Espectro: real x ideal ({n_periodos} períodos inteiros)")
ax.grid(True)
ax.legend()
salvar(fig, PASTA / "fig2_fft.png")

# 4b) Só as harmônicas (a triangular ideal tem apenas as ímpares, caindo com
#     1/h²; as pares ideais ficam no piso)
h = np.arange(1, N_HARMONICAS + 1)
bins = h * n_periodos
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.vlines(h * f_est, PISO_DB, em_db(x_real[bins], ref), colors="b", lw=2)
ax.plot(h * f_est, em_db(x_real[bins], ref), "bo", label="Real")
ax.plot(h * f_est, em_db(x_ideal[bins], ref), "rx", ms=10, mew=2, label="Ideal")
ax.set_xticks(h * f_est, [f"{k}f" for k in h])
ax.set(xlim=(0, (N_HARMONICAS + 1) * f_est), ylim=(PISO_DB, 5),
       xlabel="Frequência (Hz)", ylabel="Magnitude (dB)",
       title="Harmônicas: real x ideal")
ax.grid(True)
ax.legend()
salvar(fig, PASTA / "fig3_harmonicas.png")


# ===================== 5) ERROS ESTÁTICOS (requisito 5) =====================
# MÉTODO DO HISTOGRAMA: com entrada triangular, a tensão passa o mesmo tempo
# em cada faixa. Um ADC ideal produz o mesmo número de ocorrências para todos
# os códigos; um código "mais largo" aparece mais.

# --- 5.1 Separação das rampas ---
# rampa = número de meio-períodos decorridos. Pares = subidas, ímpares =
# descidas. A primeira e a última rampa são incompletas e ficam de fora.
rampa = np.floor(2 * (f_est * t + fase)).astype(int)
completa = (rampa > rampa[0]) & (rampa < rampa[-1])

# --- 5.2 Histograma: histograma[k] = número de amostras do código k ---
histograma = np.bincount(codigos[completa], minlength=NCODIGOS)
amostras_por_codigo = histograma[1:COD_MAX].mean()   # média, sem as bordas

# --- 5.3 Tensões de transição reais ---
# A fração acumulada de amostras até o código k-1 é a fração da faixa
# VIN_MIN..VIN_MAX abaixo da transição k-1 -> k. Por isso offset e ganho
# dependem de VIN_MIN e VIN_MAX terem sido bem medidos.
v_trans = VIN_MIN + (VIN_MAX - VIN_MIN) * np.cumsum(histograma[:COD_MAX]) / histograma.sum()

# --- 5.4 Offset (LSB): 1ª transição real x ideal (0,5 LSB) ---
offset_v = v_trans[0] - V_TRANS_IDEAL[0]
offset_lsb = offset_v / Q

# --- 5.5 Ganho: última transição, já sem o offset, x ideal ---
ganho_lsb = ((v_trans[-1] - offset_v) - V_TRANS_IDEAL[-1]) / Q
q_real = (v_trans[-1] - v_trans[0]) / (COD_MAX - 1)   # degrau médio real
ganho_pct = (q_real / Q - 1) * 100

# --- 5.6 DNL (códigos 1..254): largura de cada degrau / degrau médio - 1 ---
# Usar o degrau médio real (e não Q) remove offset e ganho da DNL.
dnl = np.diff(v_trans) / q_real - 1

# --- 5.7 INL (transições 1..255): distância até a reta entre os extremos ---
reta = v_trans[0] + np.arange(COD_MAX) * q_real
inl = (v_trans - reta) / q_real

# --- 5.8 Códigos ausentes (1..254; 0 e 255 são bordas e não contam) ---
ausentes = np.flatnonzero(histograma[1:COD_MAX] == 0) + 1

# --- 5.9 Monotonicidade ---
# Nas subidas, o código não pode cair. Quedas de até TOL_QUEDA_LSB são ruído.
# Só compara amostras vizinhas da MESMA rampa.
# (Código ausente NÃO torna o ADC não monotônico: é outro erro.)
na_subida = (rampa % 2 == 0) & completa
cod_subida = codigos[na_subida]
mesma_rampa = np.diff(rampa[na_subida]) == 0
quedas = int(np.sum((np.diff(cod_subida) < -TOL_QUEDA_LSB) & mesma_rampa))
monotonico = quedas == 0

# --- Figura: histograma ---
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.bar(np.arange(NCODIGOS), histograma, width=1)
ax.set(xlim=(-1, NCODIGOS), xlabel="Código", ylabel="Ocorrências",
       title=f"Histograma de códigos ({ausentes.size} ausentes, "
             f"média {amostras_por_codigo:.1f} amostras/código)")
ax.grid(True)
salvar(fig, PASTA / "fig4_histograma.png")

# --- Figura: DNL e INL ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
ax1.vlines(np.arange(1, COD_MAX), 0, dnl)
ax1.set(ylabel="DNL (LSB)", title=f"DNL: {dnl.min():.2f} a {dnl.max():.2f} LSB")
ax1.grid(True)
ax2.plot(np.arange(1, NCODIGOS), inl, "r")
ax2.set(xlim=(0, COD_MAX), xlabel="Código", ylabel="INL (LSB)",
        title=f"INL: {inl.min():.2f} a {inl.max():.2f} LSB")
ax2.grid(True)
salvar(fig, PASTA / "fig5_dnl_inl.png")

# --- Figura: curva de transferência no início (offset) e no fim (ganho) ---
# Cada degrau começa na tensão de transição v_trans[k-1] e vai para o código k.
cod_ini = np.arange(0, N_BORDA + 1)
cod_fim = np.arange(COD_MAX - N_BORDA + 1, NCODIGOS)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
ax1.step(np.r_[0, v_trans[cod_ini[1:] - 1]], cod_ini, "b", where="post", label="Real")
ax1.step(np.r_[0, V_TRANS_IDEAL[cod_ini[1:] - 1]], cod_ini, "r--", where="post",
         label="Ideal")
ax1.set(xlabel="Entrada (V)", ylabel="Código", title=f"Offset: {offset_lsb:.2f} LSB")
ax1.grid(True)
ax1.legend(loc="upper left")

ax2.step(v_trans[cod_fim - 1], cod_fim, "b", where="post", label="Real")
ax2.step(V_TRANS_IDEAL[cod_fim - 1], cod_fim, "r--", where="post", label="Ideal")
ax2.set(xlabel="Entrada (V)", title=f"Ganho: {ganho_lsb:.2f} LSB ({ganho_pct:.2f}%)")
ax2.grid(True)
salvar(fig, PASTA / "fig6_offset_ganho.png")


# ===================== 6) RESUMO (requisitos 5 e 6) ========================
# A taxa "medida" compara a frequência nominal do gerador com a encontrada
# pelo ajuste: se o gerador está certo, a diferença vem da taxa real do
# Arduino. O tempo de conversão NÃO é medido aqui: isso é feito no Arduino
# (ex.: alternar um pino na ISR do ADC e medir no osciloscópio).
fs_medida = FS * F_ENTRADA / f_est

print(f"\n--- Erros estáticos (Q = {Q * 1e3:.2f} mV) ---")
print(f"Offset: {offset_lsb:.2f} LSB")
print(f"Ganho: {ganho_lsb:.2f} LSB ({ganho_pct:.2f}%)")
print(f"DNL: {dnl.min():.2f} a {dnl.max():.2f} LSB")
print(f"INL: {inl.min():.2f} a {inl.max():.2f} LSB")
print(f"Amostras por código (média): {amostras_por_codigo:.1f} "
      f"-> resolução da DNL ~ {1 / amostras_por_codigo:.2f} LSB")
print(f"Códigos ausentes: {ausentes.size}" + (f"  {ausentes.tolist()}" if ausentes.size else ""))
print(f"Quedas > {TOL_QUEDA_LSB} LSB nas subidas: {quedas}")
print(f"Monotônico: {'sim' if monotonico else 'não'}")

print("\n--- Características dinâmicas ---")
print(f"Clock do ADC: {F_ADC / 1e3:.0f} kHz")
print(f"Tempo de conversão (teórico): {T_CONV * 1e6:.1f} us")
print(f"Taxa máxima teórica: {FS_MAX:.0f} SPS")
print(f"Taxa de amostragem nominal: {FS:.0f} SPS")
print(f"Taxa de amostragem medida: {fs_medida:.0f} SPS")
print(f"Faixa dinâmica ({NBITS} bits): {20 * np.log10(NCODIGOS):.1f} dB")

plt.show()
