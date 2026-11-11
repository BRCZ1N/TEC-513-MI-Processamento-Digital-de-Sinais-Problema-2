"""
TEC513 - MI de Processamento de Sinais - UEFS
Problema 02 - Análise do conversor A/D do ATmega2560 (8 bits)

Etapas:
  1) Leitura das amostras (serial ou .csv)
  2) Ajuste da triangular ideal ao sinal capturado
  3) Comparação no tempo (estilo osciloscópio)
  4) Comparação na frequência (FFT)
  5) Erros estáticos: offset, ganho, DNL, INL, códigos ausentes, monotonicidade

Convenção do datasheet (seção 26.6.4): a transição ideal do código k-1 para k
fica em (k - 0,5)·Q, ou seja, a primeira transição ideal está em 0,5 LSB.

Instalação:  pip install -r requirements.txt
Execução:    python analise_adc.py
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import sawtooth


# =============================================================================
# PARÂMETROS
# =============================================================================
MODO    = "arquivo"                     # "serial" ou "arquivo"
PORTA   = "/dev/ttyUSB0"                # Windows: "COM5"
ARQUIVO = "captura_1_2.5V_50Ohms.csv"
BAUD    = 2_000_000                     # igual ao Serial.begin() do .ino

# Medidas do laboratório
VREF      = 5.0     # AVcc (V)
VIN_MIN   = 0.0     # mínimo da triangular no osciloscópio (V)
VIN_MAX   = 5.0     # máximo da triangular no osciloscópio (V)
F_ENTRADA = 100     # frequência do gerador (Hz) - só referência; a real é estimada

# Configuração do Arduino (Timer1 em CTC define a taxa de amostragem)
F_CPU, PRESCALER_TMR1, OCR1A = 16e6, 8, 38
FS = F_CPU / (PRESCALER_TMR1 * (1 + OCR1A))   # ≈ 51 282 amostras/s

# Quantização
NBITS    = 8
NCODIGOS = 2 ** NBITS          # 256 códigos
COD_MAX  = NCODIGOS - 1        # 255
Q        = VREF / NCODIGOS     # degrau (V)

# Análise
AMOSTRAS_DESCARTE = 10   # início da captura (1ª conversão é mais longa)
TOL_QUEDA_LSB     = 1    # queda tolerada nas subidas (ruído)

# Gráficos
JANELA_MS     = 30       # trecho da visão geral (ms)
MEIA_JANELA   = 8        # amostras de cada lado no zoom
F_MAX_FFT     = 1000     # eixo do espectro (Hz)
PISO_DB       = -100
N_BORDA       = 8        # degraus mostrados no início/fim da curva
PASTA = Path(__file__).resolve().parent


# =============================================================================
# FUNÇÕES AUXILIARES
# =============================================================================
def triangular(t, f, fase):
    """Triangular de VIN_MIN a VIN_MAX; sobe na 1ª metade do período."""
    forma = (sawtooth(2 * np.pi * (f * t + fase), 0.5) + 1) / 2   # 0..1
    return VIN_MIN + (VIN_MAX - VIN_MIN) * forma


def espectro(v):
    """Módulo da FFT (só frequências positivas), sem o nível DC."""
    return np.abs(np.fft.rfft(v - v.mean()))


def em_db(x, ref):
    return 20 * np.log10(x / ref + 1e-12)


def salvar(fig, nome):
    fig.tight_layout()
    fig.savefig(PASTA / nome, dpi=150)


# =============================================================================
# 1) LEITURA DOS DADOS
# =============================================================================
def ler_serial():
    """Abre a porta (reinicia o Arduino), espera a captura e lê até 1 s de silêncio."""
    import time
    import serial

    with serial.Serial(PORTA, BAUD, timeout=1) as s:
        s.reset_input_buffer()
        time.sleep(2)                       # Arduino reinicia e captura
        dados = bytearray()
        while bloco := s.read(4096):
            dados += bloco

    codigos = np.array([int(x) for x in dados.decode().split() if x.isdigit()])
    np.savetxt(PASTA / ARQUIVO, codigos, fmt="%d")
    return codigos


def ler_dados():
    if MODO == "serial":
        codigos = ler_serial()
    else:
        codigos = np.loadtxt(PASTA / ARQUIVO, delimiter=",", ndmin=1)
    return codigos.ravel().astype(int)[AMOSTRAS_DESCARTE:]


# =============================================================================
# 2) AJUSTE DA TRIANGULAR IDEAL
# =============================================================================
def ajustar_triangular(t, v_real):
    """Estima frequência e fase da triangular sem otimização numérica.

    A fase do sinal cresce em linha reta com o tempo: fase(t) = f·t + fase0.
    Basta, então, descobrir a fase de cada amostra e ajustar uma reta.

    Passo 1 (estimativa grosseira): instantes em que o sinal cruza o nível
    médio SUBINDO. Isso acontece sempre em 1/4 do período, então o
    espaçamento entre os cruzamentos dá o período e o primeiro dá a fase.

    Passo 2 (refinamento): numa rampa, a tensão é proporcional à posição no
    período. Cada amostra é convertida em fase (subida: 0 -> 0,5; descida:
    0,5 -> 1) e uma regressão linear sobre todas as amostras dá f e fase0.
    Os trechos perto dos picos (10% de cada lado) ficam de fora.
    """
    u = (v_real - VIN_MIN) / (VIN_MAX - VIN_MIN)       # tensão normalizada 0..1

    # Passo 1: cruzamentos de u = 0,5 subindo, com interpolação linear
    i = np.flatnonzero((u[:-1] < 0.5) & (u[1:] >= 0.5))
    t_cruz = t[i] + (0.5 - u[i]) / (u[i + 1] - u[i]) / FS
    periodo, t0 = np.polyfit(np.arange(t_cruz.size), t_cruz, 1)
    f, fase = 1 / periodo, 0.25 - t0 / periodo

    # Passo 2: fase de cada amostra e regressão linear fase x tempo
    rampa = np.floor(2 * (f * t + fase))               # par = subida
    sobe = rampa % 2 == 0
    fase_amostra = np.where(sobe, rampa / 2 + u / 2, (rampa + 1) / 2 - u / 2)
    usar = (u > 0.1) & (u < 0.9)
    f, fase = np.polyfit(t[usar], fase_amostra[usar], 1)
    return f, fase % 1


# =============================================================================
# 3) COMPARAÇÃO NO TEMPO
# =============================================================================
def grafico_tempo(t_ms, v_real, v_ideal, cod_ideal):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7))

    # Visão geral
    ax1.step(t_ms, v_real, "b", where="post", label="Real")
    ax1.plot(t_ms, v_ideal, "r--", label="Ideal")
    ax1.set(xlim=(0, JANELA_MS), xlabel="Tempo (ms)", ylabel="Tensão (V)",
            title="Sinal capturado x triangular ideal")

    # Zoom: primeira subida que cruza o meio da faixa
    meio = NCODIGOS // 2
    cruza = np.flatnonzero((cod_ideal[:-1] < meio) & (cod_ideal[1:] >= meio)) + 1
    centro = cruza[cruza >= MEIA_JANELA][0]
    j = slice(centro - MEIA_JANELA, centro + MEIA_JANELA + 1)

    ax2.plot(t_ms[j], v_ideal[j], "k", label="Entrada")
    ax2.plot(t_ms[j], v_real[j], "bo", label="Amostras")
    ax2.step(t_ms[j], v_real[j], "r", where="post", label="Retenção (real)")
    ax2.step(t_ms[j], cod_ideal[j] * Q, "g--", where="post", label="Ideal quantizado")
    ax2.set(xlabel="Tempo (ms)", ylabel="Tensão (V)",
            title="Zoom: entrada, amostras e retenção")

    ax1.legend(loc="upper right")
    ax2.legend(loc="upper left")
    for ax in (ax1, ax2):
        ax.grid(True)
    salvar(fig, "fig1_osciloscopio.png")


# =============================================================================
# 4) COMPARAÇÃO NA FREQUÊNCIA
# =============================================================================
def grafico_frequencia(v_real, v_ideal, f_est):
    """FFT em número inteiro de períodos (evita vazamento espectral)."""
    n_periodos = int(len(v_real) // (FS / f_est))
    M = round(n_periodos * FS / f_est)

    x_real = espectro(v_real[:M])
    x_ideal = espectro(v_ideal[:M])
    ref = x_ideal.max()                     # 0 dB = fundamental ideal
    freq = np.fft.rfftfreq(M, 1 / FS)

    # Espectro completo
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.vlines(freq, PISO_DB, em_db(x_real, ref), colors="b", label="Real")
    ax.vlines(freq, PISO_DB, em_db(x_ideal, ref), colors="r", label="Ideal")
    ax.set(xlim=(0, F_MAX_FFT), ylim=(PISO_DB, 5), xlabel="Frequência (Hz)",
           ylabel="Magnitude (dB)", title=f"Espectro ({n_periodos} períodos)")
    ax.grid(True)
    ax.legend()
    salvar(fig, "fig2_fft.png")


# =============================================================================
# 5) ERROS ESTÁTICOS (método do histograma)
# =============================================================================
def erros_estaticos(codigos, t, f_est, fase):
    """Com entrada triangular, todos os códigos deveriam aparecer igualmente.

    Código mais largo = mais ocorrências. A fração acumulada de ocorrências
    dá a tensão de cada transição dentro da faixa VIN_MIN..VIN_MAX.
    """
    # Rampas: cada meio-período recebe um número (par = subida, ímpar = descida).
    # A primeira e a última são incompletas e ficam de fora.
    rampa = np.floor(2 * (f_est * t + fase)).astype(int)
    completa = (rampa > rampa[0]) & (rampa < rampa[-1])

    hist = np.bincount(codigos[completa], minlength=NCODIGOS)

    # Tensões das transições k-1 -> k (k = 1..255): real e ideal
    v_trans = VIN_MIN + (VIN_MAX - VIN_MIN) * np.cumsum(hist[:COD_MAX]) / hist.sum()
    v_ideal = (np.arange(1, NCODIGOS) - 0.5) * Q

    # Offset e ganho (em LSB, como no datasheet)
    offset_v = v_trans[0] - v_ideal[0]
    ganho_v  = (v_trans[-1] - offset_v) - v_ideal[-1]

    # DNL e INL pela reta entre a primeira e a última transição
    q_real = (v_trans[-1] - v_trans[0]) / (COD_MAX - 1)
    dnl = np.diff(v_trans) / q_real - 1                              # códigos 1..254
    inl = (v_trans - (v_trans[0] + np.arange(COD_MAX) * q_real)) / q_real   # transições 1..255

    # Códigos ausentes (sem contar as bordas 0 e 255)
    ausentes = np.flatnonzero(hist[1:COD_MAX] == 0) + 1

    # Monotonicidade: nas subidas, o código não pode cair mais que a tolerância
    sobe = (rampa % 2 == 0) & completa
    mesma_rampa = np.diff(rampa[sobe]) == 0
    quedas = np.sum((np.diff(codigos[sobe]) < -TOL_QUEDA_LSB) & mesma_rampa)

    return dict(hist=hist, v_trans=v_trans, v_ideal=v_ideal, sobe=sobe,
                offset=offset_v / Q, ganho=ganho_v / Q,
                ganho_pct=(q_real / Q - 1) * 100,
                dnl=dnl, inl=inl, ausentes=ausentes, quedas=int(quedas))


def graficos_erros(e):
    # Histograma
    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.bar(np.arange(NCODIGOS), e["hist"], width=1)
    ax.set(xlim=(-1, NCODIGOS), xlabel="Código", ylabel="Ocorrências",
           title=f"Histograma de códigos ({e['ausentes'].size} ausentes)")
    ax.grid(True)
    salvar(fig, "fig3_histograma.png")

    # DNL e INL
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax1.vlines(np.arange(1, COD_MAX), 0, e["dnl"])
    ax1.set(ylabel="DNL (LSB)",
            title=f"DNL: {e['dnl'].min():.2f} a {e['dnl'].max():.2f} LSB")
    ax2.plot(np.arange(1, NCODIGOS), e["inl"], "r")
    ax2.set(xlim=(0, COD_MAX), xlabel="Código", ylabel="INL (LSB)",
            title=f"INL: {e['inl'].min():.2f} a {e['inl'].max():.2f} LSB")
    for ax in (ax1, ax2):
        ax.grid(True)
    salvar(fig, "fig4_dnl_inl.png")

    # Curva de transferência: início (offset) e fim (ganho)
    ini = np.arange(1, N_BORDA + 1)              # códigos 1..8
    fim = np.arange(NCODIGOS - N_BORDA, NCODIGOS)  # códigos 248..255
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    for ax, cods, titulo in ((ax1, ini, f"Offset: {e['offset']:.2f} LSB"),
                             (ax2, fim, f"Ganho: {e['ganho']:.2f} LSB")):
        ax.step(e["v_trans"][cods - 1], cods, "b", where="post", label="Real")
        ax.step(e["v_ideal"][cods - 1], cods, "r--", where="post", label="Ideal")
        ax.set(xlabel="Entrada (V)", ylabel="Código", title=titulo)
        ax.grid(True)
        ax.legend(loc="upper left")
    salvar(fig, "fig5_offset_ganho.png")


def grafico_monotonicidade(codigos, v_entrada, sobe):
    """Monotonicidade nas subidas.

    Esquerda: código x tensão de entrada. Se for monotônico, a escada só sobe.
    Direita: variação do código entre amostras vizinhas (Δcódigo). Valores
    abaixo de -TOL_QUEDA_LSB seriam quedas reais (não monotonicidade).
    """
    delta = np.diff(codigos[sobe])
    delta = delta[np.abs(delta) < NCODIGOS // 2]   # ignora o salto entre rampas

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    ax1.plot(v_entrada[sobe], codigos[sobe], "b.", ms=2)
    ax1.set(xlabel="Entrada (V)", ylabel="Código",
            title="Código x entrada (subidas)")

    valores, contagem = np.unique(delta, return_counts=True)
    ax2.bar(valores, contagem, color=np.where(valores < -TOL_QUEDA_LSB, "r", "b"))
    ax2.axvline(-TOL_QUEDA_LSB - 0.5, color="r", ls="--", label="limite de queda")
    ax2.set(xlabel="Δcódigo entre amostras (LSB)", ylabel="Ocorrências",
            title="Variação do código nas subidas")
    ax2.legend()
    for ax in (ax1, ax2):
        ax.grid(True)
    salvar(fig, "fig6_monotonicidade.png")


# =============================================================================
# PROGRAMA PRINCIPAL
# =============================================================================
codigos = ler_dados()
t = np.arange(codigos.size) / FS
v_real = codigos * Q

f_est, fase = ajustar_triangular(t, v_real)
v_ideal = triangular(t, f_est, fase)
cod_ideal = np.clip(np.round(v_ideal / Q), 0, COD_MAX)   # ADC ideal

grafico_tempo(t * 1e3, v_real, v_ideal, cod_ideal)
grafico_frequencia(v_real, v_ideal, f_est)

e = erros_estaticos(codigos, t, f_est, fase)
graficos_erros(e)
grafico_monotonicidade(codigos, v_ideal, e["sobe"])

print(f"Amostras: {codigos.size}  |  Frequência ajustada: {f_est:.3f} Hz")
print(f"Offset:   {e['offset']:.2f} LSB")
print(f"Ganho:    {e['ganho']:.2f} LSB ({e['ganho_pct']:.2f}%)")
print(f"DNL:      {e['dnl'].min():.2f} a {e['dnl'].max():.2f} LSB")
print(f"INL:      {e['inl'].min():.2f} a {e['inl'].max():.2f} LSB")
print(f"Códigos ausentes: {e['ausentes'].tolist() or 'nenhum'}")
print(f"Monotônico: {'sim' if e['quedas'] == 0 else 'não'} "
      f"({e['quedas']} quedas > {TOL_QUEDA_LSB} LSB nas subidas)")

plt.show()