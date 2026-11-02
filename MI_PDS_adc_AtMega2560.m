%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TEC513-MI de PDS-UEFS
%Problema 02
%Arquivo para teste na recepção de dados pela USB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%% LIMPA E FECHA TUDO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;							%limpa a tela do octave
clear all;				%limpa todas as variáveis do octave
close all;				%fecha todas as janelas

%%%%%%%%%%%%%%%%%%% CHAMADA DAS BIBLIOTECAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pkg load signal					%biblioteca para processamento de sinais
pkg load instrument-control		%biblioteca para comunicação serial

%%%%%%%%%%%%%%%%%%% ALOCAÇÃO DE VARIÁVEIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAX_RESULTS=2560;
  fs=50000;    % Aqui ajusta a frequência de amostragem usada no processo na ADC
amostras = 30000;  % Quantidade a amostras que irá usar para visualizar na tela
raw = [];				 % variavel para armazenar os dados cru recebido pela USB (raw)

%%%%%%%%%%%%%%%%%%% ABERTURA DA PORTA SERIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s1 = serial("COM3"); 	    		%Abre a porta serial que esta no microcontrolador
set(s1,'baudrate',115200);		%velocidade de transmissão 2Mbps
set(s1,'bytesize',8); 			  %8 bits de dados
set(s1,'parity','n'); 			  %sem paridade ('y' 'n')
set(s1,'stopbits',1); 			  %1 bit de parada (1 ou 2)
set(s1,'timeout',200);			  %tempo ocioso sem conecção 20.0 segundos
pause(1);				          		%espera 1 segundo antes de ler dado

%%%%%%%%%%%%%%%%%%% LEITURA DA MENSAGEM INICIAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;							%primeiro índice de leitura
while(1)						%espera para ler a mensagem inicial
	t(i) = fread(s1,1);			%le as amostras de uma em uma
	if (t(i)==10)				%se for lido um enter (10 em ASC2)
		break;					%sai do loop
	endif
	i = i+1;					%incrementa o índice de leitura
end
c = char(t); 					%transformando caracteres recebidos em string
printf('recebido: %s',c);		%imprime na tela do octave o que foi recebido

%%%%%%%%%%%%%%%%%%% CAPTURA DAS AMOSTRAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);						%cria uma figura
tic;							    %captura do tempo inicial
	data = fread(s1,amostras,'uint8');	%captura das amostras
  raw = cat(2,raw,data);			  %armazena o dado bruto (raw = sem processamento)
  x=char(raw);                  %converte em carateres os dados recebidos
  d=str2num(x);                 %string de carateres em números
  x=double(d);                  %os números inteiros em double se necessário paras a próximas etapas
  time=[(0:1:length(d)-1)/fs];  %tamanhos do dominio normalizado

  subplot(3,1,1);           %plotando as figuras
  plot(time,d*5/1023)   %plota as amostras interpoladas
  xlabel('t(s)');
  title('Sinal gerado x(t)');
  subplot(3,1,2);
  stem(d);		          %plota a janela de amostras
  xlabel('n');
  title('x[n]');
  subplot(3,1,3);
  stairs(d);	        	%plota a janela de amostras reguradas
  xlabel('n');
  title('x[n] segurado');
	hold on;							%mantem as amostras anteriores
toc;							  %captura do tempo final

%%%%%%%%%%%%%%%%%%% FECHA A PORTA DE COMUNICAÇÃO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fclose(s1);
%%%%%%%%%%%%%%%%%%% PÓS-PROCESSAMENTO E ANÁLISE DE ERROS %%%%%%
disp('Iniciando análise de erros estáticos...');

% 1. Parâmetros do ADC
n_bits = 8;                     % Resolução do ATmega2560
niveis = 2^n_bits;               % 1024 níveis (0 a 1023)
Vref = 5.0;                      % Tensão de referência em Volts
Q = Vref / niveis;               % Passo ideal ou LSB (Largura do degrau)

% Garantir que 'd' é um vetor coluna e remover possíveis NaNs da conversão de string
d = d(~isnan(d));

% 2. Método do Histograma para calcular DNL e INL
% Conta a ocorrência de cada código (de 0 a 1023) nas amostras recebidas
[ocorrencias, codigos_adc] = hist(d, 0:(niveis-1));

% Para não distorcer a média por saturação, removemos os códigos extremos (0 e 1023)
ocorrencias_validas = ocorrencias(2:end-1);
media_ocorrencias = mean(ocorrencias_validas);

% 3. Cálculo da DNL (Não Linearidade Diferencial) em LSB
% DNL = (Largura Real / Largura Ideal) - 1
DNL = (ocorrencias / media_ocorrencias) - 1;
DNL(1) = 0;   % Zera os extremos que sofrem "clipping" da onda
DNL(end) = 0;

% 4. Cálculo da INL (Não Linearidade Integral) em LSB
% A INL é a soma cumulativa dos erros da DNL
INL = cumsum(DNL);

% 5. Verificação da Monotonicidade
% Um conversor é monotônico se nenhum código digital for menor que o anterior
% para uma tensão maior (Na DNL, isso significa DNL > -1 LSB sempre).
if all(DNL > -1)
    disp('-> Monotonicidade: OK (O conversor é Monotônico)');
else
    disp('-> Monotonicidade: FALHA (Foi detectado DNL <= -1 LSB)');
end

% 6. Cálculo do Erro de Offset e Ganho
% Vamos reconstruir a tensão real do conversor usando o desvio da INL
tensao_ideal = (0:niveis-1)' .* Q;
tensao_real = tensao_ideal + (INL' .* Q);

% Usamos regressão linear (polyfit de grau 1: y = mx + b)
% para achar a reta real do conversor
coeficientes = polyfit(0:niveis-1, tensao_real', 1);
ganho_real = coeficientes(1);   % Inclinação real
offset_real = coeficientes(2);  % Interceptação com o eixo Y

% O erro é a diferença entre a reta real e a reta ideal teórica
erro_offset = offset_real - 0.0; % Idealmente cruzaria o zero
erro_ganho = ganho_real - Q;     % Idealmente o ganho seria 1 Q por código

printf('-> Erro de Offset: %f V\n', erro_offset);
printf('-> Erro de Ganho: %f V por transição\n', erro_ganho);
printf('-> DNL Máxima: %f LSB\n', max(abs(DNL)));
printf('-> INL Máxima: %f LSB\n', max(abs(INL)));


figure(2);

% Subplot 1: DNL
subplot(2,1,1);
plot(codigos_adc, DNL, 'b', 'LineWidth', 1);
title('Não Linearidade Diferencial (DNL)');
xlabel('Código Digital de Saída');
ylabel('DNL (LSB)');
xlim([0 1023]);
grid on;

% Subplot 2: INL
subplot(2,1,2);
plot(codigos_adc, INL, 'r', 'LineWidth', 1.5);
title('Não Linearidade Integral (INL)');
xlabel('Código Digital de Saída');
ylabel('INL (LSB)');
xlim([0 1023]);
grid on;

% Figura 3: Comparação da Curva de Transferência (Real vs Ideal)
figure(3);
plot(codigos_adc, tensao_ideal, 'k--', 'LineWidth', 2); hold on;
plot(codigos_adc, tensao_real, 'r-', 'LineWidth', 1);
title('Curva de Transferência de Tensão: Simulação Ideal vs Real');
xlabel('Código Digital');
ylabel('Tensão Equivalente (V)');
legend('Reta Ideal (Simulação)', 'Curva Real (Medida)');
grid on;
