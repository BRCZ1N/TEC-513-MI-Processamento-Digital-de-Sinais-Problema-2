//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%TEC513-PDS-UEFS
//%Problema 02
//%Arquivo modificado: Prescaler em 16 + Saída física em 8 bits
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

const byte adcPin = 0;              // A0  <=  ATMega 2560 A0 pin for ADC0

const int MAX_RESULTS = 2560;       // size for store data to be sent

// CONFIGURADO PARA BYTE: Armazena estritamente 8 bits (0 a 255)
volatile byte results [MAX_RESULTS]; 
volatile int resultNumber; 

// ADC complete ISR
ISR (ADC_vect)
  {
  if (resultNumber >= MAX_RESULTS)
    ADCSRA = 0;  // turn off ADC
  else
    // CONFIGURADO PARA 8 BITS: Captura diretamente o registrador HIGH
    results [resultNumber++] = ADCH; 
  }  // end of ADC_vect
  
EMPTY_INTERRUPT (TIMER1_COMPB_vect);

// ADC configure initialization
void ADC_init() {
  ADCSRA =  bit (ADEN) | bit (ADIE) | bit (ADIF);   // turn ADC on, want interrupt on completion
  
  // CONFIGURADO: Prescaler de 16 (Clock do ADC = 16MHz/16 = 1MHz)
  ADCSRA |= bit (ADPS2);  
  
  // CONFIGURADO: Vref em 5V (REFS0) + Justificação à Esquerda (ADLAR) para isolar os 8 bits no ADCH
  ADMUX = bit (REFS0) | bit (ADLAR) | (adcPin & 7); 
  
  ADCSRB = bit (ADTS0) | bit (ADTS2);  // Timer/Counter1 Compare Match B
  ADCSRA |= bit (ADATE);   // turn on automatic triggering
} 

// Timer sampling configure
void timer(){
  TCCR1A = 0;
  TCCR1B = 0;
  TCNT1 = 0;
  
  TCCR1B = bit (CS11) | bit (WGM12);  // CTC, prescaler of 8
  TIMSK1 = bit (OCIE1B);  
  
  // Frequência de amostragem controlada pelo Timer = 2 kHz
  OCR1A = 999;    
  OCR1B = 999;   
}

// Setup and run
void setup ()
  {
  Serial.begin (115200);
  Serial.println ();
  
  timer();
  ADC_init();
  
  // wait for buffer to fill
  while (resultNumber < MAX_RESULTS)
    { }
    
  // Envia os dados convertidos de 8 bits pelo canal Serial
  for (int i = 0; i < MAX_RESULTS; i++)
    Serial.println (results [i]);

  }  // end of setup

void loop () { }
