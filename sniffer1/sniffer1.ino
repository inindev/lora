
// Copyright (C) 2024, John Clark <inindev@gmail.com>

#include "Arduino.h"
#include "Wire.h"
#include "LoRaWan_APP.h"
#include "payload_parse.h"


#define RF_FREQUENCY1                               910300000 // Hz
#define RF_FREQUENCY2                               923300000 // Hz

#define USE_COLOR_LOGGING                           true

#define TX_OUTPUT_POWER                             14        // dBm
#define LORA_BANDWIDTH                              0         // [0: 125 kHz,
                                                              //  1: 250 kHz,
                                                              //  2: 500 kHz,
                                                              //  3: reserved]
#define LORA_SPREADING_FACTOR                       7         // [SF7..SF12]
#define LORA_CODINGRATE                             1         // [1: 4/5,
                                                              //  2: 4/6,
                                                              //  3: 4/7,
                                                              //  4: 4/8]
#define LORA_PREAMBLE_LENGTH                        8         // same for Tx and Rx
#define LORA_SYMBOL_TIMEOUT                         0         // symbols
#define LORA_PUBLIC_NETWORK_ON                      true
#define LORA_FIX_LENGTH_PAYLOAD_ON                  false

// rgb red:   tx active
// rgb green: rx complete
#ifndef LoraWan_RGB
#define LoraWan_RGB 0
#endif

void on_tx_done(void);
void on_tx_timeout(void);
void on_rx_done(uint8_t* payload, uint16_t payload_len, int16_t rssi, int8_t snr);
void sleep(void);
void flash_rgb(void);
void gpio_on();

typedef enum {
    LOWPOWER,
    RX,
    TX
} States_t;

static RadioEvents_t radio_events;

uint32_t current_freq;
int16_t tx_count;
int16_t rx_count;
int16_t last_rssi;
States_t state;
bool sleep_active;
char tx_buf[256];

void setup() {
    tx_count = 0;
    rx_count = 0;
    last_rssi = 0;
    sleep_active = false;

    Serial.begin(115200);
    while (!Serial) delay(1);
    delay(100);
    Serial.println("CubeCell ASR6502 / SX1262");

    pinMode(Vext, OUTPUT);
    digitalWrite(Vext, LOW);

    delay(100);
    flash_rgb();

    pinMode(P3_3, INPUT);
    attachInterrupt(P3_3, sleep, FALLING);

    radio_events.TxDone = on_tx_done;
    radio_events.TxTimeout = on_tx_timeout;
    radio_events.RxDone = on_rx_done;
    radio_events.RxTimeout = on_rx_timeout;
    radio_events.RxError = on_rx_error;

    Radio.Init(&radio_events);
    Radio.SetPublicNetwork(LORA_PUBLIC_NETWORK_ON);
    Radio.SetChannel(RF_FREQUENCY1); current_freq = RF_FREQUENCY1;

    Radio.SetTxConfig(MODEM_LORA,   // RadioModems_t modem
        TX_OUTPUT_POWER,            // int8_t power
        0,                          // uint32_t fdev
        LORA_BANDWIDTH,             // uint32_t bandwidth
        LORA_SPREADING_FACTOR,      // uint32_t datarate
        LORA_CODINGRATE,            // uint8_t coderate
        LORA_PREAMBLE_LENGTH,       // uint16_t preambleLen
        LORA_FIX_LENGTH_PAYLOAD_ON, // bool fixLen
        true,                       // bool crcOn
        0, 0,                       // bool freqHopOn, uint8_t hopPeriod
        false,                      // bool iqInverted
        3000);                      // uint32_t timeout

    Radio.SetRxConfig(MODEM_LORA,   // RadioModems_t modem
        LORA_BANDWIDTH,             // uint32_t bandwidth
        LORA_SPREADING_FACTOR,      // uint32_t datarate
        LORA_CODINGRATE,            // uint8_t coderate
        0,                          // uint32_t bandwidthAfc
        LORA_PREAMBLE_LENGTH,       // uint16_t preambleLen
        LORA_SYMBOL_TIMEOUT,        // uint16_t symbTimeout
        LORA_FIX_LENGTH_PAYLOAD_ON, // bool fixLen
        0,                          // uint8_t payloadLen
        true,                       // bool crcOn
        0, 0,                       // bool freqHopOn, uint8_t hopPeriod
        false,                      // bool iqInverted
        true);                      // bool rxContinuous

    state = RX;
}

void flip_radio() {
    bool invertIq = false;
    if(current_freq != RF_FREQUENCY1) {
        current_freq = RF_FREQUENCY1;
        invertIq = false;
        Serial.printf("============================================================\n");
    } else {
        current_freq = RF_FREQUENCY2;
        invertIq = true;
    }

    Radio.SetChannel(current_freq);
    Radio.SetRxConfig(MODEM_LORA,   // RadioModems_t modem
        LORA_BANDWIDTH,             // uint32_t bandwidth
        LORA_SPREADING_FACTOR,      // uint32_t datarate
        LORA_CODINGRATE,            // uint8_t coderate
        0,                          // uint32_t bandwidthAfc
        LORA_PREAMBLE_LENGTH,       // uint16_t preambleLen
        LORA_SYMBOL_TIMEOUT,        // uint16_t symbTimeout
        LORA_FIX_LENGTH_PAYLOAD_ON, // bool fixLen
        0,                          // uint8_t payloadLen
        true,                       // bool crcOn
        0, 0,                       // bool freqHopOn, uint8_t hopPeriod
        invertIq,                   // bool iqInverted
        true);                      // bool rxContinuous
}

void loop() {
    switch (state) {
    case TX:
        delay(800);
        sprintf(tx_buf, "hello %d (rssi: %d)", ++tx_count, last_rssi);
        turnOnRGB(0x100000, 0);
        Serial.printf("\r\nsending packet \"%s\", length %d\r\n", tx_buf, strlen(tx_buf));
        Radio.Send((uint8_t*)tx_buf, strlen(tx_buf));
        state = LOWPOWER;
        break;

    case RX:
        Serial.printf("\r\nRX waiting for packet (%.1f MHz)...\n", (current_freq / 1000000.0));
        Radio.Rx(0);
        state = LOWPOWER;
        break;

/*    case LOWPOWER:
        if(sleep_active) {
            Radio.Sleep();
            Wire.end();
            detachInterrupt(RADIO_DIO_1);
            turnOffRGB();
            pinMode(GPIO0, ANALOG);
            pinMode(GPIO1, ANALOG);
            pinMode(GPIO2, ANALOG);
            pinMode(GPIO3, ANALOG);
            pinMode(GPIO4, ANALOG);
            pinMode(GPIO5, ANALOG);
            pinMode(Vext,  ANALOG);
            pinMode(ADC,   ANALOG);
        }

        lowPowerHandler();
        break;
*/
    default:
        break;
    }

    Radio.IrqProcess();
}

void on_tx_done(void) {
    Serial.println("TX complete...");
    turnOnRGB(0, 0);
    state = RX;
}

void on_tx_timeout(void) {
    Radio.Sleep();
    Serial.println("TX timeout...");
    //state = TX;
}

void on_rx_done(uint8_t* payload, uint16_t payload_len, int16_t rssi, int8_t snr) {
    last_rssi = rssi;

    gpio_on();
    turnOnRGB(0x001000, 100);
    turnOnRGB(0, 0);

    Radio.Sleep();

    Serial.printf("------------------------\r\n");
    Serial.printf("[RX packet %d] rssi: %d, snr: %d, size: %d\r\n", ++rx_count, rssi, snr, payload_len);
    print_mini_report(payload, payload_len, USE_COLOR_LOGGING);

    if(get_msg_type(payload) != mtype_unconfirmed_data_up) flip_radio();
    state = RX;
}

void on_rx_timeout(void) {
    Serial.println("** on_rx_timeout **");
}

void on_rx_error(void) {
    Serial.println("!! on_rx_error !!");
}

void sleep(void) {
    delay(10);
    if(digitalRead(P3_3) == 0) {
        Serial.println("click");
        sleep_active = true;
    }
}

void flash_rgb(void) {
    for(uint32_t i = 0; i <= 30; i++) {
        turnOnRGB(i << 16, 10);
    }
    for(uint32_t i = 0; i <= 30; i++) {
        turnOnRGB(i << 8, 10);
    }
    for(uint32_t i = 0; i <= 30; i++) {
        turnOnRGB(i, 10);
    }
    turnOnRGB(0, 0);
}

void gpio_on(void) {
    pinMode(GPIO0, OUTPUT);
    pinMode(GPIO1, OUTPUT);
    pinMode(GPIO2, OUTPUT);
    pinMode(GPIO3, OUTPUT);
    pinMode(GPIO4, OUTPUT);
    pinMode(GPIO5, OUTPUT);
    digitalWrite(GPIO0, HIGH);
    digitalWrite(GPIO1, HIGH);
    digitalWrite(GPIO2, HIGH);
    digitalWrite(GPIO3, HIGH);
    digitalWrite(GPIO4, HIGH);
    digitalWrite(GPIO5, HIGH);
}

