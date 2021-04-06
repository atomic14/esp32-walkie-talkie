#include <Arduino.h>
#include <driver/i2s.h>
#include <WiFi.h>

#include "Application.h"
#include "I2SMEMSSampler.h"
#include "ADCSampler.h"
#include "I2SOutput.h"
#include "UdpTransport.h"
#include "EspNowTransport.h"
#include "OutputBuffer.h"
#include "SampleSink.h"
#include "TinyPICOIndicatorLed.h"
#include "config.h"

// i2s config for using the internal ADC
i2s_config_t adcI2SConfig = {
    .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_RX | I2S_MODE_ADC_BUILT_IN),
    .sample_rate = 16000,
    .bits_per_sample = I2S_BITS_PER_SAMPLE_16BIT,
    .channel_format = I2S_CHANNEL_FMT_ONLY_LEFT,
    .communication_format = I2S_COMM_FORMAT_I2S_LSB,
    .intr_alloc_flags = ESP_INTR_FLAG_LEVEL1,
    .dma_buf_count = 4,
    .dma_buf_len = 64,
    .use_apll = false,
    .tx_desc_auto_clear = false,
    .fixed_mclk = 0};

// i2s config for reading from both channels of I2S
i2s_config_t i2sMemsConfig = {
    .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_RX),
    .sample_rate = 16000,
    .bits_per_sample = I2S_BITS_PER_SAMPLE_32BIT,
    .channel_format = I2S_CHANNEL_FMT_ONLY_LEFT,
    .communication_format = i2s_comm_format_t(I2S_COMM_FORMAT_I2S),
    .intr_alloc_flags = ESP_INTR_FLAG_LEVEL1,
    .dma_buf_count = 4,
    .dma_buf_len = 64,
    .use_apll = false,
    .tx_desc_auto_clear = false,
    .fixed_mclk = 0};

// i2s microphone pins
i2s_pin_config_t i2s_mic_pins = {
    .bck_io_num = I2S_MIC_SERIAL_CLOCK,
    .ws_io_num = I2S_MIC_LEFT_RIGHT_CLOCK,
    .data_out_num = I2S_PIN_NO_CHANGE,
    .data_in_num = I2S_MIC_SERIAL_DATA};

// i2s speaker pins
i2s_pin_config_t i2s_speaker_pins = {
    .bck_io_num = I2S_SPEAKER_SERIAL_CLOCK,
    .ws_io_num = I2S_SPEAKER_LEFT_RIGHT_CLOCK,
    .data_out_num = I2S_SPEAKER_SERIAL_DATA,
    .data_in_num = I2S_PIN_NO_CHANGE};

Application::Application()
{
  m_output_buffer = new OutputBuffer(300 * 16);
#ifdef USE_I2S_MIC_INPUT
  m_input = new I2SMEMSSampler(i2s_mic_pins);
#else
  m_input = new ADCSampler(ADC_UNIT_1, ADC1_CHANNEL_7);
#endif
  m_output = new I2SOutput();
#ifdef USE_ESP_NOW
  m_transport = new EspNowTransport(m_output_buffer);
#else
  m_transport = new UdpTransport(m_output_buffer);
#endif
  m_indicator_led = new TinyPICOIndicatorLed();
}

void Application::begin()
{
  // show a flashing indicator that we are trying to connect
  m_indicator_led->set_default_color(0);
  m_indicator_led->set_is_flashing(true, 0xff0000);
  m_indicator_led->begin();
  // bring up WiFi
  WiFi.mode(WIFI_STA);
  // but don't connect if we're using ESP NOW
#ifndef USE_ESP_NOW
  WiFi.begin(WIFI_SSID, WIFI_PSWD);
  if (WiFi.waitForConnectResult() != WL_CONNECTED)
  {
    Serial.println("Connection Failed! Rebooting...");
    delay(5000);
    ESP.restart();
  }
  // this has a dramatic effect on packet RTT
  WiFi.setSleep(WIFI_PS_NONE);
  Serial.print("My IP Address is: ");
  Serial.println(WiFi.localIP());
#else
  WiFi.disconnect();
#endif
  Serial.print("My MAC Address is: ");
  Serial.println(WiFi.macAddress());
  // do any setup of the transport
  m_transport->begin();
  // connected so show a solid green light
  m_indicator_led->set_default_color(0x00ff00);
  m_indicator_led->set_is_flashing(false, 0x00ff00);
  // setup the transmit button
  pinMode(GPIO_TRANSMIT_BUTTON, INPUT_PULLDOWN);
  // start both the input and output I2S devices
  Serial.println("Starting I2S Output");
  m_output->start(I2S_NUM_0, i2s_speaker_pins, m_output_buffer);
  Serial.println("Starting I2S Input");
  m_input->start(I2S_NUM_1, i2sMemsConfig, m_transport);
}

void Application::loop()
{
  unsigned long current_time = millis();
  // run the application state machine
  bool transmit_pushed = digitalRead(GPIO_TRANSMIT_BUTTON);
  switch (m_current_state)
  {
  // the application is waiting for the user to push the transmit button or for audio packets to arrive
  case IDLE:
  {
    // have we started receiving packets?
    if (current_time - m_transport->get_last_packet_received() < 200)
    {
      m_indicator_led->set_is_flashing(true, 0xff0000);
      m_current_state = RECEIVING;
    }
    else
    {
      if (transmit_pushed)
      {
        m_indicator_led->set_is_flashing(true, 0xff0000);
        m_current_state = TRANSMITTING;
        m_transport->should_send(true);
      }
    }
  }
  break;
  // the user has pushed the transmit button
  case TRANSMITTING:
  {
    if (!transmit_pushed)
    {
      m_indicator_led->set_is_flashing(false, 0);
      m_current_state = IDLE;
      m_transport->should_send(false);
      Serial.println("Finished transmitting");
    }
  }
  break;
  // we are receiving audio packets
  case RECEIVING:
  {
    // have we stopped receiving packets?
    if (current_time - m_transport->get_last_packet_received() > 200)
    {
      // haven't received any packets - switch back to idle
      m_indicator_led->set_is_flashing(false, 0);
      m_current_state = IDLE;
      Serial.println("Finished receiving");
    }
  }
  break;
  }
}
