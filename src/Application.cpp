#include <Arduino.h>
#include <driver/i2s.h>

#include "Application.h"
#include "I2SMEMSSampler.h"
#include "I2SOutput.h"
#include "UdpTransport.h"

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
    .channel_format = I2S_MIC_CHANNEL,
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

// Task to write samples to our server
void samples_task(void *param)
{
  Application *application = (Application *)param;
  const TickType_t xMaxBlockTime = pdMS_TO_TICKS(100);
  while (true)
  {
    // wait for some samples to save
    uint32_t ulNotificationValue = ulTaskNotifyTake(pdTRUE, xMaxBlockTime);
    if (ulNotificationValue > 0)
    {
      application->service();
    }
  }
}

Application::Application()
{
  output = new I2SOutput();
  input = new I2SMEMSSampler(i2s_mic_pins);
  transport = new UdpTransport();
}

void Application::begin()
{
  TaskHandle_t samples_task_handle;
  xTaskCreatePinnedToCore(samples_task, "I2S Writer Task", 4096, this, 1, &samples_task_handle, 1);

  Serial.println("Starting I2S Output");
  output->start(I2S_NUM_1, i2s_speaker_pins);

  Serial.println("Starting I2S Input");
  input->start(I2S_NUM_0, i2sMemsConfig, 1436 / 2, samples_task_handle);

  transport->begin(output);
}

void Application::service()
{
  transport->send_audio(reinterpret_cast<uint8_t *>(input->getCapturedAudioBuffer()), input->getBufferSizeInBytes());
}