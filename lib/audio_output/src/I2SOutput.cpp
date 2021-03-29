
#include <Arduino.h>
#include <FreeRTOS.h>
#include <driver/i2s.h>
#include "I2SOutput.h"

void I2SOutput::start(i2s_port_t i2s_port, i2s_pin_config_t &i2s_pins, OutputBuffer *output_buffer)
{
    // i2s config for writing both channels of I2S
    i2s_config_t i2s_config = {
        .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_TX),
        .sample_rate = 16000,
        .bits_per_sample = I2S_BITS_PER_SAMPLE_16BIT,
        .channel_format = I2S_CHANNEL_FMT_RIGHT_LEFT,
        .communication_format = (i2s_comm_format_t)(I2S_COMM_FORMAT_I2S),
        .intr_alloc_flags = ESP_INTR_FLAG_LEVEL1,
        .dma_buf_count = 2,
        .dma_buf_len = 64,
        .use_apll = false,
        .tx_desc_auto_clear = true,
        .fixed_mclk = 0};
    //install and start i2s driver
    i2s_driver_install(i2s_port, &i2s_config, 4, &m_i2s_queue);
    // set up the i2s pins
    i2s_set_pin(i2s_port, &i2s_pins);
    // clear the DMA buffers
    i2s_zero_dma_buffer(i2s_port);
    // start sending data
    Output::start(i2s_port, output_buffer);
}
