
#include "I2SOutput.h"

I2SOutput::I2SOutput(i2s_port_t i2s_port, i2s_pin_config_t &i2s_pins) : Output(i2s_port), m_i2s_pins(i2s_pins)
{
}

void I2SOutput::start(uint32_t sample_rate)
{
    // i2s config for writing both channels of I2S
    i2s_config_t i2s_config = {
        .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_TX),
#if ESP_IDF_VERSION >= ESP_IDF_VERSION_VAL(4, 4, 1)
        .sample_rate = sample_rate, 
#else
        .sample_rate = (int)sample_rate, 
#endif
        .bits_per_sample = I2S_BITS_PER_SAMPLE_16BIT,
        .channel_format = I2S_CHANNEL_FMT_RIGHT_LEFT,
#if ESP_IDF_VERSION >= ESP_IDF_VERSION_VAL(4, 2, 0)
        .communication_format = (i2s_comm_format_t)(I2S_COMM_FORMAT_STAND_I2S),
#else
        .communication_format = (i2s_comm_format_t)(I2S_COMM_FORMAT_I2S),
#endif
        .intr_alloc_flags = ESP_INTR_FLAG_LEVEL1,
        .dma_buf_count = 2,
        .dma_buf_len = 1024,
        .use_apll = false,
        .tx_desc_auto_clear = true,
        .fixed_mclk = 0,
#if ESP_IDF_VERSION >= ESP_IDF_VERSION_VAL(4, 4, 1)
        .mclk_multiple = I2S_MCLK_MULTIPLE_DEFAULT, // Unused
        .bits_per_chan = I2S_BITS_PER_CHAN_DEFAULT // Use bits per sample
#endif
        };
    //install and start i2s driver
    i2s_driver_install(m_i2s_port, &i2s_config, 0, NULL);
    // set up the i2s pins
    i2s_set_pin(m_i2s_port, &m_i2s_pins);
    // clear the DMA buffers
    i2s_zero_dma_buffer(m_i2s_port);

    i2s_start(m_i2s_port);
}
