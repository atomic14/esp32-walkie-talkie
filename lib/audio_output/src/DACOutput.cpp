#include "DACOutput.h"

#if CONFIG_IDF_TARGET_ESP32

void DACOutput::start(uint32_t sample_rate)
{
    // i2s config for writing both channels of I2S
    i2s_config_t i2s_config = {
        .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_TX | I2S_MODE_DAC_BUILT_IN),
#if ESP_IDF_VERSION >= ESP_IDF_VERSION_VAL(4, 4, 1)
        .sample_rate = sample_rate, 
#else
        .sample_rate = (int)sample_rate, 
#endif
        .bits_per_sample = I2S_BITS_PER_SAMPLE_16BIT,
        .channel_format = I2S_CHANNEL_FMT_RIGHT_LEFT,
#if ESP_IDF_VERSION >= ESP_IDF_VERSION_VAL(4, 2, 0)
        .communication_format = (i2s_comm_format_t)(I2S_COMM_FORMAT_STAND_MSB),
#else
        .communication_format = (i2s_comm_format_t)(I2S_COMM_FORMAT_I2S_MSB),        
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
    i2s_driver_install(I2S_NUM_0, &i2s_config, 0, NULL);
    // enable the DAC channels
    i2s_set_dac_mode(I2S_DAC_CHANNEL_BOTH_EN);
    // clear the DMA buffers
    i2s_zero_dma_buffer(I2S_NUM_0);

    i2s_start(I2S_NUM_0);
}

#endif
