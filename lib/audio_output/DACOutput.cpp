
#include <Arduino.h>
#include "driver/i2s.h"
#include <math.h>
#include <SPIFFS.h>
#include <FS.h>

#include "SampleSource.h"
#include "DACOutput.h"

// number of frames to try and send at once (a frame is a left and right sample)
#define NUM_FRAMES_TO_SEND 128

void i2sWriterTask(void *param)
{
    DACOutput *output = (DACOutput *)param;
    int availableBytes = 0;
    int buffer_position = 0;
    Frame_t frames[128];
    while (true)
    {
        // wait for some data to be requested
        i2s_event_t evt;
        if (xQueueReceive(output->m_i2sQueue, &evt, portMAX_DELAY) == pdPASS)
        {
            if (evt.type == I2S_EVENT_TX_DONE)
            {
                size_t bytesWritten = 0;
                do
                {
                    if (availableBytes == 0)
                    {
                        // get some frames from the wave file - a frame consists of a 16 bit left and right sample
                        output->m_sample_generator->getFrames(frames, NUM_FRAMES_TO_SEND);
                        // how many bytes do we now have to send
                        availableBytes = NUM_FRAMES_TO_SEND * sizeof(uint32_t);
                        // reset the buffer position back to the start
                        buffer_position = 0;
                    }
                    // do we have something to write?
                    if (availableBytes > 0)
                    {
                        // write data to the i2s peripheral
                        i2s_write(I2S_NUM_0, buffer_position + (uint8_t *)frames,
                                  availableBytes, &bytesWritten, portMAX_DELAY);
                        availableBytes -= bytesWritten;
                        buffer_position += bytesWritten;
                    }
                } while (bytesWritten > 0);
            }
        }
    }
}

void DACOutput::start(SampleSource *sample_generator)
{
    m_sample_generator = sample_generator;
    // i2s config for writing both channels of I2S
    i2s_config_t i2sConfig = {
        .mode = (i2s_mode_t)(I2S_MODE_MASTER | I2S_MODE_TX | I2S_MODE_DAC_BUILT_IN),
        .sample_rate = 16000,
        .bits_per_sample = I2S_BITS_PER_SAMPLE_16BIT,
        .channel_format = I2S_CHANNEL_FMT_RIGHT_LEFT,
        .communication_format = (i2s_comm_format_t)(I2S_COMM_FORMAT_I2S_MSB),
        .intr_alloc_flags = ESP_INTR_FLAG_LEVEL1,
        .dma_buf_count = 4,
        .dma_buf_len = 64};

    //install and start i2s driver
    i2s_driver_install(I2S_NUM_0, &i2sConfig, 4, &m_i2sQueue);
    // enable the DAC channels
    i2s_set_dac_mode(I2S_DAC_CHANNEL_BOTH_EN);
    // clear the DMA buffers
    i2s_zero_dma_buffer(I2S_NUM_0);
    // start a task to write samples to the i2s peripheral
    TaskHandle_t writerTaskHandle;
    xTaskCreate(i2sWriterTask, "i2s Writer Task", 4096, this, 1, &writerTaskHandle);
}
