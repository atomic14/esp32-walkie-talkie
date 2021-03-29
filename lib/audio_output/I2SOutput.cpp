
#include <Arduino.h>
#include <FreeRTOS.h>
#include <driver/i2s.h>
#include "I2SOutput.h"

// number of frames to try and send at once (a frame is a left and right sample)
const int NUM_FRAMES_TO_SEND = 128;
// number of ms to build up before sending samples to the amplifier - assumes 16KHz sample rate
const int SAMPLES_TO_BUFFER = 300 * 16;

void i2s_writer_task(void *param)
{
    I2SOutput *output = (I2SOutput *)param;
    int available_bytes = 0;
    int buffer_position = 0;
    int16_t frames[2 * NUM_FRAMES_TO_SEND];
    while (true)
    {
        // wait for some data to be requested
        i2s_event_t evt;
        if (xQueueReceive(output->m_i2s_queue, &evt, portMAX_DELAY) == pdPASS)
        {
            if (evt.type == I2S_EVENT_TX_DONE)
            {
                size_t bytes_written = 0;
                do
                {
                    if (available_bytes == 0)
                    {
                        // pull more samples from the circular buffer to send
                        // we convert them to stereo for the AMP board that will
                        // output (L+R)/2
                        for (int i = 0; i < NUM_FRAMES_TO_SEND; i++)
                        {
                            int16_t sample = output->pop_sample();
                            frames[i * 2] = sample;
                            frames[i * 2 + 1] = sample;
                        }
                        available_bytes = NUM_FRAMES_TO_SEND * 2 * sizeof(int16_t);
                    }
                    // write data to the i2s peripheral
                    i2s_write(output->m_i2s_port, buffer_position + (uint8_t *)frames,
                              available_bytes, &bytes_written, portMAX_DELAY);
                    available_bytes -= bytes_written;
                    buffer_position += bytes_written;
                } while (bytes_written > 0);
            }
        }
    }
}

I2SOutput::I2SOutput()
{
    // reading and writing to the beginning of the buffer
    m_read_head = 0;
    m_write_head = 0;
    // make sufficient space for our bufferring and incoming data
    m_buffer_size = 2 * SAMPLES_TO_BUFFER;
    m_buffer = (int16_t *)malloc(sizeof(int16_t) * m_buffer_size);
}

void I2SOutput::push_samples(int16_t *samples, int count)
{
    // copy the samples into the buffer wrapping around as needed
    for (int i = 0; i < count; i++)
    {
        m_buffer[m_write_head] = samples[count];
        m_write_head = (m_write_head + 1) % m_buffer_size;
    }
}

int16_t I2SOutput::pop_sample()
{
    // have we buffered enough data?
    if (m_read_head == (m_write_head + SAMPLES_TO_BUFFER) % m_buffer_size)
    {
        // don't have enough data buffered so just write 0s
        return 0;
    }
    int16_t sample = m_buffer[m_read_head];
    m_read_head = (m_read_head + 1) % m_buffer_size;
    return sample;
}

void I2SOutput::start(i2s_port_t i2s_port, i2s_pin_config_t &i2s_pins)
{
    m_i2s_port = i2s_port;
    m_i2s_pins = i2s_pins;

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
    i2s_driver_install(m_i2s_port, &i2s_config, 4, &m_i2s_queue);
    // set up the i2s pins
    i2s_set_pin(m_i2s_port, &m_i2s_pins);
    // clear the DMA buffers
    i2s_zero_dma_buffer(m_i2s_port);
    // start a task to write samples to the i2s peripheral
    TaskHandle_t writer_task_handle;
    xTaskCreate(i2s_writer_task, "i2s Writer Task", 4096, this, 1, &writer_task_handle);
}
