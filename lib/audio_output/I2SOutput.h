#ifndef __i2s_output_h__
#define __i2s_output_h__

#include <Arduino.h>
#include <driver/i2s.h>

/**
 * Base Class for both the ADC and I2S sampler
 **/
class I2SOutput
{
private:
    i2s_port_t m_i2s_port;
    i2s_pin_config_t m_i2s_pins;
    // I2S write task
    TaskHandle_t m_i2s_writer_task_handle;
    // i2s writer queue
    QueueHandle_t m_i2s_queue;
    // src of samples for us to play
    int m_read_head;
    int m_write_head;
    int m_buffer_size;
    int16_t *m_buffer;

public:
    I2SOutput();
    void start(i2s_port_t i2sPort, i2s_pin_config_t &i2sPins);
    void push_samples(int16_t *samples, int count);
    int16_t pop_sample();
    friend void i2s_writer_task(void *param);
};

#endif