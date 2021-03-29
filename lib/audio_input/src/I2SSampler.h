#ifndef __sampler_base_h__
#define __sampler_base_h__

#include <Arduino.h>
#include <driver/i2s.h>

class SampleSink;
/**
 * Base Class for both the ADC and I2S sampler
 **/
class I2SSampler
{
private:
    // I2S reader task
    TaskHandle_t m_i2s_reader_task_handle = NULL;
    // i2s reader queue
    QueueHandle_t m_i2sQueue = NULL;
    // i2s port
    i2s_port_t m_i2sPort = I2S_NUM_0;

protected:
    // the input buffer to store samples in
    SampleSink *m_sample_sink;

    virtual void configureI2S() = 0;
    virtual void processI2SData(uint8_t *i2sData, size_t bytesRead) = 0;
    i2s_port_t getI2SPort()
    {
        return m_i2sPort;
    }

public:
    void start(i2s_port_t i2sPort, i2s_config_t &i2sConfig, SampleSink *sample_sink);
    void stop();
    friend void i2s_reader_task(void *param);
};

#endif