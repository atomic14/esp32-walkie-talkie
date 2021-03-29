
#include <Arduino.h>
#include "I2SSampler.h"
#include "driver/i2s.h"
#include "SampleSink.h"

void i2s_reader_task(void *param)
{
    I2SSampler *sampler = (I2SSampler *)param;
    while (true)
    {
        // wait for some data to arrive on the queue
        i2s_event_t evt;
        if (xQueueReceive(sampler->m_i2sQueue, &evt, portMAX_DELAY) == pdPASS)
        {
            if (evt.type == I2S_EVENT_RX_DONE)
            {
                size_t bytesRead = 0;
                do
                {
                    // read data from the I2S peripheral
                    uint8_t i2sData[1024];
                    // read from i2s
                    i2s_read(sampler->getI2SPort(), i2sData, 1024, &bytesRead, 10);
                    // process the raw data
                    sampler->processI2SData(i2sData, bytesRead);
                } while (bytesRead > 0);
            }
        }
    }
}

void I2SSampler::start(i2s_port_t i2sPort, i2s_config_t &i2sConfig, SampleSink *sample_sink)
{
    m_i2sPort = i2sPort;
    m_sample_sink = sample_sink;

    //install and start i2s driver
    i2s_driver_install(m_i2sPort, &i2sConfig, 4, &m_i2sQueue);
    // set up the I2S configuration from the subclass
    configureI2S();
    if (!m_i2s_reader_task_handle)
    {
        // start a task to read samples from the ADC
        xTaskCreatePinnedToCore(i2s_reader_task, "i2s Reader Task", 4096, this, 1, &m_i2s_reader_task_handle, 1);
    }
}

void I2SSampler::stop()
{
    // stop the i2S driver
    i2s_driver_uninstall(m_i2sPort);
    // NOTE this leaves the task running - there's not really a clean way of terminating tasks
}
