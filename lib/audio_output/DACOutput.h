#ifndef __sampler_base_h__
#define __sampler_base_h__

#include <Arduino.h>
#include "driver/i2s.h"

class SampleSource;

/**
 * Base Class for both the ADC and I2S sampler
 **/
class DACOutput
{
private:
    // I2S write task
    TaskHandle_t m_i2sWriterTaskHandle;
    // i2s writer queue
    QueueHandle_t m_i2sQueue;
    // src of samples for us to play
    SampleSource *m_sample_generator;

public:
    void start(SampleSource *sample_generator);

    friend void i2sWriterTask(void *param);
};

#endif