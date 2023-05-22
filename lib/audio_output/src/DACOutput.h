#pragma once
#include <freertos/FreeRTOS.h>
#include <driver/i2s.h>

#include "Output.h"

/**
 * Base Class for both the ADC and I2S sampler
 **/
class DACOutput : public Output
{
public:
    DACOutput(i2s_port_t i2s_port) : Output(i2s_port) {}
    void start(uint32_t sample_rate);
    virtual int16_t process_sample(int16_t sample)
    {
        // DAC needs unsigned 16 bit samples
        return sample + 32768;
    }
};
