#ifndef __i2s_output_h__
#define __i2s_output_h__

#include <Arduino.h>
#include <driver/i2s.h>
#include "Output.h"

class OutputBuffer;
/**
 * Base Class for both the ADC and I2S sampler
 **/
class I2SOutput : public Output
{
public:
    void start(i2s_port_t i2sPort, i2s_pin_config_t &i2sPins, OutputBuffer *output_buffer);
    virtual int16_t process_sample(uint8_t sample)
    {
        // I2S needs signed 16 bit samples
        int16_t processed = sample;
        return (processed - 128) << 5;
    }
};

#endif