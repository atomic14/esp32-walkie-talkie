#ifndef __dac_output_h__
#define __dac_output_h__

#include <Arduino.h>
#include "Output.h"
#include <driver/i2s.h>

class OutputBuffer;
/**
 * Base Class for both the ADC and I2S sampler
 **/
class DACOutput : public Output
{
public:
    void start(OutputBuffer *output_buffer);
    virtual int16_t process_sample(uint8_t sample)
    {
        // DAC needs unsigned 16 bit samples
        return sample << 8;
    }
};

#endif