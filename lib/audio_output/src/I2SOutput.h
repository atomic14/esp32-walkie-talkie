#pragma once

#include "Output.h"

/**
 * Base Class for both the ADC and I2S sampler
 **/
class I2SOutput : public Output
{
private:
    i2s_pin_config_t m_i2s_pins;

public:
    I2SOutput(i2s_port_t i2s_port, i2s_pin_config_t &i2s_pins);
    void start(uint32_t sample_rate);
};
