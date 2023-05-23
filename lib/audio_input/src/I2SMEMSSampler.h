#pragma once

#include "I2SSampler.h"

class I2SMEMSSampler : public I2SSampler
{
private:
    i2s_pin_config_t m_i2sPins;
    bool m_fixSPH0645;
    int32_t *m_raw_samples;
    int m_raw_samples_size;

protected:
    void configureI2S();
    
public:
    I2SMEMSSampler(
        i2s_port_t i2s_port,
        i2s_pin_config_t &i2s_pins,
        i2s_config_t i2s_config,
        int raw_samples_size,
        bool fixSPH0645 = false);
    ~I2SMEMSSampler();
    virtual int read(int16_t *samples, int count);
};
