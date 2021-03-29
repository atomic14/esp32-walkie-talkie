#ifndef __i2s_sampler_h__
#define __i2s_sampler_h__

#include "I2SSampler.h"
// #include "SampleFilter.h"

class I2SMEMSSampler : public I2SSampler
{
private:
    i2s_pin_config_t m_i2sPins;
    // SampleFilter sampleFilter;
    bool m_fixSPH0645;

protected:
    void configureI2S();
    void processI2SData(uint8_t *i2sData, size_t bytesRead);

public:
    I2SMEMSSampler(i2s_pin_config_t &i2sPins, bool fixSPH0645 = false);
};

#endif