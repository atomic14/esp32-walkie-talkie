#include <Arduino.h>
#include "driver/i2s.h"
#include "driver/adc.h"
#include "ADCSampler.h"

ADCSampler::ADCSampler(adc_unit_t adcUnit, adc1_channel_t adcChannel)
{
    m_adcUnit = adcUnit;
    m_adcChannel = adcChannel;
}

void ADCSampler::configureI2S()
{
    //init ADC pad
    i2s_set_adc_mode(m_adcUnit, m_adcChannel);
    // enable the adc
    i2s_adc_enable(getI2SPort());
}

/**
 * Process the raw data that have been read from the I2S peripherals into samples
 **/
void ADCSampler::processI2SData(uint8_t *i2sData, size_t bytesRead)
{
    uint16_t *rawSamples = (uint16_t *)i2sData;
    for (int i = 0; i < bytesRead / 2; i++)
    {
        addSample((2048 - (rawSamples[i] & 0xfff)) * 15);
    }
}
