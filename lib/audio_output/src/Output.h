#pragma once

#include <freertos/FreeRTOS.h>
#include <driver/i2s.h>

/**
 * Base Class for both the DAC and I2S output
 **/
class Output
{
protected:
  i2s_port_t m_i2s_port = I2S_NUM_0;

public:
  Output(i2s_port_t i2s_port);
  virtual void start(int sample_rate) = 0;
  void stop();
  // override this in derived classes to turn the sample into
  // something the output device expects - for the default case
  // this is simply a pass through
  virtual int16_t process_sample(int16_t sample) { return sample; }
  void write(int16_t *samples, int count);
};
