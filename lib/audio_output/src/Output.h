#ifndef __output_h__
#define __output_h__

#include <Arduino.h>
#include <driver/i2s.h>

class OutputBuffer;
/**
 * Base Class for both the DAC and I2S output
 **/
class Output
{
protected:
  i2s_port_t m_i2s_port = I2S_NUM_0;
  // I2S write task
  TaskHandle_t m_i2s_sample_task_handle = NULL;
  // i2s writer queue
  QueueHandle_t m_i2s_queue = NULL;
  // output buffer
  OutputBuffer *m_output_buffer = NULL;

  void start(i2s_port_t i2sPort, OutputBuffer *output_buffer);

public:
  void stop();
  // override this in derived classes to turn the 8 bit PCM sample into
  // something the output device expects
  virtual int16_t process_sample(uint8_t sample) = 0;
  // give the task function access to our data
  friend void i2s_writer_task(void *param);
};

#endif