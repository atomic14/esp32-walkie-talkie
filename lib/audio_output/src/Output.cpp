
#include <Arduino.h>
#include <FreeRTOS.h>
#include <driver/i2s.h>
#include "Output.h"
#include "OutputBuffer.h"

// number of frames to try and send at once (a frame is a left and right sample)
const int NUM_FRAMES_TO_SEND = 256;

void i2s_writer_task(void *param)
{
  Output *output = (Output *)param;
  int available_bytes = 0;
  int buffer_position = 0;
  // the raw samples
  uint8_t *raw_samples = (uint8_t *)malloc(NUM_FRAMES_TO_SEND);
  // this will contained the prepared samples for sending to the I2S device
  int16_t *frames = (int16_t *)malloc(2 * sizeof(int16_t) * NUM_FRAMES_TO_SEND);
  while (true)
  {
    // wait for some data to be requested
    i2s_event_t evt;
    if (xQueueReceive(output->m_i2s_queue, &evt, portMAX_DELAY) == pdPASS)
    {
      if (evt.type == I2S_EVENT_TX_DONE)
      {
        size_t bytes_written = 0;
        do
        {
          if (available_bytes == 0)
          {
            // fill the output buffer with more data
            output->m_output_buffer->remove_samples(raw_samples, NUM_FRAMES_TO_SEND);
            for (int i = 0; i < NUM_FRAMES_TO_SEND; i++)
            {
              int16_t sample = output->process_sample(raw_samples[i]);
              frames[i * 2] = sample;
              frames[i * 2 + 1] = sample;
            }
            available_bytes = NUM_FRAMES_TO_SEND * 2 * sizeof(uint16_t);
            buffer_position = 0;
          }
          // write data to the i2s peripheral
          i2s_write(output->m_i2s_port, buffer_position + (uint8_t *)frames,
                    available_bytes, &bytes_written, portMAX_DELAY);
          available_bytes -= bytes_written;
          buffer_position += bytes_written;
        } while (bytes_written > 0);
      }
    }
  }
}

void Output::start(i2s_port_t i2s_port, OutputBuffer *output_buffer)
{
  m_i2s_port = i2s_port;
  m_output_buffer = output_buffer;
  // start a task to write samples to the i2s peripheral
  if (m_i2s_sample_task_handle == NULL)
  {
    xTaskCreate(i2s_writer_task, "i2s Writer Task", 4096, this, 1, &m_i2s_sample_task_handle);
  }
}

void Output::stop()
{
  // stop the i2S driver
  i2s_driver_uninstall(m_i2s_port);
  // NOTE this leaves the task running - there's not really a clean way of terminating tasks
}