#pragma once

#include <Arduino.h>
#include <freertos/FreeRTOS.h>

/**
 * @brief Circular buffer for 8 bit unsigned PCM samples
 * 
 */
class OutputBuffer
{
private:
  // how many samples should we buffer before outputting data?
  int m_number_samples_to_buffer;
  // where are we reading from
  int m_read_head;
  // where are we writing to
  int m_write_head;
  // keep track of how many samples we have
  int m_available_samples;
  // the total size of the buffer
  int m_buffer_size;
  // are we currently buffering samples?
  bool m_buffering;
  // the sample buffer
  uint8_t *m_buffer;
  // thread safety
  SemaphoreHandle_t m_semaphore;

public:
  OutputBuffer(int number_samples_to_buffer) : m_number_samples_to_buffer(number_samples_to_buffer)
  {
    // create a semaphore and make it available for locking
    m_semaphore = xSemaphoreCreateBinary();
    xSemaphoreGive(m_semaphore);
    // set reading and writing to the beginning of the buffer
    m_read_head = 0;
    m_write_head = 0;
    m_available_samples = 0;
    // we'll start off buffering data as we have no samples yet
    m_buffering = true;
    // make sufficient space for the bufferring and incoming data
    m_buffer_size = 3 * number_samples_to_buffer;
    m_buffer = (uint8_t *)malloc(m_buffer_size);
    memset(m_buffer, 0, m_buffer_size);
    if (!m_buffer)
    {
      Serial.println("Failed to allocate buffer");
    }
  }

  // we're adding samples that are 8 bit as they are coming from the transport
  void add_samples(const uint8_t *samples, int count)
  {
    xSemaphoreTake(m_semaphore, portMAX_DELAY);
    // check if there is still room in the buffer
    if (m_available_samples + count <= m_buffer_size)
    {
      // copy the samples into the buffer wrapping around as needed
      for (int i = 0; i < count; i++)
      {
        m_buffer[m_write_head] = samples[i];
        m_write_head = (m_write_head + 1) % m_buffer_size;
      }
      m_available_samples += count;
    }
    xSemaphoreGive(m_semaphore);
  }

  // convert the samples to 16 bit as they are going to the output
  void remove_samples(int16_t *samples, int count)
  {
    xSemaphoreTake(m_semaphore, portMAX_DELAY);
    for (int i = 0; i < count; i++)
    {
      samples[i] = 0;
      // if we have no samples and we aren't already buffering then we need to start buffering
      if (m_available_samples == 0 && !m_buffering)
      {
        Serial.println("Buffering");
        m_buffering = true;
        samples[i] = 0;
      }
      // are we buffering?
      if (m_buffering && m_available_samples < m_number_samples_to_buffer)
      {
        // just return 0 as we don't have enough samples yet
        samples[i] = 0;
      }
      else
      {
        // we've buffered enough samples so no need to buffer anymore
        m_buffering = false;
        // just send back the samples we've got and move the read head forward
        int16_t sample = m_buffer[m_read_head];
        samples[i] = (sample - 128) << 5;
        m_read_head = (m_read_head + 1) % m_buffer_size;
        m_available_samples--;
      }
    }
    xSemaphoreGive(m_semaphore);
  }

  void flush()
  {
    // flush all samples in the outputbuffer
    xSemaphoreTake(m_semaphore, portMAX_DELAY);
    m_read_head = 0;
    m_write_head = 0;
    m_available_samples = 0;
    xSemaphoreGive(m_semaphore);
  }
};
