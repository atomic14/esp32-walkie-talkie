#pragma once

#include <FreeRTOS.h>

/**
 * @brief Destination for audio data
 * 
 */
class SampleSink
{
public:
  virtual void add_sample(int16_t sample) = 0;
};

/**
 * @brief Double buffered samples - we can be processing one 
 * set of samples while the other samples are being read.
 * 
 */
class DoubleBuffer : public SampleSink
{
private:
  // double buffer so we can be capturing samples while sending data
  uint8_t *m_audioBuffer1;
  uint8_t *m_audioBuffer2;
  // current position in the audio buffer
  int32_t m_audioBufferPos = 0;
  // current audio buffer
  uint8_t *m_currentAudioBuffer;
  // buffer containing samples that have been captured already
  uint8_t *m_capturedAudioBuffer;
  // size of the audio buffer in samples
  int32_t m_bufferSize;
  // processing task handle - notified every time we reach the
  // end of a buffer
  TaskHandle_t m_processing_task_handle;

public:
  DoubleBuffer(int buffer_size, TaskHandle_t processing_task_handle = NULL)
  {
    m_processing_task_handle = processing_task_handle;
    m_bufferSize = buffer_size;
    m_audioBuffer1 = (uint8_t *)malloc(m_bufferSize);
    m_audioBuffer2 = (uint8_t *)malloc(m_bufferSize);

    m_currentAudioBuffer = m_audioBuffer1;
    m_capturedAudioBuffer = m_audioBuffer2;
  }
  int32_t get_buffer_size()
  {
    return m_bufferSize;
  };
  uint8_t *get_captured_audio_buffer()
  {
    return m_capturedAudioBuffer;
  }
  void add_sample(int16_t sample)
  {
    // add the sample to the current audio buffer converting to unsigned 8 bit PCM
    m_currentAudioBuffer[m_audioBufferPos] = (sample + 32768) >> 8;
    m_audioBufferPos++;
    // have we filled the buffer with data?
    if (m_audioBufferPos == m_bufferSize)
    {
      // swap to the other buffer
      std::swap(m_currentAudioBuffer, m_capturedAudioBuffer);
      // reset the buffer position
      m_audioBufferPos = 0;
      // tell the writer task to save the data
      xTaskNotify(m_processing_task_handle, 1, eIncrement);
    }
  }
};