#include "Transport.h"

Transport::Transport(OutputBuffer *output_buffer, size_t buffer_size)
{
  m_output_buffer = output_buffer;
  m_buffer_size = buffer_size;
  m_buffer = (uint8_t *)malloc(m_buffer_size);
  m_index = 0;
}

void Transport::add_sample(int16_t sample)
{
  m_buffer[m_index] = (sample + 32768) >> 8;
  m_index++;
  // have we reached a full packet?
  if (m_index == m_buffer_size)
  {
    if (m_should_send)
    {
      send();
    }
    m_index = 0;
  }
}