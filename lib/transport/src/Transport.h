#pragma once
#include <stdlib.h>
#include <stdint.h>

class OutputBuffer;

class Transport
{
protected:
  // audio buffer for samples we need to send
  uint8_t *m_buffer = NULL;
  int m_buffer_size = 0;
  int m_index = 0;
  // should we send
  bool m_should_send = false;
  // when was the last time we recieved a packed?
  unsigned long m_last_packet_received = 0;

  OutputBuffer *m_output_buffer = NULL;

  virtual void send() = 0;

public:
  Transport(OutputBuffer *output_buffer, size_t buffer_size);
  void add_sample(int16_t sample);
  void should_send(bool should_send)
  {
    m_should_send = should_send;
  }
  unsigned long get_last_packet_received()
  {
    return m_last_packet_received;
  }
  virtual bool begin() = 0;
};