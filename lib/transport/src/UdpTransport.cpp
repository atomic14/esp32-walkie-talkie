#include <Arduino.h>
#include <AsyncUDP.h>
#include "UdpTransport.h"
#include "OutputBuffer.h"

const int MAX_UDP_SIZE = 1436;

UdpTransport::UdpTransport(OutputBuffer *output_buffer) : Transport(output_buffer, MAX_UDP_SIZE)
{
}

unsigned long last_packet;
bool UdpTransport::begin()
{
  udp = new AsyncUDP();
  last_packet = millis();
  if (udp->listen(8192))
  {
    udp->onPacket([this](AsyncUDPPacket packet)
                  {
                    // our packets contain unsigned 8 bit PCM samples
                    // so we can push them straight into the output buffer
                    if ((packet.length() > this->m_header_size) && (packet.length() <= MAX_UDP_SIZE) && (memcmp(packet.data(), this->m_buffer, this->m_header_size) == 0)) 
                    {
                      this->m_output_buffer->add_samples(packet.data() + m_header_size, packet.length() - m_header_size);
                    }
                  });
    return true;
  }
  Serial.println("Failed to listen");
  return false;
}

void UdpTransport::send()
{
  udp->broadcast(m_buffer, m_index);
}
