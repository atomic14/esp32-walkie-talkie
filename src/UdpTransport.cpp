#include <Arduino.h>
#include <AsyncUDP.h>
#include "UdpTransport.h"
#include "I2SOutput.h"

const int MAX_UDP_SIZE = 1436;

UdpTransport::UdpTransport()
{
}

bool UdpTransport::begin(I2SOutput *output)
{
  udp = new AsyncUDP();
  if (udp->listenMulticast(IPAddress(239, 1, 2, 3), 8192, 25, TCPIP_ADAPTER_IF_STA))
  {
    udp->onPacket([output](AsyncUDPPacket packet) {
      String src = packet.remoteIP().toString();
      int count = packet.length() / sizeof(int16_t);
      int16_t *samples = reinterpret_cast<int16_t *>(packet.data());
      output->push_samples(samples, count);
    });
    return true;
  }
  Serial.println("Failed to listen");
  return false;
}

void UdpTransport::send_audio(uint8_t *data, int length)
{
  if (length > MAX_UDP_SIZE)
  {
    Serial.println("Dropping bytes from message");
  }
  udp->write(data, std::min(MAX_UDP_SIZE, length));
}