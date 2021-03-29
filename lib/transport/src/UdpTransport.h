#pragma once

#include "Transport.h"

class OutputBuffer;
class AsyncUDP;

class UdpTransport : public Transport
{
private:
  AsyncUDP *udp;

protected:
  void send();

public:
  UdpTransport(OutputBuffer *output_buffer);
  bool begin() override;
};