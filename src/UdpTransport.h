#pragma once

class I2SOutput;
class AsyncUDP;

class UdpTransport
{
private:
  AsyncUDP *udp;

public:
  UdpTransport();
  bool begin(I2SOutput *audio_output);
  void send_audio(uint8_t *data, int bytes);
};