#pragma once

class I2SOutput;
class I2SSampler;
class UdpTransport;

class Application
{
private:
  I2SOutput *output;
  I2SSampler *input;
  UdpTransport *transport;

  void service();

public:
  Application();
  void begin();
  friend void samples_task(void *param);
};