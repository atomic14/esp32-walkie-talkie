#pragma once

class I2SOutput;
class I2SSampler;
class Transport;
class OutputBuffer;
class IndicatorLed;

class Application
{
private:
  I2SOutput *m_output;
  I2SSampler *m_input;
  Transport *m_transport;
  IndicatorLed *m_indicator_led;
  OutputBuffer *m_output_buffer;

public:
  Application();
  void begin();
  void loop();
};