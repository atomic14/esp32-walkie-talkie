#pragma once

class I2SOutput;
class I2SSampler;
class Transport;
class OutputBuffer;
class IndicatorLed;

typedef enum
{
  IDLE,         // the application is waiting for the user to push the transmit button or for audio packets to arrive
  TRANSMITTING, // the user has pushed the transmit button
  RECEIVING,    // we are receiving audio packets
} Application_State_t;

class Application
{
private:
  I2SOutput *m_output;
  I2SSampler *m_input;
  Transport *m_transport;
  IndicatorLed *m_indicator_led;
  OutputBuffer *m_output_buffer;

  Application_State_t m_current_state = IDLE;

  void service();

public:
  Application();
  void begin();
  void loop();
  friend void samples_task(void *param);
};