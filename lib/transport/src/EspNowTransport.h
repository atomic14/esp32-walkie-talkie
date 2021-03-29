#pragma once

#include "Transport.h"

class OutputBuffer;

class EspNowTransport: public Transport {
protected:
  void send();
public:
  EspNowTransport(OutputBuffer *output_buffer);
  virtual bool begin() override;
  friend void receiveCallback(const uint8_t *macAddr, const uint8_t *data, int dataLen);
};