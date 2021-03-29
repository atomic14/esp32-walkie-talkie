#pragma once

#include "IndicatorLed.h"

class TinyPICO;
class TinyPICOIndicatorLed : public IndicatorLed
{
private:
  TinyPICO *m_tp = NULL;

protected:
  void set_led_rgb(uint32_t color);

public:
  TinyPICOIndicatorLed();
};