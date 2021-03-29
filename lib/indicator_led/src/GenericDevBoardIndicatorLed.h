#pragma once

#include "IndicatorLed.h"

class GenericDevBoardIndicatorLed : public IndicatorLed
{
protected:
  void set_led_rgb(uint32_t color);

public:
  GenericDevBoardIndicatorLed();
};