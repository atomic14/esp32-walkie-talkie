#include "Arduino.h"
#include "GenericDevBoardIndicatorLed.h"

#ifdef LED_BUILTIN
const uint8_t BUILT_IN_LED = LED_BUILTIN;
#else
const uint8_t BUILT_IN_LED = GPIO_NUM_2;
#endif

GenericDevBoardIndicatorLed::GenericDevBoardIndicatorLed()
{
  pinMode(BUILT_IN_LED, OUTPUT);
}

// we don't really have any colors so just use the built in LED
void GenericDevBoardIndicatorLed::set_led_rgb(uint32_t color)
{
  if (color == 0)
  {
    digitalWrite(BUILT_IN_LED, LOW);
  }
  else
  {
    digitalWrite(BUILT_IN_LED, HIGH);
  }
}
