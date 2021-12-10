#include <Arduino.h>
#include <freertos/FreeRTOS.h>
#include "IndicatorLed.h"

void update_indicator_task(void *param)
{
  IndicatorLed *indicator = reinterpret_cast<IndicatorLed *>(param);
  while (true)
  {
    if (indicator->m_is_flashing)
    {
      indicator->set_led_rgb(indicator->m_flash_color);
      vTaskDelay(100);
    }
    indicator->set_led_rgb(indicator->m_default_color);
    vTaskDelay(100);
  }
}

void IndicatorLed::begin()
{
  TaskHandle_t task_handle;
  xTaskCreate(update_indicator_task, "Indicator LED Task", 4096, this, 0, &task_handle);
}

void IndicatorLed::set_is_flashing(bool is_flashing, uint32_t flash_color)
{
  m_is_flashing = is_flashing;
  m_flash_color = flash_color;
}
void IndicatorLed::set_default_color(uint32_t color)
{
  m_default_color = color;
}
