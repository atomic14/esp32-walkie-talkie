#include <Arduino.h>

extern "C"
{
#include <SP_ENC.h>
#include <SP_DEC.h>
}

void encoder_task(void *param)
{
  while (true)
  {
    int16_t speech[160];
    for (int i = 0; i < 160; i++)
    {
      speech[i] = 0;
    }
    long long start = millis();
    for (int i = 0; i < 1000; i++)
    {
      int16_t encoded[20];
      speechEncoder(speech, encoded);
      //speechDecoder(encoded, speech);
    }
    long long end = millis();
    Serial.printf("Processing time = %lld\n", end - start);
    Serial.println("Finished");
    vTaskDelay(100);
  }
}

// esp_task_wdt_init(100, false);
// TaskHandle_t taskHandle;
// xTaskCreatePinnedToCore(encoder_task, "Encoder", 8192, NULL, 10, &taskHandle, 1);
