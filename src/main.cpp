#include <Arduino.h>
#include "Application.h"

// our application
Application *application;

void setup()
{
  Serial.begin(115200);
  // start up the application
  application = new Application();
  application->begin();
  Serial.println("Application started");
}

void loop()
{
  // nothing to do - the application is doing all the work
  vTaskDelay(pdMS_TO_TICKS(1000));
}
