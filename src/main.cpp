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
  application->loop();
  delay(50);
}
