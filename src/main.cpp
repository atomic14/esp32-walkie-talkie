#include <Arduino.h>
#include <WiFi.h>
#include <driver/i2s.h>
#include "I2SOutput.h"
#include "I2SMEMSSampler.h"
#include "UdpTransport.h"

void setup()
{
  Serial.begin(115200);
  WiFi.mode(WIFI_STA);
  WiFi.begin("CMGResearch", "02087552867");
  if (WiFi.waitForConnectResult() != WL_CONNECTED)
  {
    Serial.println("Connection Failed! Rebooting...");
    delay(5000);
    ESP.restart();
  }
  Serial.println("Started");
}

void loop()
{
  delay(1000);
}
