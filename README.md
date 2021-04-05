# Overview

We've made a Walkie-Talkie using the ESP32.

[Explanatory video](https://www.youtube.com/watch?v=d_h38X4_eQQ)

[![Demo Video](https://img.youtube.com/vi/d_h38X4_eQQ/0.jpg)](https://www.youtube.com/watch?v=d_h38X4_eQQ)

Audio data is transmitted over either UDP broadcast or ESP-NOW. So the Walkie-Talkie will even work without a WiFi network!

I'm using my own microphone board (available on Tindie: https://www.tindie.com/products/21519/) but the code will work equally well with any I2S microphone (e.g. the INMP441) and you can easily modify it to use the built-in ADC for analogue microphones.

For output, I'm using an I2S amplifier breakout board which I'm using the drive a 4ohm speaker. Once again, you can modify the code to use the built-in DAC for output which will let you use headphones or an analogue amplifier board.

I've got a great series of videos on ESP32 Audio which are a great resource for anyone who wants to learn more about audio on the ESP32 which you can find here: https://www.youtube.com/playlist?list=PL5vDt5AALlRfGVUv2x7riDMIOX34udtKD

For this project I've 3D printed a case - you can access the Fusion 360 project here: https://a360.co/2PXgAUS

I've also created a custom PCB - you can access the schematic here: https://easyeda.com/chris_9044/esp32-walkie-talkie

The boards were manufactured by PCBWay and as always they've done a really great job. You can order the boards directly from PCBWay here: https://www.pcbway.com/project/shareproject/ESP32_Audio_Board_For_Walkie_Talkie.html

And you can help support the channel by using my referral link: https://www.pcbway.com/setinvite.aspx?inviteid=403566 for other PCBs.

However, you can also easily wire this up on breadboard - that's how I prototyped it. Everything is I2S based so it's just straightforward jumper wires.

# Setup

Everything is configured from the `src/config.h` file. To use UDP Broadcast comment out the line:

```
#define USE_ESP_NOW
```

Make sure you update the WiFi SSID and Password:

```
// WiFi credentials
#define WIFI_SSID << YOUR_SSID >>
#define WIFI_PSWD << YOUR_PASSWORD >>
```

The pins for the microphone and the amplifier board are all setup in the same `config.h` file.

# Building and Running

I'm using PlatformIO for this project so you will need to have that installed. Open up the project and connect your ESP32. You should be able to just hit build and run.

Obviously, you'll need two ESP32 boards and components to do anything :)
