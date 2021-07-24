
#include "Output.h"
#include <esp_log.h>
#include <driver/i2s.h>

static const char *TAG = "OUT";

// number of frames to try and send at once (a frame is a left and right sample)
const int NUM_FRAMES_TO_SEND = 256;

Output::Output(i2s_port_t i2s_port) : m_i2s_port(i2s_port)
{
  // this will contain the prepared samples for sending to the I2S device
  m_frames = (int16_t *)malloc(2 * sizeof(int16_t) * NUM_FRAMES_TO_SEND);
}

Output::~Output()
{
  free(m_frames);
}

void Output::stop()
{
  // stop the i2S driver
  i2s_stop(m_i2s_port);
  i2s_driver_uninstall(m_i2s_port);
}

void Output::write(int16_t *samples, int count)
{
  int sample_index = 0;
  while (sample_index < count)
  {
    int samples_to_send = 0;
    for (int i = 0; i < NUM_FRAMES_TO_SEND && sample_index < count; i++)
    {
      int sample = process_sample(samples[sample_index]);
      m_frames[i * 2] = sample;  // left channel
      m_frames[i * 2 + 1] = sample;  // right channel
      samples_to_send++;
      sample_index++;
    }
    // write data to the i2s peripheral
    size_t bytes_written = 0;
    i2s_write(m_i2s_port, m_frames, samples_to_send * sizeof(int16_t) * 2, &bytes_written, portMAX_DELAY);
    if (bytes_written != samples_to_send * sizeof(int16_t) * 2)
    {
      ESP_LOGE(TAG, "Did not write all bytes");
    }
  }
}
