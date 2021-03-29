#ifndef __sampler_base_h__
#define __sampler_base_h__

#include <Arduino.h>
#include <driver/i2s.h>

/**
 * Base Class for both the ADC and I2S sampler
 **/
class I2SSampler
{
private:
    // double buffer so we can be capturing samples while sending data
    int16_t *m_audioBuffer1;
    int16_t *m_audioBuffer2;
    // current position in the audio buffer
    int32_t m_audioBufferPos = 0;
    // current audio buffer
    int16_t *m_currentAudioBuffer;
    // buffer containing samples that have been captured already
    int16_t *m_capturedAudioBuffer;
    // size of the audio buffers in bytes
    int32_t m_bufferSizeInBytes;
    // size of the audio buffer in samples
    int32_t m_bufferSizeInSamples;
    // I2S reader task
    TaskHandle_t m_readerTaskHandle;
    // writer task
    TaskHandle_t m_writerTaskHandle;
    // i2s reader queue
    QueueHandle_t m_i2sQueue;
    // i2s port
    i2s_port_t m_i2sPort;

protected:
    void addSample(int16_t sample);
    virtual void configureI2S() = 0;
    virtual void processI2SData(uint8_t *i2sData, size_t bytesRead) = 0;
    i2s_port_t getI2SPort()
    {
        return m_i2sPort;
    }

public:
    int32_t getBufferSizeInBytes()
    {
        return m_bufferSizeInBytes;
    };
    int16_t *getCapturedAudioBuffer()
    {
        return m_capturedAudioBuffer;
    }
    void start(i2s_port_t i2sPort, i2s_config_t &i2sConfig, int32_t bufferSizeInSamples, TaskHandle_t writerTaskHandle);

    friend void i2sReaderTask(void *param);
};

#endif