/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"

#define SOUND_SPEED 346.13
#define p0 1.1839
#define N_CH 4

#define visualizer_W 285
#define visualizer_H 90

#define azi_resolution 314
#define ele_resolution 157


//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent   : public AudioAppComponent, private Timer
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent();

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (Graphics& g) override;
    void resized() override;
    
    void paintDoA ();
    
    //==============================================================================
    std::complex<float>** getComplexFFTBuffer(float** fftBuffer, size_t fftSize);
    
    std::complex<float>* getPressure(std::complex<float>** complexFFTBuffer, size_t numSamples);
    std::complex<float>** getVelocityVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    
    float** getIntensityVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    float* getEnergyVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    
    void getDOA(std::complex<float>** complexFFTBuffer, size_t numSamples);
    void getDifuseness(std::complex<float>** complexFFTBuffer, size_t numSamples, int dt = 10);
    
    void cart2Sph(float** doa, size_t numSamples);

    //==============================================================================
    void timerCallback() override;
    
    
private:
    //==============================================================================
    // Audio Variables
    float** audioBuffer;
    float* bufferEnergy;
    
    dsp::FFT forwardFFT;
    float** fftBuffer;
    
    float** doa;
    float* diffuseness;
    
    float* azimuth;
    float* elevation;
    float* radius;
    
    float bufferAzi;
    float bufferEle;
    float bufferRad;
    
    int bufferCounter;
    float timeStepAziSin;
    float timeStepAziCos;
    float timestepAzi;
    float timestepEle;
    float timestepRad;
    
    float resultAzi;
    float resultEle;
    float resultRad;
    
    // Interface Variables
    ToggleButton outputBttn;
    Slider energyThr;
    Slider numBuffers;
    Slider noiseThr;
    Slider diffusenesThr;
    
    Label energyLabel;
    Label numBuffersLabel;
    Label noiseLabel;
    Label diffusenesLabel;

    Random randNum;
    
    bool nextBlockReady = false;
    float point_x;
    float point_y;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
