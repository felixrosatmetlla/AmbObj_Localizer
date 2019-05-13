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

#define noiseLevel 1e-5

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent   : public AudioAppComponent
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
    
    //==============================================================================
    std::complex<float>** getComplexFFTBuffer(float** fftBuffer, size_t numSamples);
    
    std::complex<float>* getPressure(std::complex<float>** complexFFTBuffer, size_t numSamples);
    std::complex<float>** getVelocityVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    
    float** getIntensityVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    float* getEnergyVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    
    void getDOA(std::complex<float>** complexFFTBuffer, size_t numSamples);
    void getDifuseness(std::complex<float>** complexFFTBuffer, size_t numSamples, int dt = 10);
    
    void cart2Sph(float** doa, size_t numSamples);

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
    float timestepAzi;
    float timestepEle;
    float timestepRad;
    
    // Interface Variables
    ToggleButton outputBttn;
    Slider energyThr;
    
    Random randNum;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
