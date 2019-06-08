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
    void createSlider(Slider &slider, int x, int y, int width, int height, double minRange, double maxRange, double value, Slider::SliderStyle style, Slider::TextEntryBoxPosition textboxPos, bool textboxRead, int textboxWidth, int textboxHeight, bool isSkewed, double interval = 0, double skewFactor = 0);
    
    void createLabel(Label &label ,const String &labelText, NotificationType notification, Component *attachComp, bool onLeft = false);
    
    void createToggleButton(ToggleButton &toggleButton, int x, int y, int width, int height, const String &buttonText);
    
    //==============================================================================
    void paint (Graphics& g) override;
    void resized() override;
    
    void paintDoA ();
    
    //==============================================================================
    void getComplexFFTBuffer(float** fftBuffer, size_t fftSize);
    
    void getPressure(std::complex<float>** complexFFTBuffer, size_t numSamples);
    void getVelocityVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    
    void getIntensityVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    void getEnergyVector(std::complex<float>** complexFFTBuffer, size_t numSamples);
    
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
    
    std::complex<float>** complexFFTBuffer;
    std::complex<float>* p_k;
    std::complex<float>** u_k;
    float** i_k;
    
    float* s1;
    float* s2;
    float* e_k;
    
    float* i_mean;
    float* i_norm;
    
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
