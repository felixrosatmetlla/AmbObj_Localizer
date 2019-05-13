/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent():forwardFFT(9)
{
    
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (800, 600);
    
    outputBttn.setBounds(0, 0, 100, 100);
    outputBttn.setButtonText("Output Audio");
    outputBttn.changeWidthToFitText();

    addAndMakeVisible(outputBttn);
    
    energyThr.setBounds(0, 150, 100, 120);
    energyThr.setRange(0, 2);
    energyThr.setSliderStyle(Slider::SliderStyle::LinearVertical);
    energyThr.setValue(0.5);
    energyThr.setTextBoxStyle(Slider::TextBoxBelow, true, 70, 20);
    addAndMakeVisible(energyThr);
    
    // Some platforms require permissions to open input channels so request that here
    if (RuntimePermissions::isRequired (RuntimePermissions::recordAudio)
        && ! RuntimePermissions::isGranted (RuntimePermissions::recordAudio))
    {
        RuntimePermissions::request (RuntimePermissions::recordAudio,
                                     [&] (bool granted) { if (granted)  setAudioChannels (2, 2); });
    }
    else
    {
        // Specify the number of input and output channels that we want to open
        MainComponent::setAudioChannels (N_CH, N_CH);
    }
}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    // This function will be called when the audio device is started, or when
    // its settings (i.e. sample rate, block size, etc) are changed.

    // You can use this function to initialise any resources you might need,
    // but be careful - it will be called on the audio thread, not the GUI thread.

    // For more details, see the help for AudioProcessor::prepareToPlay()
    
    size_t num_frames = samplesPerBlockExpected;
    size_t sample_size_bytes = sizeof(float);
    //size_t num_channels = 4;
    //size_t num_samples = num_frames * num_channels;
    
    audioBuffer = (float**)calloc(N_CH, sizeof(float*));
    for(int i = 0; i < N_CH; i++ ){
        audioBuffer[i] = (float*)calloc(num_frames,sample_size_bytes);
    }
    
    fftBuffer = (float**)calloc(N_CH, sizeof(float*));
    for(int i = 0; i < N_CH; i++ ){
        fftBuffer[i] = (float*)calloc(num_frames*2,sample_size_bytes);
    }
    
    
    bufferEnergy = (float*)calloc(N_CH,sample_size_bytes);
    
    doa = (float**)calloc(N_CH-1, sizeof(float*));
    for(int i = 0; i < N_CH-1; i++ ){
        doa[i] = (float*)calloc(num_frames*2,sample_size_bytes);
    }
    
    diffuseness = (float*)calloc(num_frames,sample_size_bytes);
    
    azimuth = (float*)calloc(num_frames,sample_size_bytes);
    elevation = (float*)calloc(num_frames,sample_size_bytes);
    radius = (float*)calloc(num_frames,sample_size_bytes);
    
    bufferCounter = 0;

}

void MainComponent::getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill)
{
    // Your audio-processing code goes here!

    // For more details, see the help for AudioProcessor::getNextAudioBlock()

    // Right now we are not producing any data, in which case we need to clear the buffer
    // (to prevent the output of random noise)
    bufferEnergy[0] = 0;
    bufferEnergy[1] = 0;
    bufferEnergy[2] = 0;
    bufferEnergy[3] = 0;
    
    bufferAzi = 0;
    bufferEle = 0;
    bufferRad = 0;
    
    // Get input audio from microphone
    for(int ch = 0; ch < bufferToFill.buffer->getNumChannels(); ch++ ){
        for(int sample = 0; sample < bufferToFill.buffer->getNumSamples(); sample++ ){
            audioBuffer[ch][sample] = *bufferToFill.buffer->getWritePointer(ch, sample);
            //std::cout << "ch_" << ch << ": " << audioBuffer[ch][sample] << std::endl;
            audioBuffer[ch][sample] += (randNum.nextFloat() - 0.5) * noiseLevel;
            bufferEnergy[ch] = bufferEnergy[ch] + pow(*bufferToFill.buffer->getReadPointer(ch, sample),2);
        
        }
        bufferEnergy[ch] = bufferEnergy[ch]/bufferToFill.buffer->getNumSamples();
        
        std::memcpy(fftBuffer[ch], audioBuffer[ch], sizeof(float)*bufferToFill.buffer->getNumSamples());
        
    }

    if(bufferEnergy[0] > energyThr.getValue() || bufferEnergy[1] > energyThr.getValue()){
        for(int ch=0; ch<4; ch++){
            forwardFFT.performRealOnlyForwardTransform(fftBuffer[ch], true);
        }
        
        /*for(int ch=0; ch<4; ch++){
            for(int sample = 0; sample<bufferToFill.buffer->getNumSamples(); sample++){
                std::cout << "fft: " << fftBuffer[ch][sample] << std::endl;
                
            }
        }*/
        
        std::complex<float>** complexFFTBuffer = getComplexFFTBuffer(fftBuffer, bufferToFill.buffer->getNumSamples());
        
        /*for(int ch=0; ch<4; ch++){
            for(int sample = 0; sample<bufferToFill.buffer->getNumSamples(); sample++){
                std::cout << "Complex fft: " << complexFFTBuffer[ch][sample] << std::endl;
                
            }
        }*/

        
        getDOA(complexFFTBuffer, bufferToFill.buffer->getNumSamples());
        
        /*for(int sample = 0; sample<bufferToFill.buffer->getNumSamples(); sample++){
            std::cout << "azi: " << azimuth[sample] << std::endl;
            std::cout << "ele: " << elevation[sample] << std::endl;
            //std::cout << "radius: " << radius[sample] << std::endl;
        }*/
        
        getDifuseness(complexFFTBuffer, bufferToFill.buffer->getNumSamples());
        
        /*for(int sample = 0; sample<bufferToFill.buffer->getNumSamples(); sample++){
            std::cout << "difuseness: " << diffuseness[sample] << std::endl;

        }*/
        
        int counter = 0;
        
        for(int sample = 0; sample < bufferToFill.buffer->getNumSamples(); sample++){
            if(diffuseness[sample] > 0.9){
                if(!std::isnan(azimuth[sample])){
                    
                    bufferAzi += azimuth[sample];
                    bufferEle += elevation[sample];
                    bufferRad += radius[sample];
                    
                    counter++;
                }
            }
        }
        
        bufferAzi = bufferAzi/counter;
        bufferEle = bufferEle/counter;
        bufferRad = bufferRad/counter;
        
        
        if(bufferCounter == 100){
            timestepAzi = timestepAzi/bufferCounter;
            timestepEle = timestepEle/bufferCounter;
            timestepRad = timestepRad/bufferCounter;
            
            std::cout << "azimuth: " << timestepAzi << std::endl;
            std::cout << "elevation: " << timestepEle << std::endl;
            std::cout << "radius: " << timestepRad << std::endl;
            
            bufferCounter = 0;
            timestepAzi = 0;
            timestepEle = 0;
            timestepRad = 0;
        }
        else{
            timestepAzi += bufferAzi;
            timestepEle += bufferEle;
            timestepRad += bufferRad;
            
        }
        
        
        bufferCounter++;
        //std::cout << "azi: " << bufferAzi << std::endl;
        //std::cout << "ele: " << bufferEle << std::endl;
        //std::cout << "rad: " << bufferRad << std::endl;
        /*// FFT working check
        for(int ch=0; ch<2; ch++){
             forwardFFT.performRealOnlyInverseTransform(fftBuffer[ch]);
             }
         
             for(int ch = 0; ch < bufferToFill.buffer->getNumChannels(); ch++ ){
             for(int sample = 0; sample < bufferToFill.buffer->getNumSamples(); sample++ ){
                 *bufferToFill.buffer->getWritePointer(ch, sample) = fftBuffer[ch][sample];
             }
         }*/
    }
    else{
        bufferToFill.clearActiveBufferRegion();
    }
    
    
    if(outputBttn.getToggleState() == false){
        bufferToFill.clearActiveBufferRegion();
    }
    
}

void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
    
    for(int ch = 0; ch < N_CH; ch++){
        free(audioBuffer[ch]);
    }
    free(audioBuffer);
    
    for(int ch = 0; ch < N_CH; ch++){
        free(fftBuffer[ch]);
    }
    free(fftBuffer);
    
    free(bufferEnergy);
    
    for(int ch = 0; ch < N_CH-1; ch++){
        free(doa[ch]);
    }
    free(doa);
    
    free(diffuseness);
    
    free(azimuth);
    free(elevation);
    free(radius);


}

//==============================================================================
std::complex<float>** MainComponent::getComplexFFTBuffer(float** fftBuffer, size_t numSamples)
{
    std::complex<float>** complexFFTBuffer = (std::complex<float>**)calloc(numSamples,sizeof(std::complex<float>*));
    for(int i = 0; i < N_CH; i++ ){
        complexFFTBuffer[i] = (std::complex<float>*)calloc(numSamples,sizeof(std::complex<float>));
    }
    
    for(int channel = 0; channel < N_CH; channel++){
        for(int sample = 0, fftsample = 0; sample < numSamples/2; sample++, fftsample+=2){
            complexFFTBuffer[channel][sample].real(fftBuffer[channel][fftsample]);
            complexFFTBuffer[channel][sample].imag(fftBuffer[channel][fftsample+1]);
            //std::cout << "complex: " << complexFFTBuffer[channel][sample] << std::endl;
        }
    }
    
    return complexFFTBuffer;
}

std::complex<float>* MainComponent::getPressure(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    std::complex<float>* p_k = (std::complex<float>*)calloc(numSamples,sizeof(std::complex<float>));
    std::memcpy(p_k, complexFFTBuffer[0], sizeof(std::complex<float>)*numSamples);
    
    /*for(int sample = 0; sample<numSamples; sample++){
        std::cout << "p_k: " << p_k[sample] << std::endl;
    }*/
    
    return p_k;
}

std::complex<float>** MainComponent::getVelocityVector(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    float scale = -1.0/(MathConstants<float>::sqrt2 * p0 * SOUND_SPEED);
    
    std::complex<float>** u_k = (std::complex<float>**)calloc(N_CH-1, sizeof(std::complex<float>*));
    for(int i = 0; i < N_CH-1; i++ ){
        u_k[i] = (std::complex<float>*)calloc(numSamples,sizeof(std::complex<float>));
        
    }
    
    for(int i = 1; i < N_CH; i++ ){
        std::memcpy(u_k[i-1], complexFFTBuffer[i], sizeof(std::complex<float>)*numSamples);
        
        for(int j = 0; j < numSamples; j++){
            u_k[i-1][j] =  scale * u_k[i-1][j];
            //std::cout << "u_k: " << u_k[i-1][j] << std::endl;
        }
    }
    
    return u_k;
}

float** MainComponent::getIntensityVector(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    std::complex<float>* p_k = getPressure(complexFFTBuffer, numSamples);
    std::complex<float>** u_k = getVelocityVector(complexFFTBuffer, numSamples);
    
    float** i_k = (float**)calloc(N_CH-1, sizeof(float*));
    for(int i = 0; i < N_CH-1; i++ ){
        i_k[i] = (float*)calloc(numSamples,sizeof(float));
        
        for(int j = 0; j < numSamples; j++){
            std::complex<float> pk_uk = p_k[j] * std::conj(u_k[i][j]);
            i_k[i][j] = 0.5 * pk_uk.real();
            
            //std::cout << "i_k: " << i_k[i][j] << std::endl;

        }
    }
    
    return i_k;
}

float* MainComponent::getEnergyVector(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    std::complex<float>* p_k = getPressure(complexFFTBuffer, numSamples);
    std::complex<float>** u_k = getVelocityVector(complexFFTBuffer, numSamples);
    
    float* s1 = (float*)calloc(numSamples,sizeof(float));
    float* s2 = (float*)calloc(numSamples,sizeof(float));
    
    float* energy = (float*)calloc(numSamples,sizeof(float));

    for(int sample = 0; sample < numSamples; sample++){
        float u_norm = std::sqrt(u_k[0][sample] * u_k[0][sample] + u_k[1][sample] * u_k[1][sample] + u_k[2][sample] * u_k[2][sample]).real();
        s1[sample] = std::pow(u_norm, 2);
        
        s2[sample] = std::pow(std::abs(p_k[sample]), 2);
        
        energy[sample] = ((p0/4.0) * s1[sample]) + ((1.0/(4 * p0 * std::pow(SOUND_SPEED, 2))) * s2[sample]);
    }
    
    return energy;
    
}

void MainComponent::getDOA(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    float** i_k = getIntensityVector(complexFFTBuffer, numSamples);
    
    
    float* i_norm = (float*)calloc(numSamples,sizeof(float));
    
    for(int sample = 0; sample < numSamples; sample++){
        i_norm[sample] = std::sqrt(std::pow(i_k[0][sample], 2) + std::pow(i_k[1][sample], 2) + std::pow(i_k[2][sample], 2));
        //std::cout << "i_norm: " << i_norm[sample] << std::endl;
    }
    
    for(int channel = 0; channel < N_CH-1; channel++){
        for(int sample = 0; sample < numSamples; sample++){
            doa[channel][sample] = -i_k[channel][sample]/i_norm[sample];
            //std::cout << "doa: " << doa[channel][sample] << std::endl;

        }
    }
    
    cart2Sph(doa, numSamples);

}

void MainComponent::cart2Sph(float** doa, size_t numSamples){
    
    for(int sample = 0; sample < numSamples; sample++){
        float x_y = (doa[0][sample] * doa[0][sample]) + (doa[1][sample] * doa[1][sample]);
        
        radius[sample] = std::sqrt((doa[0][sample] * doa[0][sample]) + (doa[1][sample] * doa[1][sample]) + (doa[2][sample] * doa[2][sample]));
        
        azimuth[sample] = std::atan2(doa[1][sample], doa[0][sample]);
        
        elevation[sample] = std::atan2(doa[2][sample], std::sqrt(x_y));
    }
    
}

void MainComponent::getDifuseness(std::complex<float>** complexFFTBuffer, size_t numSamples, int dt)
{
    float** i_k = getIntensityVector(complexFFTBuffer, numSamples);
    float* e_k = getEnergyVector(complexFFTBuffer, numSamples);
    
    float* i_mean = (float*)calloc(N_CH-1,sizeof(float));
    float e_mean = 0;
    
    for(int i = dt/2; i < numSamples - dt/2; i++){
        for(int ch = 0; ch < N_CH-1; ch++){
            for(int sample = 0; sample < dt; sample++){
                i_mean[ch] += i_k[ch][i + sample];
                e_mean += e_k[i + sample];
            }
            i_mean[ch] = i_mean[ch]/dt;
            
            float num = std::sqrt(std::pow(i_mean[0], 2) + std::pow(i_mean[1], 2) + std::pow(i_mean[2], 2));
            float den = SOUND_SPEED * e_mean;
            
            diffuseness[i] = 1 - (num/den);
        }
    }
    
    for(int i = 0; i < dt/2; i++){
        diffuseness[i] = diffuseness[dt/2];
    }
    
    for(int i = (int)numSamples - dt/2; i < numSamples; i++){
        diffuseness[i] = diffuseness[numSamples - dt/2 -1];
    }
}

//==============================================================================
void MainComponent::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    // You can add your drawing code here!
}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
}