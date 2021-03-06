/*
  ==============================================================================

    This file was auto-generated!

  ==============================================================================
*/

#include "MainComponent.h"

//==============================================================================
MainComponent::MainComponent():forwardFFT(9), data("Data")
{
    
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (960, 530);
    startTimerHz(60);
    
    createToggleButton(outputBttn, 0, 0, 100, 100, "Output Audio");
    
    createLabel(energyLabel, "Energy Thr.", dontSendNotification, &energyThr);
    createSlider(energyThr, 25, 100, 100, 120, 0, 2, 0.5, Slider::SliderStyle::RotaryVerticalDrag, Slider::TextBoxBelow, true, 70, 20, false);

    createLabel(numBuffersLabel, "Number Buffers", dontSendNotification, &numBuffers);
    createSlider(numBuffers, 150, 100, 100, 120, 1, 100, 10, Slider::SliderStyle::RotaryVerticalDrag, Slider::TextBoxBelow, true, 70, 20, false, 1);

    createLabel(noiseLabel, "Noise Level", dontSendNotification, &noiseThr);
    createSlider(noiseThr, 150, 300, 100, 120, 0, 1e-1, 1e-6, Slider::SliderStyle::RotaryVerticalDrag, Slider::TextBoxBelow, true, 70, 20, true, 0, 1e-4);
    
    
    createLabel(diffusenesLabel, "Diffuseness Thr.", dontSendNotification, &diffusenesThr);
    createSlider(diffusenesThr, 25, 300, 100, 120, 0, 1, 0.9, Slider::SliderStyle::RotaryVerticalDrag, Slider::TextBoxBelow, true, 70, 20, false);
    
    numBuffers.onValueChange = [this] { bufferCounter = 0; };
    
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
    
    String dataToString = data.createDocument(String());
    data.writeToFile(File("/Users/felixrosatmetlla/Desktop/TFG/AmbObj_Localizer/AmbObjLoc/AppData/appdata.xml"), String());
}

void MainComponent::createToggleButton(ToggleButton &toggleButton, int x, int y, int width, int height, const String &buttonText){
    
    toggleButton.setBounds(x, y, width, height);
    toggleButton.setButtonText("Output Audio");
    toggleButton.changeWidthToFitText();
    
    addAndMakeVisible(toggleButton);
}

void MainComponent::createLabel(Label &label ,const String &labelText, NotificationType notification, Component *attachComp, bool onLeft){
    
    label.setText(labelText, notification);
    label.attachToComponent(attachComp, onLeft);
    
    addAndMakeVisible(label);
}

void MainComponent::createSlider(Slider &slider, int x, int y, int width, int height, double minRange, double maxRange, double value, Slider::SliderStyle style, Slider::TextEntryBoxPosition textboxPos, bool textboxRead, int textboxWidth, int textboxHeight, bool isSkewed, double interval, double skewFactor){
    
    slider.setBounds(x, y, width, height);
    slider.setRange(minRange, maxRange, interval);
    slider.setValue(value);
    slider.setSliderStyle(style);
    slider.setTextBoxStyle(textboxPos, textboxRead, textboxWidth, textboxHeight);
    
    if(isSkewed){
        slider.setSkewFactorFromMidPoint(skewFactor);
    }
    addAndMakeVisible(slider);
    
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
        fftBuffer[i] = (float*)calloc(forwardFFT.getSize()*2,sample_size_bytes);
    }
    
    bufferEnergy = (float*)calloc(N_CH,sample_size_bytes);
    
    complexFFTBuffer = (std::complex<float>**)calloc(N_CH,sizeof(std::complex<float>*));
    for(int i = 0; i < N_CH; i++ ){
        complexFFTBuffer[i] = (std::complex<float>*)calloc(forwardFFT.getSize(),sizeof(std::complex<float>));
    }
    
    p_k = (std::complex<float>*)calloc(forwardFFT.getSize() + 1,sizeof(std::complex<float>));
    
    u_k = (std::complex<float>**)calloc(N_CH-1, sizeof(std::complex<float>*));
    for(int i = 0; i < N_CH-1; i++ ){
        u_k[i] = (std::complex<float>*)calloc(forwardFFT.getSize() + 1,sizeof(std::complex<float>));
        
    }
    
    i_k = (float**)calloc(N_CH-1, sizeof(float*));
    for(int i = 0; i < N_CH-1; i++ ){
        i_k[i] = (float*)calloc(forwardFFT.getSize() + 1,sizeof(float));
    }
    
    s1 = (float*)calloc(forwardFFT.getSize() + 1,sizeof(float));
    s2 = (float*)calloc(forwardFFT.getSize() + 1,sizeof(float));

    e_k = (float*)calloc(forwardFFT.getSize() + 1,sizeof(float));

    i_norm = (float*)calloc(forwardFFT.getSize() + 1,sizeof(float));

    doa = (float**)calloc(N_CH-1, sizeof(float*));
    for(int i = 0; i < N_CH-1; i++ ){
        doa[i] = (float*)calloc(num_frames*2,sample_size_bytes);
    }
    
    i_mean = (float*)calloc(N_CH-1,sizeof(float));
    diffuseness = (float*)calloc(num_frames,sample_size_bytes);
    
    azimuth = (float*)calloc(num_frames,sample_size_bytes);
    elevation = (float*)calloc(num_frames,sample_size_bytes);
    radius = (float*)calloc(num_frames,sample_size_bytes);
    
    bufferCounter = 0;
    
    timeStepAziSin = 0;
    timeStepAziCos = 0;
    timestepEle = 0;
    timestepAzi = 0;
    timestepRad = 0;

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
    
    float aziSin = 0;
    float aziCos = 0;
    
    
    
    // Get input audio from microphone
    for(int ch = 0; ch < bufferToFill.buffer->getNumChannels(); ch++ ){
        for(int sample = 0; sample < bufferToFill.buffer->getNumSamples(); sample++ ){
            audioBuffer[ch][sample] = *bufferToFill.buffer->getWritePointer(ch, sample);
            //std::cout << "ch_" << ch << ": " << audioBuffer[ch][sample] << std::endl;
            audioBuffer[ch][sample] += (randNum.nextFloat() - 0.5) * noiseThr.getValue();
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
            for(int sample = 0; sample<bufferToFill.buffer->getNumSamples()*2; sample++){
                if(sample == bufferToFill.buffer->getNumSamples()/2+1){
                    std::cout << "size: " << bufferToFill.buffer->getNumSamples()/2+1 << std::endl;
                    std::cout << "fft: " << fftBuffer[ch][sample] << std::endl;
                }
                if(sample == bufferToFill.buffer->getNumSamples()+1){
                    std::cout << "size: " << bufferToFill.buffer->getNumSamples() << std::endl;
                    std::cout << "fft: " << fftBuffer[ch][sample] << std::endl;
                }
                if(sample == bufferToFill.buffer->getNumSamples()*2-7){
                    std::cout << "size: " << bufferToFill.buffer->getNumSamples()*2 << std::endl;
                    std::cout << "fft: " << fftBuffer[ch][sample] << std::endl;
                }
                //std::cout << "fft: " << fftBuffer[ch][sample] << std::endl;
                
            }
        }*/
        
        getComplexFFTBuffer(fftBuffer, forwardFFT.getSize() + 1);
        
        /*for(int ch=0; ch<4; ch++){
            for(int sample = 0; sample<bufferToFill.buffer->getNumSamples(); sample++){
                std::cout << "Complex fft: " << complexFFTBuffer[ch][sample] << std::endl;
                
            }
        }*/

        
        getDOA(complexFFTBuffer, forwardFFT.getSize() + 1);
        
        /*for(int sample = 0; sample<forwardFFT.getSize() + 1; sample++){
            //std::cout << "azi: " << azimuth[sample] << std::endl;
            //std::cout << "ele: " << elevation[sample] << std::endl;
            //std::cout << "radius: " << radius[sample] << std::endl;
        }*/
        
        getDifuseness(complexFFTBuffer, forwardFFT.getSize() + 1);
        
        /*for(int sample = 0; sample<bufferToFill.buffer->getNumSamples(); sample++){
            std::cout << "difuseness: " << diffuseness[sample] << std::endl;

        }*/
        
        int counter = 0;
        
        for(int sample = 0; sample < forwardFFT.getSize() + 1; sample++){
            if(diffuseness[sample] > diffusenesThr.getValue()){
                if(!std::isnan(azimuth[sample])){
                    
                    aziSin += std::sin(azimuth[sample]);
                    aziCos += std::cos(azimuth[sample]);
                    
                    //bufferAzi += azimuth[sample];
                    bufferEle += elevation[sample];
                    bufferRad += radius[sample];
                    
                    counter++;
                }
            }
        }
        
        //aziSin = aziSin/counter;
        //aziCos = aziCos/counter;
        
        bufferAzi = std::atan2(aziSin, aziCos);
        if(bufferAzi < 0){
            bufferAzi += MathConstants<float>::twoPi;
        }
        //std::cout << "azimuth: " << bufferAzi << std::endl;

        //bufferAzi = bufferAzi/counter;
        bufferEle = bufferEle/counter;
        bufferRad = bufferRad/counter;
        
        
        if(bufferCounter == numBuffers.getValue()){
            
            timestepAzi = std::atan2(timeStepAziSin, timeStepAziCos);
            if(timestepAzi < 0){
                timestepAzi += MathConstants<float>::twoPi;
            }
            
            timestepEle = timestepEle/bufferCounter;
            timestepRad = timestepRad/bufferCounter;
            
            resultAzi = timestepAzi;
            resultEle = timestepEle;
            resultRad = timestepRad;
            
            nextBlockReady = true;
            MainComponent::writeXMLData (resultAzi, resultEle);
            std::cout << "azimuth: " << resultAzi << std::endl;
            std::cout << "elevation: " << resultEle << std::endl;
            //std::cout << "radius: " << resultRad << std::endl;
            
            bufferCounter = 0;
            
            timeStepAziSin = 0;
            timeStepAziCos = 0;
            
            timestepAzi = 0;
            timestepEle = 0;
            timestepRad = 0;
            
            
        }
        else{
            
            timeStepAziSin += std::sin(bufferAzi);
            timeStepAziCos += std::cos(bufferAzi);
            
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
        free(complexFFTBuffer[ch]);
    }
    free(complexFFTBuffer);
    
    free(p_k);
    
    for(int ch = 0; ch < N_CH-1; ch++){
        free(u_k[ch]);
    }
    free(u_k);
    
    for(int ch = 0; ch < N_CH-1; ch++){
        free(i_k[ch]);
    }
    free(i_k);
    
    free(s1);
    free(s2);
    free(e_k);

    free(i_norm);
    for(int ch = 0; ch < N_CH-1; ch++){
        free(doa[ch]);
    }
    free(doa);
    
    free(i_mean);
    free(diffuseness);
    
    free(azimuth);
    free(elevation);
    free(radius);


}

//==============================================================================
void MainComponent::getComplexFFTBuffer(float** fftBuffer, size_t fftSize)
{
    
    for(int channel = 0; channel < N_CH; channel++){
        for(int sample = 0, fftsample = 0; sample < fftSize; sample++, fftsample+=2){
            complexFFTBuffer[channel][sample].real(fftBuffer[channel][fftsample]);
            complexFFTBuffer[channel][sample].imag(fftBuffer[channel][fftsample+1]);
            //std::cout << "complex: " << complexFFTBuffer[channel][sample] << std::endl;
        }
    }
    
}

void MainComponent::getPressure(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    
    std::memcpy(p_k, complexFFTBuffer[0], sizeof(std::complex<float>)*numSamples);
    
    /*for(int sample = 0; sample<numSamples; sample++){
        std::cout << "p_k: " << p_k[sample] << std::endl;
    }*/
    
}

void MainComponent::getVelocityVector(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    float scale = -1.0/(MathConstants<float>::sqrt2 * p0 * SOUND_SPEED);
    
    
    
    for(int i = 1; i < N_CH; i++ ){
        std::memcpy(u_k[i-1], complexFFTBuffer[i], sizeof(std::complex<float>)*numSamples);
        
        for(int j = 0; j < numSamples; j++){
            u_k[i-1][j] =  scale * u_k[i-1][j];
            //std::cout << "u_k: " << u_k[i-1][j] << std::endl;
        }
    }
    
}

void MainComponent::getIntensityVector(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    getPressure(complexFFTBuffer, numSamples);
    getVelocityVector(complexFFTBuffer, numSamples);
    
    for(int i = 0; i < N_CH-1; i++ ){
        for(int j = 0; j < numSamples; j++){
            std::complex<float> pk_uk = p_k[j] * std::conj(u_k[i][j]);
            i_k[i][j] = 0.5 * pk_uk.real();
            
            //std::cout << "i_k: " << i_k[i][j] << std::endl;

        }
    }
    
}

void MainComponent::getEnergyVector(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    getPressure(complexFFTBuffer, numSamples);
    getVelocityVector(complexFFTBuffer, numSamples);
    
    for(int sample = 0; sample < numSamples; sample++){
        float u_norm = std::sqrt(u_k[0][sample] * u_k[0][sample] + u_k[1][sample] * u_k[1][sample] + u_k[2][sample] * u_k[2][sample]).real();
        s1[sample] = std::pow(u_norm, 2);
        
        s2[sample] = std::pow(std::abs(p_k[sample]), 2);
        
        e_k[sample] = ((p0/4.0) * s1[sample]) + ((1.0/(4 * p0 * std::pow(SOUND_SPEED, 2))) * s2[sample]);
    }
    
    
}

void MainComponent::getDOA(std::complex<float>** complexFFTBuffer, size_t numSamples)
{
    getIntensityVector(complexFFTBuffer, numSamples);
    
    
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
        
        //std::cout << "z: " << doa[2][sample] << std::endl;
        //std::cout << "x_y: " << std::sqrt(x_y) << std::endl;
        elevation[sample] = std::atan2(doa[2][sample], std::sqrt(x_y));
    }
    
}

void MainComponent::getDifuseness(std::complex<float>** complexFFTBuffer, size_t numSamples, int dt)
{
    getIntensityVector(complexFFTBuffer, numSamples);
    getEnergyVector(complexFFTBuffer, numSamples);
    
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
void MainComponent::writeXMLData (float azimuth, float elevation) {
    
    XmlElement* result = new XmlElement("Result");
    result->setAttribute("Azimuth", azimuth);
    result->setAttribute("Elevation", elevation);
    data.addChildElement(result);
}

//==============================================================================
void MainComponent::timerCallback() {
    
    if(nextBlockReady){
        paintDoA();
        nextBlockReady = false;
        repaint();
    }
        
}

void MainComponent::paintDoA() {
    
}

//==============================================================================
void MainComponent::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    // You can add your drawing code here!
    Rectangle<float> rectArea (visualizer_W, visualizer_H, azi_resolution*2, ele_resolution*2);
    g.drawRect (rectArea, 2.0f);
    
    Line<float> vertLine (visualizer_W, visualizer_H + ele_resolution, visualizer_W + azi_resolution*2, visualizer_H + ele_resolution);
    g.drawLine(vertLine, 0.5f);
    
    Line<float> horLine (visualizer_W + azi_resolution, visualizer_H ,visualizer_W + azi_resolution , visualizer_H + ele_resolution*2);
    g.drawLine(horLine, 0.5f);

    
    Rectangle<float> pointArea (5, 5);
    
    if(resultAzi > MathConstants<float>::pi){
        point_x = visualizer_W + (resultAzi - MathConstants<float>::twoPi)*100 + azi_resolution;
    }
    else{
        point_x = visualizer_W + (resultAzi)*100 + azi_resolution;
    }
    
    point_y = visualizer_H + (-1)*(resultEle*100) + ele_resolution;
    
    
    //std::cout << "x: " << point_x << std::endl;
    //std::cout << "y: " << point_y << std::endl;
    Point<float> pointDoA (point_x, point_y);
    pointArea.setCentre (pointDoA);
    
    g.setColour (Colours::cornflowerblue);
    g.fillRect (pointArea);

}

void MainComponent::resized()
{
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.
}
