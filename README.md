# Analyzing-calcium-imaging-in-MATLAB
<p align="center"> <b> Background information </b> </p>

Most people are comfortable with the idea of rcording from the brain using electrical activity; however, there is a second method of recording from the brain that takes advantage of changes in intracellular calcium that occur when a neuron is active. In recent years scientists have developed a family of proteins, known as calcium indicators, that enable optical recording of individual neurons. As the name suggests, a calcium indicator is sensitive to the amount of calcium that is present in a cell expressing the protein. When the concentration of calcium inside of the cell increases the fluorescence activity of the cell also increases, and when the concentration of calcium decreases so too does the fluorescence. This enables researchers like me to image the fluorescence activity of neurons near the surface of the brain. Below is the average projection of one of these videos. The black lines are blood vessels in the field of view, and the round, white circles are individual neurons. 

![Example field of view](https://github.com/RedingtonE/Analyzing-calcium-imaging-in-MATLAB/blob/master/data/Blackandwhitefieldofview.png)

Before analyzing data from these videos I correct for movement of the brain using registration algorithms and reduce background noise using bandpass and homomorphic filtering. Once the videos have been sufficiently cleaned up I can identify neurons in the field of view using constrained non-negative matrix factorization (CNMF). Below is an example field of view in black and white with the individual neurons identified using CNMF overlaid in green. 

![Example field of view](https://github.com/RedingtonE/Analyzing-calcium-imaging-in-MATLAB/blob/master/data/examplefieldofview.png)

<p align="center"> <b> Extracting fluorescence activity in this repository </b> </p>

The code in this repository represents a small piece of the workflow that's used to generate data from videos and then analyze it to generate novel insights into how the brain functions. This code assumes that the input video has already been registered and filtered, and that a matrix of input neurons (masks) has already been generated. Using this input the function genDfF will isolate average fluorescence and background traces. These are then used to normalize the fluorescence into the deltaF/F. An example of the output calcium activity is provided in the data folder in the file 'exampledata.mat'. Below is an example deltaF/F trace: 

![Example calcium trace](https://github.com/RedingtonE/Analyzing-calcium-imaging-in-MATLAB/blob/master/data/examplecalciumtransients_codeoutput.png)

This trace shows the typical characteristics of a deltaF/F trace. First, the baseline is noisy compared to the signal. This is because the emission of light is a poisson process. As a result, the dim signals acquired in calcium imaging will have low SNR. Second, the trace contains calcium transients, which are characterized by a sudden increase in fluorescence intensity followed by a slow, exponential decay to baseline. An efficient way of identifying calcium transients is estimating the background of the signal, and then setting a threshold. The function peakLocation.m will set a threshold using the median absolute deviation of the deltaF/F trace, and then identify calcium transients that meet certain temporal criteria. 

![Thresholding example](https://github.com/RedingtonE/Analyzing-calcium-imaging-in-MATLAB/blob/master/data/examplecalciumtransients_codeoutput_threshold.png)

In the above figure an orange line has been drawn to show the threshold that is used to select calcium transients. 

![Thresholding example](https://github.com/RedingtonE/Analyzing-calcium-imaging-in-MATLAB/blob/master/data/examplecalciumtransients_codeoutput_threshold_identifiedPeaks_andtroughs.png)

Finally the peaks and troughs of the calcium transients are identified and saved for future work. Once these calcium transients are identified they can be used to address neuroscience questions. In my repository on decoding neural activity I use these transient times to train a bayesian decoder how to decode a mouse's position in space from neural activity. 

