# Reward salience but not spatial attention dominates the value representation in the orbitofrontal cortex
## Download dataset and run
Download the dataset from http://doi.org/10.5281/zenodo.7090240, unzip it and copy it to the reward_salience_in_OFC directory    
Navigate to the reward_salience_in_OFC directory and run main  
All the main figures in the paper will be reproduced
## Data description
* bhvD: task conditions in monkey D. Each cell contains task information for each recording session.  
* bhvG: task conditions in monkey G.
* ensembleDLPFC: DLPFC neuronal ensemble. 
* ensembleOFC: OFC neuronal ensemble.
* exampleDLPFC: example DLPFC neuron.
  * norm_fr: 1 by 2 cell array. Neuronal responses aligned to stimulus onset and luminance change.
  * aveStim_norm_fr: trial_num by 1 vector. Neuronal responses averaged during the whole stimulus period.
* exampleOFC: example OFC neuron.
* populationDLPFC: population DLPFC activity.
  * negIdx: negative-tuned DLPFC neurons' index.
  * posIdx: positive-tuned DLPFC neurons' index.
  * singleNeuronFrVar: neuronal responses of individual DLPFC neurons.
  * singleNeuronRegResults: linear regression results of individual DLPFC neurons.
* populationOFC: population OFC activity.
