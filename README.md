# Script for the CHI 2023 paper entitled:
## Are You Human?  Investigating the Perceptions and Evaluations of Virtual Versus Human Instagram Influencers

From the authors: Anika Nissen, Colin Conrad, and Aaron Newman

The raw EEG and self-reported data to this paper can be downloaded from the following URL:
<url>TODO INPUT URL HERE</url>

Through this link, you can download the data collected in this study:
* `Study 1 Raw Data.csv` will provide you the raw collected data in Study 1
* `Study 2/Raw Behavioral Data` containts the logfiles with the responses of the participants created in the PsychoPy Paradigm
* `Study 2/total_behavioral_data_new.csv` provides the cleaned and aggregated self-reported data from the logfiles
* `Study 2/Raw EEG Data` contains the raw EEG data collected through AsaLab and converted with Matlab to the eeglab data format
* `Study 2/preprocessing logs` contains the automatically excluded ICA components from the batch script
* `Study 2/Epochs`contains the created epochs from the batch preprocessing script
* `Study 2/EEG-Results` contains the completely preprocessed EEG data that is calculated as means and peaks and exported as csv data to be analyzed in R


In this repository, you can find the figure and jupyter notebooks used to run the data analyses from the paper.
Please find short descriptions of each study what each of the scripts does here:

### Behavioral Data: Study 1 (Folder Study 1)
*Description:* Collected self-reported data through an online questionnaire created in qualtrics (here called `Survey`).

* Have a look at the questions asked and images shown in the survey, as well as the logic used in the survey: `01 Survey Virtual Influencers in PDF` and `01 Survey Virtual Influencers in Qualtrics-format`
* Look at the computed results from the collected data: `03 Computed Results in jamovi` if you want to open it directly in jamovi or `03 Computed Results in PDF` if you just want to get a more detailed look of the results


### Behavioral and EEG Data: Study 2 (Folder Study 2)
*Description:* In study 2, we first collected demographics and Instagram use behavior with a questionnaire (here called `Survey EEG`). After that, the EEG preparation took place. In a booth, participants were shown an event-related paradigm (here called `Paradigm EEG`)

Study Design:
* Have a look at the questions asked and images shown in the survey, as well as the logic used in the survey: `01 Survey EEG Instagram Influencers in PDF` and `01 Survey EEG Instagram Influencers in Qualtrics-format`
* Have a look at the ERP Paradigm in which EEG was measured. Files of the PsychoPy Experiment can be found in the folder `02 Paradigm EEG ERP`

Self-Reported Data:
* To get insights into the demographics of the recruited sample for this study, please take a look at `03 EEG-Demographics` to see them in a PDF file, or open the `03 EEG-Demographics in jamovi` file in jamovi.
* The group analyses of the self-reported data with Mixed-Effects Models can be seen in the `04 stats_Instagram_Behavioral_in_R` notebook

EEG Data:
* EEG data was preprocessed with a batch script and subsequent individual inspection and correction of the ICA of the batch script in `05 EEG-Batch-Preprocessing`
* After preprocessing, we took a look at the grand averages across subjects in the `06 stats_Instagram_grand_average_analysis` notebook
* To run further mixed-effects models for the ERP components of interest (N400, LPP), we have to prepare the EEG data first in `07 stats_Instagram_eeg_preparation`
* Finally, we rand mixed-effects models for just the EEG, as well as with the self-reported data integrated in the scripts `08 stats_Instagram_LPP_analysis_vis` for the LPP and in the `08 stats_Instagram_N400_analysis_vis` for the N400 component
* `NCIL_functions.R`is a helper script for the analyses in the prior step


