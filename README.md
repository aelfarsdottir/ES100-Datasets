## Code

The python code processes the sensor data into training and validation sets. 
The MATLAB code then uses those training and validation sets to develop and evaluate models.

## Preprocessed data
Using the python code, individual sensor files were resampled, interpolated, filled, and concatenated to form these preprocessed sets for training and validation of models in MATLAB.

##### DATASET4
Dataset4 was used to develop second iteration models that were eventually used in MPC trials. The fullrow.csv file (50 Mb) is available via Google Drive, as well as the MATLAB workspace structures containing N4SID identified models, and tabulated results of fit to estimation and validation data.

fullrow.csv: https://drive.google.com/open?id=1T7BlYssDevQvOC4LX8SQvB7P3z9-h01K
workspace_031319_essentials.mat: https://drive.google.com/open?id=1WQeSVsohE-RZ6Zgh-qnToohKr_5el82F

##### ALLDATA
All data includes data from August 28, 2018 to March 20, 2019. The preprocessed 1 minute and 5 second interpolations are too large (71 Mb and 838 Mb respectively) to upload, and are shared via Google Drive instead.

alldata_interp1m.csv: https://drive.google.com/open?id=1e-QuYJjgpj2YBUQusM8oF1NkYlpnDJSt
alldata_interp5s.csv: https://drive.google.com/open?id=15afoVYa50Q930uU2TrUMe8beEolen0e7
workspace_032219_ALLDATA2.mat: https://drive.google.com/open?id=1J0LB3c7PPuanKzsbyqyDubOezQsQPdxb
