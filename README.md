# Kaggle DSB 2016

Solution based on R and primarily using EBImage for image processing. Developed on Mac OS
and Windows. No other 3rd p required except required R packages.

See 
* https://www.kaggle.com/c/second-annual-data-science-bowl/model
* https://www.kaggle.com/wiki/ModelSubmissionBestPractices
* https://github.com/perdoperdo/dsb
* https://www.blogger.com/blogger.g?blogID=3288168447132353117#allpages


## Initialization

Run playlist.R to read all image files and meta information, creating imagelist.csv, which drives
all subsequent steps. Since this reads all images for meta info, it takes a few hours to run for the
train/validate data.

## Segmentation

Run extract.R to segment all images from both train and validate/test dataset and create two files
segments-<dataset>.csv with details of all segments of all images. Just for display purpuses, the
segmented images are also saved in a folder "segmented". This is a very slow process and can take
a couple of days for the 700 images in the train/validate set.

## Classification

Run "classify.R" to train the system interactively on what the LV segments are. The results are saved
in a file "segments-predict.csv". This only needs to be done (and can only be done) on the train set.
After having classified the first few dozen images, the system will start to build models and highlight
(in green) the predicted segments.

## Predictions

Run "predict.R" to predict the systole and diastole volumes on the validate/test dataset. This involves
imputing missing data and outlier detection and can take a couple of hours.

## Feedback cycle

Models for segment prediction - if available - will be used in the segmentation process. Therefore it is
useful to re-do the segmentation (and prediction) after doing a first batch of manual classification.

