# LumarAnalysisIRlab
## Installation

```R
# Install via devtools package to clone the repository into R session 

devtools::install_github("IRlabDroso/LumarAnalysisIRlab")

# Or with Tar.gz source file 

install.packages("Path_to_file", type ="source", repos = NULL, dependencies = TRUE)
```

## Usage

Each function of this package are depending on another one except the first in the pipeline. User can run each functions separately or all at once calling `pipeline()`.

```R
#### Running all functions separetely ####

# Split the csv 
setwd("path/to/your/exp_info.csv+raw.csv")
csv_split = Split_CSV(getwd())

# format the table 
csv_DF = format_csv(csv_input = csv_split)

# Identify the pics and overlapp them
odor_1 = Finding_pics(csv_DF$rolling_mean, exp_odorant="water 1")

# Compute z score on dataset
odor_1$pulse_csv = Z_score_calculation(odor_1$pulse_csv)

# Generate all the plots
PlotTrace(odor_1,combined=F,z_score=T)
PlotTraceCondition(odor_1,groupBy="pulse",z_score=T)
PlotResume(odor_1,groupby="antenna",z_score=T)
PlotResume(odor_1,groupby="",z_score=T)

#### Run all at once ####
setwd("path/to/your/exp_info.csv+raw.csv")
pipeline(z_score=T)
```
## Testing 

For testing purposing we included mock data to make run the package and understand the outputs. The data is accessible in the `data/.` folder or calling it directly in R with `data()` function. The data contain a fake experiment where 3 odors as been tested on 3 different olfactory receptors (OR1, OR2, OR3). To test the package with this data this how to proceed:
```R
library(LumarAnalysisIRlab)
setwd("Path/to/new/directory/")
data(Exp_infos)
data(raw_data_mock)

write.csv(exp_info,file="Exp_Infos.csv",row.names=F)
write.csv(raw,file="raw_data_mock.csv",row.names=F)

pipeline(z_score=T)

```


