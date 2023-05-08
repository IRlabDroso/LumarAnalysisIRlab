# LumarAnalysisIRlab
## Installation

```R
# Install via devtools package to clone the repository into R session 

devtools::install_github("IRlabDroso/LumarAnalysisIRlab")
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
odor_1$pulse_csv = z_score_calculation(odor_1$pulse_csv)

# Generate all the plots
PlotTrace(odor_1,combined=F,z_score=T)
PlotTraceCondition(odor_1,groupBy="pulse",z_score=T)
PlotResume(odor_1,groupBy="antenna",z_score=T)
PlotResume(odor_1,groupBy="",z_score=T)

#### Run all at once ####
setwd("path/to/your/exp_info.csv+raw.csv")
pipeline(z_score=T)
```
## Testing 

For testing purposing we included mock data to make run the package and understand the outputs. The data is accessible in the `data/.` folder. The data contain a fake experiment 
where 3 odorant (Cider, Acetophenone, Hexanal) as been tested on 3 different olfactory receptors (OR1, OR2, OR3). To test the package with this data this how to proceed: Firstly download both `.csv` files. Then execute the next script: 
```R
csv_split = Split_CSV(Path


```


