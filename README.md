# Data-code

## 1. Project Overview

This project investigates the lagged response of sea surface chlorophyll-a (Chla) to environmental factors.

* **R (v4.4.1)** is primarily used to run statistical models and generate quantitative results.
* **Python (v3.12.7)** is used to visualize annual mean concentrations and lag timing results.
* **MATLAB (v24.2)** is used to generate final visualizations, particularly for latitudinal distributions.

---

## 2. System Requirements

* Operating System: Windows/macOS/Linux (tested on macOS 14.6)
* Python 3.12.7 with necessary packages (`numpy`, `pandas`, `matplotlib`, `seaborn`, `notebook`, etc.)
* R 4.4.1 with required packages (`mgcv`, `dlnm`, `mvmeta`, `data.table`, `doParallel`, etc.)
* MATLAB 2024a (R24.2) with standard visualization toolboxes
* Recommended: A multi-core CPU (number of cores should be selected based on your computer’s configuration)

---

## 3. Installation

* **Python**:

  * Install Python 3.12.7 from [python.org](https://www.python.org)
  * Install required packages with pip:

    ```bash
    pip install numpy pandas matplotlib seaborn notebook
    ```
  * Recommended: use **Visual Studio Code (VSCode)** with the **Python** and **Jupyter** extensions to run `.ipynb` files interactively.

* **R**:

  * Install R 4.4.1 from [CRAN](https://cran.r-project.org)
  * Install required packages using:

    ```r
    install.packages(c("mgcv", "dlnm", "dplyr", "ggplot2", "mvmeta", "data.table", "doParallel", "cowplot", "parallel", "readxl", "patchwork", "rstatix"))
    ```
 
* **MATLAB**:

  * Install MATLAB 2024a (R24.2)
  * No additional toolboxes required beyond the default installation

**Typical install time on a standard desktop computer:**

* Python environment: ~5–10 minutes  
* R environment: ~5–10 minutes  
* MATLAB setup: ~15–30 minutes (depending on installer and licensing, one-time)
* Times may vary depending on network speed and whether the software has been previously installed.
---

## 4. Demo

The data/ folder contains a small real-world dataset (2–3°E, 40–41°S) for demonstration. Simulated data used in statistical testing is generated internally within the R scripts.

### Execution Order:

### 1. R Scripts (Modeling and Analysis)
Run the following scripts **in order** using R:

1. **`Code1_Find lag onset time and duration.R`**  
   - **Input:** `data/1.Original data of each single point/` (includes `chl_01.csv`, ..., `chl_25.csv`)
   - **Output:** Lag timing/duration data and response plots at individual marine locations  
   - **Note:** This step generates the dataset used by downstream scripts.

2. **`Code2_1_Games-Howell test.R`**  
   - **Input:** Simulated data (generated internally)  
   - **Output:** Pairwise comparison plots between environmental drivers using Games-Howell test

3. **`Code3_Meta analysis.R`**  
   - **Input:** `data/3.Summary data/` (includes `Zone1.csv`, ..., `Zone4.csv`) 
   - **Output:** Meta-analysis lag effect plots and summary values at regional scale

4. **`Code4_Forest plot.R`**  
   - **Input:** `data/4.Forest meta data/` (filtered meta-analysis results, includes `Forest_meta_neg.xlsx` and `Forest_meta_pos.xlsx`)  
   - **Output:** Forest plots showing maximum and minimum lag effect sizes across factors

---

### 2. Python Notebooks (Visualization)
Open and run the following `.ipynb` notebooks in **JupyterLab** or **VSCode**:

1. **`Code2_2_the lag onset visualization.ipynb`**  
   - **Input:** `data/2.Single point results/` (generated from Code1_Find lag onset time and duration.R, includes `lagday_po4.csv`, ..., `lagday_ssw.csv`)  
   - **Output:** Spatial visualization of lag onset time

2. **`Code2_3_the lag duration visualization.ipynb`**  
   - **Input:** `data/2.Single point results/` (generated from Code1_Find lag onset time and duration.R, includes `duration_po4.csv`, ..., `durtion_ssw.csv`)  
   - **Output:** Spatial visualization of lag duration

3. **`Code5_Global distributions of annual mean factors.ipynb`**  
   - **Input:** `data/5.Annual average data/` (includes `Chl_data.mar`, `Lat.mat` and `Lon.mat`)  
   - **Output:** Global maps showing annual mean concentrations of selected environmental factors

---

### 3. MATLAB Script (Final Visualization)
Run the following MATLAB script:

1. **`Code2_4_Latitudinal_lag_onset_time_distribution.m`**  
   - **Input:** `data/2.Single point results/` (includes `lagday_po4.csv`, ..., `lagday_ssw.csv`) 
   - **Output:** Latitudinal visualization of mean lag onset time

2. **`Code2_5_Latitudinal_lag_duration_distribution.m`**  
   - **Input:** `data/2.Single point results/` (includes `duration_po4.csv`, ..., `durtion_ssw.csv`) 
   - **Output:** Latitudinal visualization of mean lag duration

### 4. Expected Runtime

* R scripts: \~5–8 minutes (in total)
* Python notebooks: \~2–3 minutes
* MATLAB: \~1 minute

---

## 5. Instructions for Use

To apply the workflow to your own data:

1. Replace the dataset in `data/` with your own data (matching format is required)
2. Adjust file paths inside each script if needed
3. Follow the prescribed execution order:

   * First run all **R scripts** to generate intermediate modeling outputs
   * Then run **Python notebooks** to visualize annual means and lag information
   * Finally, run the **MATLAB script** to generate latitudinal distribution figures

---

## 6. Reproducibility

This workflow is fully reproducible with the provided subset of real data included in the `data/` folder. 

- **Software versions used:**  
  - R: 4.4.1  
  - Python: 3.12.7  
  - MATLAB: R2024a

- **Data availability:**  
  A small subset of real data (from 2–3°E, 40–41°S) is shared in the `data/` folder for demonstration. The full original dataset is not publicly available. Simulated data used for some statistical tests is generated internally within the R scripts.

- **Reproducibility:**  
  All code, analysis steps, and output formats are identical and fully reproducible using the provided data subset, allowing validation of the entire analysis pipeline.

