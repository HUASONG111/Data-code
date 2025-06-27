# Data-code

## 1. Project Overview

This project investigates the lagged response of sea surface chlorophyll-a (Chla) to environmental factors.

* **R (v4.4.1)** is primarily used to run statistical models and generate quantitative results.
* **Python (v3.12.7)** is used to visualize annual mean concentrations and lag timing results.
* **MATLAB (v24.2)** is used to generate final visualizations, particularly for latitudinal distributions.

---

## 2. System Requirements

* Operating System: Windows/macOS/Linux (tested on macOS 14.6)
* Python 3.12.7 with necessary packages (`numpy`, `pandas`, `matplotlib`, `seaborn`, `notebook`)
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
    install.packages(c("mgcv", "dlnm", "mvmeta", "data.table", "doParallel"))
    ```

* **MATLAB**:

  * Install MATLAB 2024a (R24.2)
  * No additional toolboxes required beyond the default installation

---

## 4. Demo

A small simulated dataset is provided in the `data/` folder.

### Execution Order:

**1. R Scripts (Modeling and Analysis)**
Run the following R scripts in order:

* `Find lag onset time and duration.R` – identify lag timing and duration
* `Lagged response across four Longhurst biomes. .R` – model across biomes
* `Forest plot.R` – generate forest plot
* `Meta analysis.R` – perform meta-analysis
* `Games-Howell test.R` – statistical testing (Games-Howell)

**2. Python Notebooks (Visualization)**
Open and run the following `.ipynb` notebooks in VSCode or Jupyter:

* `Global distributions of annual mean factors.ipynb` – plot annual mean maps
* `Visualization of the lag onset time and duration.ipynb` – visualize lag patterns

**3. MATLAB Script (Final Visualization)**
Run the following MATLAB script:

* `Latitudinal_distributions.m` – visualize latitudinal distributions of lag metrics

### Expected Runtime:

* R scripts: \~5 minutes to **up to 12 hours**, depending on the dataset size and spatial extent (e.g., global scale)
* Python notebooks: \~2–3 minutes
* MATLAB: \~1 minute

---

## 5. Instructions for Use

To apply the workflow to your own data:

1. Replace the simulated dataset in `data/` with your own data (matching format is required)
2. Adjust file paths inside each script if needed
3. Follow the prescribed execution order:

   * First run all **R scripts** to generate intermediate modeling outputs
   * Then run **Python notebooks** to visualize annual means and lag information
   * Finally, run the **MATLAB script** to generate latitudinal distribution figures

---

## 6. Reproducibility

All quantitative results in the associated manuscript can be fully reproduced using this workflow and the sample dataset provided in the `data/` folder.

* Ensure consistent software versions (R 4.4.1, Python 3.12.7, MATLAB R2024a)
* Use fixed seeds in R for reproducibility of parallel models
* Outputs from each stage are saved and reused by subsequent scripts

---

## 7. Contact

For questions, suggestions, or technical issues, please contact:
**\[Your Name]**
Email: **[your.email@example.com](mailto:your.email@example.com)**
