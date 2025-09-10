# VeritasRI: Your Solution for Clinical Reference Intervals 📊

## Say Goodbye to Manual Data Analysis.
Have you ever spent hours manually sifting through clinical data, running the same analyses over and over, or struggling to interpret complex results? We get it. That's why we built VeritasRI, an intuitive R Shiny application designed to solve these real-world problems in clinical and medical research. This isn't just another statistical tool; it's a way to transform your workflow from a chore into a precise, automated process.

VeritasRI takes the pain out of estimating clinical **reference intervals (RIs)**. By leveraging advanced statistical methods like **refineR** and **Gaussian Mixture Models (GMM)**, it helps you accurately determine reference ranges even when your data is messy, mixed, or contains hidden subpopulations.

<img width="1920" height="887" alt="1" src="https://github.com/user-attachments/assets/5c87bab8-2d16-426b-ba51-47e99faa6969" />

> "Veritas filia temporis." (Truth is the daughter of time.) — Latin Proverb

In clinical research, uncovering the truth about a patient's health depends on accurate reference intervals. **Veritas** (truth) is at the core of our mission, and VeritasRI is the tool that helps you find that truth efficiently and reliably, turning complex data into clear, actionable insights.

---

## ✨ What Can VeritasRI Do for You?

### 1. Core Analysis, Made Simple 📈
The app's main tab guides you through the essential steps of estimating RIs. Just upload your Excel file, select your value, age, and gender columns, and let the app handle the rest. You can filter by demographics, choose data transformation models, and visualize your results with clear, interactive plots that show estimated RIs and confidence intervals.

### 2. Uncover Hidden Subpopulations 🔬
Your data often holds secrets. Our **GMM Clustering** feature helps you find them. By using the powerful **mclust** package, VeritasRI can detect distinct subpopulations within your dataset based on value and age. The app then provides detailed summaries and color-coded scatter plots, so you can easily see and understand these hidden groups.

### 3. Get Results Faster with Parallel Analysis 🚀
Need to run multiple analyses for different age groups or genders? No problem. The **Parallel Analysis** tab uses the **future** package to run multiple **refineR** analyses at the same time, leveraging all your CPU cores to give you results in a fraction of the time. You can view all your results together in a combined summary or explore each group individually.

<img width="1091" height="1101" alt="Scherm_afbeelding 2025-08-26 om 16 02 16" src="https://github.com/user-attachments/assets/7664dae8-d948-452d-ac71-337a0b1fbeb5" />

---

## 🧠 The Brains Behind the App

We've integrated cutting-edge methods to ensure you get the most reliable results:

* **refineR:** A robust statistical method that handles mixed populations, so you don't have to manually exclude "pathological" data.
* **Box-Cox & Modified Box-Cox Transformations:** These smart transformations help normalize non-Gaussian data, ensuring your analyses are valid.
* **Gaussian Mixture Models (GMM):** The go-to method for identifying distinct subpopulations in your data.
* **Parallel Computing:** The `future` package turbocharges your analysis by running multiple tasks concurrently.

---

## 🛠️ Get Started in Minutes

1.  **Install R and RStudio:** If you don't have them, download R from [CRAN](https://cran.r-project.org/) and RStudio from [posit](https://posit.co/download/rstudio-desktop/).

2.  **Install Packages:** Run this command in your R console to get all the necessary packages:
    ```R
    install.packages(c(
      "shiny", "refineR", "readxl", "tidyverse", "mclust", "moments", "shinyjs", "car",
      "shinyFiles", "shinyWidgets", "bslib", "ggplot2", "future", "future.apply", "bsicons", "remotes"
    ))
    ```

3.  **Clone & Run:** Clone the repository and run the app from RStudio.
    ```bash
    git clone [https://github.com/yakubinaweed/VeritasRI.git](https://github.com/yakubinaweed/VeritasRI.git)
    ```
    Open the project folder in RStudio, then open the `app.R` file and click the **'Run App'** button in the top-right corner.

---

## 📝 A Quick Guide to Usage

* **Analyze:** Choose your analysis tab (Main, GMM, or Parallel).
* **Upload:** Just drop your Excel file (`.xlsx`) into the app.
* **Select:** Tell the app which columns contain your values, age, and gender, and select your settings.
* **Press Analyze:** Click the button to start the analysis.
* **Review:** Dive into the plots and summaries to understand your results.

---

## 🤝 We Welcome Your Contributions!

VeritasRI is an open-source project. If you have ideas, suggestions, or bug reports, please feel free to open an issue or submit a pull request.

[**View the VeritasRI GitHub Repository**](https://github.com/yakubinaweed/VeritasRI)

---

## 📄 License

This project is licensed under the MIT License—see the `LICENSE` file for details.
