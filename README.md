# Genom Classify

## Project Overview
This project processes genomic data from bacterial and archaeal organisms retrieved from the NCBI database, cleans the data, and trains Random Forest models using varying k-mer sizes and chunk sizes. The pipeline evaluates model performance and saves results for further analysis.

## Setup Instructions

### 1. Create a Conda Environment
Create a new Conda environment named `genom_classify` with Python 3.12:
```bash
conda create -n genom_classify python=3.12 -y
```

### 2. Activate the Environment
Activate the environment:
```bash
conda activate genom_classify
```

### 3. Install Dependencies
Install the required dependencies listed in the `requirements.txt` file:
```bash
pip install -r requirements.txt
```

## Running the Project

### Step 1: Data Processing
To process the genomic data, run the following command from the root of the project:
```bash
python data_processing.py
```

This script performs the following tasks:
- Retrieves more than 100 bacterial and archaeal organisms' genomes from the NCBI database.
- Filters out genomes with invalid sequences (e.g., sequences containing letters other than A, C, G, or T).
- Randomly selects 10 organisms from the filtered data and saves them to memory for use in the modeling pipeline.

### Step 2: Model Training and Evaluation
Run the main pipeline using:
```bash
python main.py
```

This script trains Random Forest models using varying configurations of:
- **k-mer sizes**: [2, 3, 4, 5, 6, 7, 8]
- **chunk sizes**: [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]

For each combination of parameters, the pipeline:
1. Trains a Random Forest model.
2. Saves the trained model to the `models/` directory.
3. Generates test gene CSV files in the `test_genes/` directory.
4. Outputs a `grid_search_results.csv` file with accuracy metrics for each configuration.

### Recommendation
Running the full pipeline with all k-mer sizes and chunk sizes takes a significant amount of time. For initial testing, modify `main.py` to restrict:
- **k-mer sizes**: [2, 3, 4]
- **chunk sizes**: [500, 1000, 1500]

This will significantly reduce runtime while allowing you to validate the pipeline functionality.

## Project Structure
- `data_processing.py`: Handles data retrieval, filtering, and selection.
- `main.py`: Main script for model training and evaluation.
- `models/`: Directory where trained models are saved.
- `test_genes/`: Directory where test gene CSV files are saved.
- `grid_search_results.csv`: Contains accuracy metrics for each model configuration.
- `requirements.txt`: Lists all Python dependencies required for the project.

## Future Steps
Once the pipeline is validated and running efficiently, you can expand the parameter grid or experiment with different machine learning models and hyperparameters to improve performance.

