import os
import pandas as pd
from src.model_training import evaluate_and_save_model
from src.utils import load_dataframe

print("loading data...")
# load the data
df = load_dataframe("df_organisms_selection.pkl")

chunk_sizes = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
ks = [2, 3, 4, 5, 6, 7, 8]
# chunk_sizes = [500, 1000, 1500]
# ks = [2, 3, 4]

rf_params = {
    "n_estimators": 100,       # 100 estimators
    "max_depth": 15,           # limit the depth of trees
    "max_features": 'sqrt',    # use the square root of the number of features
    "min_samples_split": 10,   # minimum number to split a node
    "min_samples_leaf": 5,     # minimum number of samples per leaf
    "random_state": 42
}

# create directories to save the results
os.makedirs("models", exist_ok=True)
os.makedirs("test_genes", exist_ok=True)

# train the models for all parameter combinations
results = []
for chunk_size in chunk_sizes:
    for k in ks:
        metrics = evaluate_and_save_model(
            df,
            chunk_size,
            k,
            rf_params,
            overlap=0,
            model_dir="models",
            test_gene_dir="test_genes"
        )
        results.append(metrics)

# save the results to a CSV file
results_df = pd.DataFrame(results)
results_df.to_csv("grid_search_results.csv", index=False)
print("results saved")
