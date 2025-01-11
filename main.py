import os
import pandas as pd
from src.model_training import evaluate_and_save_model
from src.utils import load_dataframe

print("loading data...")
# load data
df = load_dataframe("df_organisms_selection.pkl")

chunk_sizes = [500, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
ks = [2, 3, 4, 5, 6, 7, 8]
rf_params = {"n_estimators": 100, "max_depth": None, "random_state": 42}

# create directories for saving outputs
os.makedirs("models", exist_ok=True)
os.makedirs("test_genes", exist_ok=True)

# train models for all parameter combinations
results = []
for chunk_size in chunk_sizes:
    for k in ks:
        accuracy = evaluate_and_save_model(
            df,
            chunk_size,
            k,
            rf_params,
            overlap=0,
            model_dir="models",
            test_gene_dir="test_genes"
        )
        results.append({"chunk_size": chunk_size, "k": k, "accuracy": accuracy})

# save grid search results
results_df = pd.DataFrame(results)
results_df.to_csv("grid_search_results.csv", index=False)
