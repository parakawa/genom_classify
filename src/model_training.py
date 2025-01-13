import os
import joblib
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import accuracy_score

def evaluate_and_save_model(df, chunk_size, k, rf_params, overlap, model_dir, test_gene_dir):
    """
    train and save a RandomForest model for a specific parameter combination,
    and also save the test genes in their raw form.
    """
    from .sequence_operations import split_gens_over_all_genomes
    from .kmer_signatures import signature_for_all_genes

    print(f"training model for chunk_size={chunk_size}, k={k}...")

    # extract genes from genomes
    X_raw, y = split_gens_over_all_genomes(df, chunk_size=chunk_size, overlap=overlap)
    print(f"generated genes: {len(X_raw)}")

    # sign genes with k-mers
    X_signature = signature_for_all_genes(X_raw, k)
    print("generated signatures.")

    # split into training and testing sets
    X_train_sig, X_test_sig, y_train, y_test, X_train_raw, X_test_raw = train_test_split(
        X_signature, y, X_raw, test_size=0.2, random_state=42
    )

    # train the RandomForest model
    rf_model = RandomForestClassifier(**rf_params)
    rf_model.fit(X_train_sig, y_train)

    # evaluation with 5-fold cross-validation
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)  # 5-fold cross-validation
    cv_scores = cross_val_score(rf_model, X_train_sig, y_train, cv=cv, scoring='accuracy')
    
    # mean accuracy from cross-validation
    cv_mean_accuracy = np.mean(cv_scores)
    cv_std_accuracy = np.std(cv_scores)
    
    # evaluate the model on the test set
    y_pred = rf_model.predict(X_test_sig)
    test_accuracy = accuracy_score(y_test, y_pred)
    print(f"test set accuracy: {test_accuracy}")
    print(f"cross-validation mean accuracy: {cv_mean_accuracy} ± {cv_std_accuracy}")

    # save the model
    model_path = os.path.join(model_dir, f"model_chunk{chunk_size}_k{k}.joblib")
    joblib.dump(rf_model, model_path)
    print(f"model saved to {model_path}")

    # select random test genes and save them
    selected_indices = np.random.choice(len(X_test_raw), size=10, replace=False)
    selected_genes = [X_test_raw[i] for i in selected_indices]
    selected_labels = [y_test[i] for i in selected_indices]
    
    test_genes_df = pd.DataFrame({
        "gene": selected_genes,
        "label": selected_labels
    })
    test_gene_path = os.path.join(test_gene_dir, f"test_genes_chunk{chunk_size}_k{k}.csv")
    test_genes_df.to_csv(test_gene_path, index=False)
    print(f"test genes saved to {test_gene_path}")

    # return the metrics to save in the report
    return {
        "chunk_size": chunk_size,
        "k": k,
        "test_accuracy": test_accuracy,
        "cv_mean_accuracy": cv_mean_accuracy,
        "cv_std_accuracy": cv_std_accuracy
    }
