#!/usr/bin/env python3 -u

###################### To run this script ##################################

# python scripts/lassoCV_SS.py -e embeddings/esm2_650m/nonviral/BLAT_ECOLX_Ranganathan2015.pt -m data/nonviral/metadata/BLAT_ECOLX_Ranganathan2015.csv -o experiments/lassoCV_SS/esm2_650m/nonviral/BLAT_ECOLX_Ranganathan2015.csv    
    

################ imports #####################
import os
import torch
import hashlib
import argparse
import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.linear_model import LassoCV
from sklearn.preprocessing import StandardScaler
from scipy.stats import spearmanr
import random

# to ignore the convergence warnings and Rho computeation warnings. 
# Use with caution. Only use when you are sure that the model is working fine.
import warnings
warnings.filterwarnings('ignore') 
from sklearn.exceptions import ConvergenceWarning


###################### Define Functions #######################
def make_seed(dataset_id, rep):
    s = f"{dataset_id}_{rep}"
    h = hashlib.md5(s.encode()).hexdigest()[:8]
    return int(h, 16)


def split_data(df, seed, test_size=0.2):
    """
    This function site-stratified splits a dataframe into train and validation data given a seed by mutation site.
    """
    # Extract unique sites and shuffle
    df = df.copy()
    df["site"] = df["mutant"].str.extract(r'(\d+)').astype(int)
    sites = df["site"].unique().tolist()
    rng = np.random.default_rng(seed)
    rng.shuffle(sites)
    
    # Populates validation set up to the fraction corresponding to 20%
    cumsum_rows = 0
    target_val_rows = len(df) * test_size
    val_sites = []
    for site in sites:
        site_row_count = (df["site"] == site).sum()
        if cumsum_rows < target_val_rows:
            val_sites.append(site)
            cumsum_rows += site_row_count
        else:
            break
    
    # Split dataframe
    val_df = df[df["site"].isin(val_sites)]
    train_df = df[~df["site"].isin(val_sites)]
    # Sanity check
    # if not any(site in train_df['site'].unique() for site in val_df['site'].unique()):
    #     print('All sites are unique for train or test.')
        
    return train_df, val_df




def data_prep(path_compressed_embed_file, path_meta_data):
    '''Run regression on compressed embeddings'''
    scaler = StandardScaler()
    meta_data = pd.read_csv(path_meta_data)
    meta_data['target'] = scaler.fit_transform(meta_data['target'].to_frame()).squeeze()
    meta_data = meta_data.query("mutant != 'WT'")

    # load and merge the data with features
    embed = torch.load(path_compressed_embed_file, weights_only=True)
    embed_df = pd.DataFrame.from_dict(embed).T.reset_index()
    embed_df.rename(columns={'index': 'ID'}, inplace=True)

    return  meta_data, embed_df



def run_regression(meta_data, embed_df, ds_name):
    '''this version computes y_pred for train and test sets'''
    # Initialize lists for storing results
    folds, num_nonzero_coefs = [], []
    r2s_train, maes_train, rmses_train = [], [], []
    r2s_test, maes_test, rmses_test = [], [], []
    rhos_train, rhos_test = [], []

    for fold, rep in enumerate([1, 2, 3, 4, 5], start=1):
        seed = make_seed(ds_name, rep)
        
        train, test = split_data(meta_data, seed)
        
        train_data = train.merge(embed_df, how='inner', left_on='ID', right_on='ID')
        test_data = test.merge(embed_df, how='inner', left_on='ID', right_on='ID')

        y_train = train_data['target']
        y_test = test_data['target']

        X_train = train_data.iloc[:, train.shape[1]:]
        X_test = test_data.iloc[:, test.shape[1]:]

    
        # Define and train the regression model
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=ConvergenceWarning)
            model = LassoCV(max_iter=1000, n_jobs=-1)
            model.fit(X_train, y_train)

            # get the number of non-zero coefficients
            coeficients = model.coef_
            num_nonzero_coef = np.sum(coeficients != 0)

            # Make predictions
            y_pred_train = pd.DataFrame(model.predict(X_train))
            y_pred_test = pd.DataFrame(model.predict(X_test))

            # Evaluate the model
            r2_train = metrics.r2_score(y_train, y_pred_train)
            mae_train = metrics.mean_absolute_error(y_train, y_pred_train)
            mse_train = metrics.mean_squared_error(y_train, y_pred_train)
            rmse_train = np.sqrt(mse_train)
            rho_train, p_value_train = spearmanr(y_train, y_pred_train)

            r2_test = metrics.r2_score(y_test, y_pred_test)
            mae_test = metrics.mean_absolute_error(y_test, y_pred_test)
            mse_test = metrics.mean_squared_error(y_test, y_pred_test)
            rmse_test = np.sqrt(mse_test)
            rho_test, p_value_test = spearmanr(y_test, y_pred_test)

            # Append results
            r2s_train.append(r2_train)
            maes_train.append(mae_train)
            rmses_train.append(rmse_train)

            r2s_test.append(r2_test)
            maes_test.append(mae_test)
            rmses_test.append(rmse_test)

            rhos_train.append(rho_train)
            rhos_test.append(rho_test)

            folds.append(fold)
            num_nonzero_coefs.append(num_nonzero_coef)

            # Return the collected results
            print(f"Results:  fold {fold}, Seed: {seed}, r2_train: {r2_train:.3f}, r2_test: {r2_test:.3f}, Num coefs: {num_nonzero_coef}")
    print(f"Results:  r2_train: {np.mean(r2s_train):.2f}, r2_test: {np.mean(r2s_test):.2f}, Num coefs: {np.mean(num_nonzero_coefs):.2f}")
    return r2s_train, maes_train, rmses_train, r2s_test, maes_test, rmses_test, rhos_train, rhos_test, folds, num_nonzero_coefs


def save_results(r2s_train, maes_train, rmses_train, r2s_test, maes_test, rmses_test, rhos_train, rhos_test, folds, num_nonzero_coefs):
    # Create dictionary for results
    res_dict = {
        #"Model": ['Lasso'] * 5,
        "Fold": folds,
        "R2_score_train": r2s_train,
        "MAE_score_train": maes_train,
        "RMSE_score_train": rmses_train,
        "R2_score_test": r2s_test,
        "MAE_score_test": maes_test,
        "RMSE_score_test": rmses_test,
        "rho_score_train": rhos_train,
        "rho_score_test": rhos_test,
        "num_nonzero_coefs": num_nonzero_coefs
    }

    # Convert results to DataFrame
    results = pd.DataFrame(res_dict).reset_index(drop=True)
    return results





############################# Run Predictions #############################

def main():
    parser = argparse.ArgumentParser(description="Run regression for different target datasets and layers")
    parser.add_argument("-e", "--embed", type=str, help="Path to the input file")
    parser.add_argument("-m", "--metadata", type=str, help="Target name in the metadata")
    parser.add_argument("-o", "--output", type=str, help="Path to the output file")
    args = parser.parse_args()
    
    # Define the target name and output file
    path_compressed_embed_file = args.embed
    path_meta_data = args.metadata
    output = args.output

    output_dir = os.path.dirname(output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
   
    # dataset name to create unique id with a hash function
    ds_name = output.split('/')[-1].split('.csv')[0] 
    
    print(f'Loading dataset {ds_name}!')
    meta_data, embed_df = data_prep(path_compressed_embed_file, path_meta_data)

    print('Fitting the model!')
    r2s_train, maes_train, rmses_train, r2s_test, maes_test, rmses_test, rhos_train, rhos_test, folds, num_nonzero_coefs = run_regression(meta_data, embed_df, ds_name)
    results = save_results(r2s_train, maes_train, rmses_train, r2s_test, maes_test, rmses_test, rhos_train, rhos_test, folds, num_nonzero_coefs)
    results.to_csv(output)
    print(f'Process Finished!')

   
if __name__ == "__main__":
    main()