import os
import argparse
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

import torch
from torch.optim import AdamW
from torchmetrics.regression import R2Score, SpearmanCorrCoef
from torch.utils.data import Dataset, DataLoader
torch.set_float32_matmul_precision('medium')

import pytorch_lightning as pl
from pytorch_lightning.loggers import CSVLogger
from pytorch_lightning.callbacks import EarlyStopping

from scripts.ESM2ne_partial import LoadFromPretrained



#python scripts/LitESM2ne_partial_trainer.py -i data/cellular/metadata/BLAT_ECOLX_Ranganathan2015.csv -o results/fineTune/test/partial/
#python scripts/LitESM2ne_partial_trainer.py -i data/cellular/metadata/YAP1_HUMAN_Fields2012_singles.csv -o results/fineTune/test/partial/

# larger sequences
# CVB3_POLG_Mattenberger2185
#python scripts/LitESM2ne_partial_trainer.py -i data/viral/metadata/CVB3_POLG_Mattenberger.csv -o results/fineTune/test/partial/ 
# BRCA1_HUMAN_RING	1863
#python scripts/LitESM2ne_partial_trainer.py -i data/viral/metadata/BRCA1_HUMAN_RING.csv -o results/fineTune/test/partial/ 



def parse_args():
    parser = argparse.ArgumentParser(description="Training script arguments.")
    parser.add_argument("-i", "--data", type=str, required=True, help="Path to input CSV file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file path.")
    parser.add_argument("--checkpoint", type=str, default="esm2_t33_650M_UR50D", help="Model checkpoint name or full path.")
    parser.add_argument("--num_classes", type=int, default=1, help="Number of classes (1 for regression).")
    parser.add_argument("--epochs", type=int, default=50, help="Number of epochs.")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size.")
    parser.add_argument("--device", type=str, default=("cuda" if torch.cuda.is_available() else "cpu"), help="Device for training.")
    parser.add_argument("--learning_rate", type=float, default=1e-4, help="Learning rate.")
    parser.add_argument("--weight_decay", type=float, default=1e-2, help="Weight decay.")
    parser.add_argument("--dropout", type=float, default=0.2, help="CLS dropout probability.")
    parser.add_argument("--split_strategy", type=str, default="pool_split", help="Split strategy to use.")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")
    return parser.parse_args()



class MyDataset(Dataset):
    """This class just loads the data and return a dataset object, returning the sequences and targets.
    Without any tokenization or padding.
    This will be handled later in the collate_fn.
    """
    def __init__(self, data_file):
        self.scaler = StandardScaler()
        self.data = data_file
        self.sequences = [(d['ID'], d['sequence']) for _, d in self.data.iterrows()]
        self.targets = self.scaler.fit_transform(self.data['target'].to_frame()).squeeze()
     
    def __len__(self):
        return len(self.targets)

    def __getitem__(self, idx):
        sequences = self.sequences[idx]
        target = torch.tensor(self.targets[idx], dtype=torch.float).unsqueeze(-1)
        
        return sequences, target


def split_data(df, seed, test_size=0.2):
    """
    This function randomly splits a dataframe into train, validation data given a seed by mutation site.
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




class MyDataModule(pl.LightningDataModule):
    def __init__(self, tokenizer, args):
        super().__init__()
        self.args = args
        self.alphabet = tokenizer
        self.split_strategy = args.split_strategy # 'pool_split' or 'site_split'
        
    def setup(self, stage=None):
        data = pd.read_csv(args.data)
        data = data.query("mutant != 'WT'")

        if self.split_strategy == 'site_split':
            train_df, val_df = split_data(data, seed=args.seed, test_size=0.2)
        elif self.split_strategy == 'pool_split':
            train_df, val_df = train_test_split(data, random_state=args.seed, test_size=0.2)
  
        self.train_dataset = MyDataset(train_df)
        self.val_dataset = MyDataset(val_df)
       
    def collate_fn(self, batch):
        """This function will be used to collate the data into a batch.
        It will handle the tokenization and padding of the sequences.
        """
        # batch: a list of ( (ID, sequence), target )
        seqs = [item[0] for item in batch]
        targets = [item[1] for item in batch]
    
        # Tokenize the sequences
        batch_converter = self.alphabet.get_batch_converter()
        batch_labels, batch_strs, batch_tokens = batch_converter(seqs)
        #batch_lens = (batch_tokens != self.alphabet.padding_idx).sum(1)
        targets = torch.stack(targets) # Stack the targets to a shape (batch_size, 1)
        return batch_tokens, targets

    def train_dataloader(self):
        return DataLoader(
            self.train_dataset,
            batch_size=self.args.batch_size,
            shuffle=True,
            collate_fn=self.collate_fn,
            num_workers=4,
            pin_memory=True,
            drop_last=True,
        )

    def val_dataloader(self):
        return DataLoader(
            self.val_dataset,
            batch_size=self.args.batch_size,
            collate_fn=self.collate_fn,
            num_workers=4,
            pin_memory=True,
            drop_last=True,
        )


################ pytorch lightning model ######################
class LitModel(pl.LightningModule):
    def __init__(self, args):
        super().__init__()
        self.args = args
        self.lr = args.learning_rate
        self.weight_decay = args.weight_decay
        
        # Define the model
        self.model_loader = LoadFromPretrained(
            checkpoint=args.checkpoint, num_classes=args.num_classes, hidden_dropout=args.dropout)
        self.model, self.alphabet = self.model_loader.get_model_details()

        # metrics and loss function
        self.loss_fn = torch.nn.MSELoss()
        self.train_r2 = R2Score()
        self.val_r2 = R2Score()

        self.train_rho = SpearmanCorrCoef()
        self.val_rho = SpearmanCorrCoef()
    
    def forward(self, batch_tokens):
        outputs = self.model(batch_tokens) # type: ignore
        return outputs['logits'] 


    def training_step(self, batch):
        """This function will be called for each batch during training.
        It will compute the loss and log it.
        """
        batch_tokens, targets = batch
        preds = self(batch_tokens)
        loss = self.loss_fn(preds, targets)
        self.train_r2(preds, targets)
        self.train_rho(preds, targets)
        # Log the loss by epoch
        self.log("train_loss", loss, on_step=False, on_epoch=True, prog_bar=True, logger=True, sync_dist=True)
        self.log("train_r2", self.train_r2, on_step=False, on_epoch=True, prog_bar=True, logger=True, sync_dist=True)
        self.log("train_rho", self.train_rho, on_step=False, on_epoch=True, prog_bar=True, logger=True, sync_dist=True)
        return loss
 
    def validation_step(self, batch):
        batch_tokens, targets = batch
        preds = self(batch_tokens)
        loss = self.loss_fn(preds, targets)
        self.val_r2(preds, targets)
        self.val_rho(preds, targets)
        self.log("val_loss", loss, on_step=False, on_epoch=True, prog_bar=True, logger=True, sync_dist=True)
        self.log("val_r2", self.val_r2, on_step=False, on_epoch=True, prog_bar=True, logger=True, sync_dist=True)
        self.log("val_rho", self.val_rho, on_step=False, on_epoch=True, prog_bar=True, logger=True, sync_dist=True)
        return loss # this loss is not used, but I could return something else and modify.

    def configure_optimizers(self):
        #return AdamW((p for p in self.parameters() if p.requires_grad), lr=self.lr, weight_decay=self.weight_decay)
        return AdamW(self.parameters(), lr=self.lr, weight_decay=self.weight_decay)




def main():
    args = parse_args()
    dataset_name = args.data.split("/")[-1].split(".csv")[0]
    print(f"Model checkpoint: {args.checkpoint}")
    print(f"Dataset name: {dataset_name}")

    
    # ########################## Training setup ##########################
    logger = CSVLogger(
            save_dir=args.output,
            name=f"{dataset_name}",
            version=args.seed)
        
    trainer = pl.Trainer(
        accelerator="gpu",
        #strategy="ddp",
        devices=[1],
        max_epochs= args.epochs,
        enable_checkpointing=False,
        logger=logger,
        callbacks=[EarlyStopping(monitor="val_loss", patience=10, mode="min")],
    )

    print("Loading model and dataset...")
    model = LitModel(args)
    datamodule = MyDataModule(model.alphabet, args)

    # Train
    trainer.fit(model, datamodule)

  
  
  

if __name__ == '__main__':
    args = parse_args()
    main()
