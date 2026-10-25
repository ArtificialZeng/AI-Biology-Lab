import os
import sys
# Add project root to sys.path
#sys.path.append(os.path.join(os.path.abspath('..') , 'ESM2ne'))
from scripts.ESM2ne_Lora import Load_from_pretrained

import argparse
import numpy as np
import pandas as pd

import torch
from torch.optim import AdamW
from torchmetrics.regression import R2Score, SpearmanCorrCoef
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
torch.set_float32_matmul_precision('high')

import pytorch_lightning as pl
from pytorch_lightning.loggers import CSVLogger
from pytorch_lightning.callbacks import EarlyStopping

#python scripts/LitESM2ne_Lora_trainer.py -i data/cellular/metadata/BLAT_ECOLX_Ranganathan2015.csv -o results/fineTune/test/lora/ --checkpoint_path esm2_t33_650M_UR50D
#python scripts/LitESM2ne_Lora_trainer.py -i data/cellular/metadata/PABP_YEAST_Fields2013_singles.csv -o results/fineTune/test/lora/ --checkpoint_path esm2_t33_650M_UR50D
#python scripts/LitESM2ne_Lora_trainer.py -i data/cellular/metadata/YAP1_HUMAN_Fields2012_singles.csv -o results/fineTune/test/lora/ --checkpoint_path esm2_t33_650M_UR50D


def parse_args():
    parser = argparse.ArgumentParser(description="Training script arguments.")
    parser.add_argument("-i", "--data", type=str, required=True, help="Path to input CSV file.")
    parser.add_argument("-o", "--output", type=str, required=True, help="Output CSV file path.")
    parser.add_argument("--checkpoint_path", type=str, default="esm2_t33_650M_UR50D", help="Model checkpoint name or full path.")
    parser.add_argument("--num_classes", type=int, default=1, help="Number of classes (1 for regression).")
    parser.add_argument("--epochs", type=int, default=50, help="Number of epochs.")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size.")
    parser.add_argument("--device", type=str, default=0, help="Device for training.")
    parser.add_argument("--learning_rate", type=float, default=1e-4, help="Learning rate.")
    parser.add_argument("--weight_decay", type=float, default=1e-2, help="Weight decay.")
    parser.add_argument("--dropout", type=float, default=0.2, help="CLS dropout probability.")
    parser.add_argument("--split_strategy", type=str, default="pool_split", help="Split strategy to use.")
    parser.add_argument("--seed", type=int, default=1234, help="Random seed for reproducibility.")

    # LoRA parameters
    parser.add_argument('--lora_r', type=int, default=8, help='Rank of the low-rank decomposition.')
    parser.add_argument('--lora_alpha', type=int, default=64, help='Scaling factor for LoRA weights. alpha/rank.')
    parser.add_argument('--lora_dropout', type=float, default=0.1, help='Dropout rate for LoRA.')
    parser.add_argument('--lora_modules', nargs='*', type=str, default=["q_proj", "v_proj"], help='Modules to apply LoRA. Default is q_proj, v_proj, as implemented by Microsoft.')
    return parser.parse_args()



class MyDataset(Dataset):
    """This class just loads the data and return a dataset object, returning the sequences and targets.
    Without any tokenization or padding.
    This will be handled later in the collate_fn.
    """
    def __init__(self, data_file):
        self.data = data_file
        self.sequences = [(d['ID'], d['sequence']) for _, d in self.data.iterrows()]
        self.targets = [d['target'] for _, d in self.data.iterrows()]

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
            train_df, val_df = split_data(data, seed=args.seed, train_pct=0.8, val_pct=0.2)
        elif self.split_strategy == 'pool_split':
            train_df, val_df = train_test_split(data, test_size=0.2, random_state=args.seed)
        
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
        _, _, tokens = batch_converter(seqs)
        targets = torch.stack(targets) # Stack the targets to a shape (batch_size, 1)
        return tokens, targets

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
        
        # Define the model
        self.model_loader = Load_from_pretrained(
            checkpoint_path=args.checkpoint_path, num_classes=args.num_classes, dropout=args.dropout,
            lora_r=args.lora_r, lora_alpha=args.lora_alpha, lora_dropout=args.lora_dropout,
            lora_modules=args.lora_modules)
        self.model, self.alphabet = self.model_loader.get_model_details()
        
        # metrics and loss function
        self.loss_fn = torch.nn.MSELoss(reduction='mean')
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
        return AdamW(self.parameters(), lr=self.args.learning_rate, weight_decay=self.args.weight_decay)


def main():
    # Parse the arguments, and collect the model and dataset names
    args = parse_args()
    model_name = args.checkpoint_path
    print(f"Model name: {model_name}")
    dataset_name = args.data.split("/")[-1].split(".csv")[0]
    print(f"Dataset name: {dataset_name}")

    print("Loading model and dataset...")
    model = LitModel(args)
    datamodule = MyDataModule(model.alphabet, args)
    
    ########################## Training setup ##########################
    logger = CSVLogger(
            save_dir=os.path.join(f"{args.output}", f"{model_name}"),
            name=f"{dataset_name}",
            version=args.seed,
            )
    
    trainer = pl.Trainer(
        accelerator="gpu",
        devices=[args.device],
        max_epochs= args.epochs,
        enable_checkpointing=False,
        gradient_clip_val=1.0,  
        logger=logger,
        callbacks=[EarlyStopping(monitor="val_loss", patience=5, mode="min")],
    )
    trainer.fit(model, datamodule)

  
if __name__ == '__main__':
    args = parse_args()
    main()