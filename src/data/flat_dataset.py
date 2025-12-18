import torch
from torch.utils.data import Dataset
import numpy as np
import re

class FlattenedDataset(Dataset):
    def __init__(self, gtex, label_position=-1, scaling_type='mean'):
        self.scaling_type = scaling_type
        self.label_position = label_position

        # Extract sample labels (tumor type)
        self.labels = gtex.iloc[:, label_position].values

        # Extract mutation context names
        self.kmer_names = gtex.drop(gtex.columns[[label_position]], axis=1).columns.tolist()

        # Extract mutation labels  
        self.mut_labels = []
        for kmer in self.kmer_names:
            # Extract mutation label from format like [C>T]
            match = re.search(r'\[([ACGT]>[ACGT])\]', kmer.upper())
            mut_label = f"[{match.group(1)}]" if match else "[UNK]"
            self.mut_labels.append(mut_label)
        
        # Extract input features only (excluding label column)
        self.mut_matrix = torch.tensor(
            gtex.drop(gtex.columns[[label_position]], axis=1).values
        ).float()

        self.num_samples, self.num_muts = self.mut_matrix.shape
        # Precompute lib sizes
        if self.scaling_type == 'mean':
            self.lib_values = self.mut_matrix.mean(dim=1, keepdim=True)
        elif self.scaling_type == 'max':
            self.lib_values = self.mut_matrix.max(dim=1, keepdim=True).values
        elif self.scaling_type == 'sum':
            self.lib_values = self.mut_matrix.sum(dim=1, keepdim=True)
        else:
            raise ValueError("Invalid scaling_type")

        # Precompute the flattened index map: one row per sample–mutation pair
        self.flattened_indices = [
            (sample_idx, mut_idx)
            for sample_idx in range(self.num_samples)
            for mut_idx in range(self.num_muts)
        ]

        # onehot encode each column name
        
        # One-hot encoder for nucleotides
        nuc_to_vec = {'A': [1, 0, 0, 0],
                      'C': [0, 1, 0, 0],
                      'G': [0, 0, 1, 0],
                      'T': [0, 0, 0, 1]}
        
        # Convert k-mer string to one-hot matrix
        def onehot_kmer(kmer):
            nucleotides = re.findall(r'[ACGT]', kmer.upper())
            return np.array([nuc_to_vec[n] for n in nucleotides], dtype=np.float32) 
        
        # Dictionary of one-hot k-mers: shape (k, 4) for each
        self.kmer_onehots = {
            kmer: torch.tensor(onehot_kmer(kmer)) for kmer in self.kmer_names
        }

    def __len__(self):
        # Total number of sample–mutation pairs
        return len(self.flattened_indices)

    def __getitem__(self, idx):
        sample_idx, mut_idx = self.flattened_indices[idx] # idx are the indices for the current batch
        value = self.mut_matrix[sample_idx, mut_idx].unsqueeze(0)

        # Use precomputed lib
        lib = self.lib_values[sample_idx]
        
        kmer = self.kmer_names[mut_idx]
        onehot = self.kmer_onehots[kmer].flatten() # onehot encoding for the given mutation context (mut_idx)

        return value, lib, sample_idx, mut_idx, onehot # add return of the onehot encoded mut context