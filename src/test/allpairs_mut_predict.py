import torch
from tqdm import tqdm
from src.dgd.latent import RepresentationLayer
from itertools import product # added

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def prepare_potential_reps(sample_list):
    """
    takes a list of samples drawn from the DGD's distributions.
    The length gives the number of distributions which defines
    the dimensionality of the output tensor.
    If the list of samples is longer than 1, we will create representations
    from the combination of each GMM's samples.
    """
    return sample_list[0]
    
def learn_new_representation(dgd, 
                             data_loader,
                             test_epochs=50,
                             learning_rates=1e-2, 
                             weight_decay=0.,
                             betas=(0.5, 0.7),
                             reduction_type="sum",
                             resampling_type="mean"):
    """
    This function learns a new representation layer for the DGD.
    The new representation layer is learned by sampling new points
    from the GMMs and finding the best fitting GMM for each sample.
    The new representation layer is then optimized to minimize the
    reconstruction loss of the DGD.
    """
    def check_devices(**tensors):
        for name, t in tensors.items():
            print(f"{name} → {t.device}")

    gmm_loss = True
    mut_gmm_loss = True # added
    
    n_samples_new = len(data_loader.dataset)
    potential_reps = prepare_potential_reps([dgd.gmm.sample_new_points(resampling_type)]).to(device) # initialize reps. The function sample_new_points samples as many points as there are Gaussians (their means)
    mut_potential_reps = prepare_potential_reps([dgd.mut_gmm.sample_new_points(resampling_type)]).to(device) # added
    #print(potential_reps)
    #print(mut_potential_reps)
    
    # Get the number of combinations
    n_samples = potential_reps.shape[0]
    n_mut = mut_potential_reps.shape[0]
    
    # Repeat and tile the tensors to align all pairs
    sample_rep_exp = potential_reps.unsqueeze(1).repeat(1, n_mut, 1)      # shape: [n_samples, n_mut, sample_dim]
    mut_rep_exp = mut_potential_reps.unsqueeze(0).repeat(n_samples, 1, 1) # shape: [n_samples, n_mut, mut_dim]
    
    # Concatenate along the feature dimension
    all_pairs = torch.cat([sample_rep_exp, mut_rep_exp], dim=-1)          # shape: [n_samples, n_mut, sample_dim + mut_dim]
    
    # Flatten to shape [n_samples * n_mut, total_dim]
    all_pairs_flat = all_pairs.view(-1, sample_rep_exp.shape[-1] + mut_rep_exp.shape[-1])

    
    dgd.eval() # evaluation mode
    X_test = dgd.decoder(all_pairs_flat.to(device)) # reconstructed data. I have to add the onehot vector for each mut type here

    rep_init_values = torch.zeros((n_samples_new, all_pairs_flat.shape[-1]), device=device) # added device
    print("rep_init_values device:", rep_init_values.device) # added
    # check_devices(rep_init_values=rep_init_values, potential_reps=potential_reps)

    # this first loop is for initialization of the representations
    for mut_data, lib, sample_idx, mut_idx in data_loader: # loop through samples
        loss = torch.empty(0).to(device) # placed on the same hardware for computations
        for X in X_test:   # loop through each reconstruction
            mut_recon_loss = dgd.decoder.loss( # 
                nn_output=X.to(device), 
                target=mut_data.to(device), 
                scale=lib,
                feature_ids=None, 
                reduction="sum"
            )
            loss = torch.cat((loss, mut_recon_loss.unsqueeze(0))) # add the reconstruction loss of mutation X to the accumulated loss. mut_gmm_loss?
        best_fit_ids = torch.argmin(loss, dim=-1).detach() # removed .cpu()
        rep_init_values[i, :] = all_pairs_flat.clone()[best_fit_ids, :] # re-cloning not necessary
        
    # Split into separate reps again to match the DGD class functions in the notebooks
    Ntest=len(data_loader.dataset)
    new_rep = RepresentationLayer(n_rep=dgd.rep_dim, # set-up the representation layer with the best values found above
                                  n_sample=Ntest,
                                  value_init=rep_init_values[:,:dgd.rep_dim]).to(device)
    
    mut_new_rep = RepresentationLayer(n_rep=dgd.mut_rep_dim, # set-up the representation layer with the best values found above
                                  n_sample=Ntest,
                                  value_init=rep_init_values[:,dgd.rep_dim:]).to(device)

    
    test_rep_optimizer = torch.optim.AdamW(new_rep.parameters(), lr=learning_rates, weight_decay=weight_decay, betas=betas) # set-up optimizer
    mut_test_rep_optimizer = torch.optim.AdamW(mut_new_rep.parameters(), lr=learning_rates, weight_decay=weight_decay, betas=betas) # added
    
    for epoch in tqdm(range(test_epochs)): # loop through epochs for a progress bar #tqdm()
        test_rep_optimizer.zero_grad()
        mut_test_rep_optimizer.zero_grad()
        for mut_data, lib, sample_idx, mut_idx in data_loader: # loop through samples
            mut_recon_loss, gmm_loss, mut_gmm_loss = dgd.forward_and_loss(
                z=new_rep(index),
                mut_z=mut_new_rep(index),
                target=[mut_data.to(device)],
                scale=lib.unsqueeze(1).to(device), # dimensional alignment with unsqueeze
                gmm_loss=gmm_loss,
                mut_gmm_loss=mut_gmm_loss,
                reduction=reduction_type
            )
            loss = mut_recon_loss + gmm_loss + mut_gmm_loss # added
            loss.backward()
        test_rep_optimizer.step() # updating the rep for each epoch
        mut_test_rep_optimizer.step() # added
    
    return new_rep, mut_new_rep


    

    
    # Cartesian product
    #all_pairs = torch.cartesian_prod(potential_reps, mut_potential_reps)
    #all_pairs_flat = all_pairs.view(-1, sample_dim + mut_dim)

    # Combine into single tensor
    #concat_potential_reps = torch.cat([all_pairs[:, 0, :], all_pairs[:, 1, :]], dim=1)

    #concat_potential_reps = concat_potential_reps = torch.cat([potential_reps, mut_potential_reps], dim=1) # concatenated, added
    #print("concat_potential_reps device:", concat_potential_reps.device) # added