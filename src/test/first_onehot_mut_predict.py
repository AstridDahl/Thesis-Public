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
    # add seed
    # HERE
    
    gmm_loss = True
    mut_gmm_loss = True # added
    
    potential_reps = prepare_potential_reps([dgd.gmm.sample_new_points(resampling_type)]).to(device) # initialize reps.
    mut_potential_reps = prepare_potential_reps([dgd.mut_gmm.sample_new_points(resampling_type)]).to(device) # added. The function sample_new_points samples as many points as there are Gaussians        (their means)
    #print(potential_reps)
    #print(mut_potential_reps)
    
    # number of reps
    Nsample=data_loader.dataset.num_samples
    Nmut=data_loader.dataset.num_muts
    #print(Nsample)
    #print(Nmut)
    
    random_indices = torch.randint(0, potential_reps.shape[0], (Nsample,))
    initial_reps = potential_reps[random_indices].clone().detach()

    mut_random_indices = torch.randint(0, mut_potential_reps.shape[0], (Nmut,))
    mut_initial_reps = mut_potential_reps[mut_random_indices].clone().detach()

    dgd.eval() # evaluation mode

    new_rep = RepresentationLayer(n_rep=dgd.rep_dim, # set-up the representation layer at origin
                                  n_sample=Nsample,
                                  value_init=initial_reps.to(device))
    
    mut_new_rep = RepresentationLayer(n_rep=dgd.mut_rep_dim, # set-up the representation layer at origin
                                  n_sample=Nmut,
                                  value_init=mut_initial_reps.to(device))

    
    test_rep_optimizer = torch.optim.AdamW(new_rep.parameters(), lr=learning_rates, weight_decay=weight_decay, betas=betas) # set-up optimizer
    mut_test_rep_optimizer = torch.optim.AdamW(mut_new_rep.parameters(), lr=learning_rates, weight_decay=weight_decay, betas=betas) # added
    
    for epoch in tqdm(range(test_epochs)): # loop through epochs for a progress bar #tqdm()
        test_rep_optimizer.zero_grad()
        mut_test_rep_optimizer.zero_grad()
        for mut_data, lib, sample_idx, mut_idx, onehot in data_loader: # loop through samples
            mut_recon_loss, gmm_loss, mut_gmm_loss = dgd.forward_and_loss(
                z=new_rep(sample_idx).to(device),
                mut_z=mut_new_rep(mut_idx).to(device),
                onehot=onehot.to(device),
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