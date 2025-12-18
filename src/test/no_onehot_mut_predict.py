import torch
from tqdm import tqdm
from src.dgd.latent import RepresentationLayer

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
                             train_loader,
                             test_loader,
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

    gmm_loss = True
    mut_gmm_loss = True
    
    Nsample=test_loader.dataset.num_samples
    
    potential_reps = dgd.gmm.sample_new_points(resampling_type).to(device) # initialize reps. One mean value per component. Remove the list wrapper and prepare()

    dgd.eval() # evaluation mode


    rep_init_values = torch.zeros((Nsample, potential_reps.shape[-1]), device=device) # added device
    #print(rep_init_values.shape)

    #Nmut=train_loader.dataset.num_muts
    mut_rep = dgd.mut_train_rep.z.to(device) # use the mut reps found in the training
    #print(mut_rep)
    # this first loop is for initialization of the sample representations
    for mut_data, lib, sample_idx, mut_idx in test_loader: # loop through batches
        loss = torch.empty(0).to(device) # placed on the same hardware for computations
        for rep in potential_reps:  # loop through each potential initial representation, determined by the number of GMM components             hosen
            rep_batch = rep.unsqueeze(0).expand(len(mut_idx), -1).to(device) # giving the potential rep the same batch dimension as the             rest of the data
            X = dgd.forward(rep_batch,mut_rep[mut_idx].to(device)) # forward instead of decoder, as it combines the three input types                before the forward pass for the final output.                The loss function in decoder.py
            #print(X) # a tensor wrapped in a list
            mut_recon_loss = dgd.decoder.loss( 
                nn_output=X[0].to(device), # unwrap bc decoder outputs a list of output modules, and I only use one modality
                target=mut_data.to(device), 
                scale=lib.to(device),
                feature_ids=None, 
                reduction="mean" # sum or mean? both makes sense
            )
            loss = torch.cat((loss, mut_recon_loss.unsqueeze(0))) # add the reconstruction loss of representation X to the 1D tensor of               losses
        best_fit_ids = torch.argmin(loss, dim=-1).detach() # removed .cpu()
        rep_init_values[sample_idx[0], :] = potential_reps.clone()[best_fit_ids, :] # take the sample index once, not the whole list of          the same index, bc we only want one rep per sample           not one per row
        #print(rep_init_values) # looks right, 374 prints (batches)
    
    #print("Expected shape:", (Nsample, dgd.rep_dim))
    #print("value_init shape:", rep_init_values.shape)
    new_rep = RepresentationLayer(n_rep=dgd.rep_dim, # set-up the representation layer with the best values found above
                                  n_sample=Nsample,
                                  value_init=rep_init_values).to(device)
    
    test_rep_optimizer = torch.optim.AdamW(new_rep.parameters(), lr=learning_rates, weight_decay=weight_decay, betas=betas) # set-up optimizer
    
    for epoch in tqdm(range(test_epochs)): # loop through epochs for a progress bar #tqdm()
        test_rep_optimizer.zero_grad()
        for mut_data, lib, sample_idx, mut_idx in test_loader: # loop through samples
            #print(new_rep(sample_idx))
            #print(mut_rep[mut_idx])
            mut_recon_loss, gmm_loss, mut_gmm_loss = dgd.forward_and_loss(
                z=new_rep(sample_idx).to(device),
                mut_z=mut_rep[mut_idx].to(device), # standard indexing bc it's a torch.nn.Parameter
                target=[mut_data.to(device)],
                scale=lib.unsqueeze(1).to(device), # dimensional alignment with unsqueeze
                gmm_loss=gmm_loss,
                mut_gmm_loss=mut_gmm_loss,
                reduction=reduction_type
            )
            loss = mut_recon_loss + gmm_loss + mut_gmm_loss
            loss.backward()
        test_rep_optimizer.step() # updating the rep for each epoch
    
    return new_rep 