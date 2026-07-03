import numpy as np
import geomloss
import torch

def distance_computation(true_data, pred_data, mode='both'):
    '''Evaluate the Wasserstein and Energy distances between true and reconstructed data at a single time point.'''
    assert true_data.shape[1] == pred_data.shape[1]
    
    # Helper function to convert to DoubleTensor efficiently
    def to_double_tensor(data):
        if isinstance(data, torch.Tensor):
            return data.double() if data.dtype != torch.float64 else data
        return torch.from_numpy(data).double()
    
    true_data = to_double_tensor(true_data)
    pred_data = to_double_tensor(pred_data)
    
    # Only create the loss functions you need
    if mode == 'wasserstein':
        wasserstein_function = geomloss.SamplesLoss("sinkhorn", p=2, blur=0.05, backend="tensorized")
        return wasserstein_function(true_data, pred_data).item()
    elif mode == 'energy':
        energy_function = geomloss.SamplesLoss("energy", p=2, backend="tensorized")
        return energy_function(true_data, pred_data).item()
    else:  # mode == 'both'
        wasserstein_function = geomloss.SamplesLoss("sinkhorn", p=2, blur=0.05, backend="tensorized")
        energy_function = geomloss.SamplesLoss("energy", p=2, backend="tensorized")
        wasserstein_loss = wasserstein_function(true_data, pred_data).item()
        energy_loss = energy_function(true_data, pred_data).item()
        return wasserstein_loss, energy_loss