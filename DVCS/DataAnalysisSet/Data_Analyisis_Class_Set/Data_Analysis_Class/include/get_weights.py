import numpy
import pandas
from hep_ml import reweight
from matplotlib import pyplot as plt
import uproot
import pandas as pd
import numpy as np
import awkward as ak
from sklearn.model_selection import train_test_split
from hep_ml.metrics_utils import ks_2samp_weighted
import sys

hist_settings = {'bins': 100, 'density': True, 'alpha': 0.7}

def draw_distributions(original, target, new_original_weights):
    plt.figure(figsize=[15, 7])
    for id, column in enumerate(branches, 1):
        xlim = numpy.percentile(numpy.hstack([target[column]]), [0.01, 99.99])
        plt.subplot(2, 3, id)
        plt.hist(original[column], weights=new_original_weights, range=xlim, **hist_settings)
        plt.hist(target[column], range=xlim, **hist_settings)
        plt.title(column)
        print('KS over ', column, ' = ', ks_2samp_weighted(original[column], target[column], 
                                                           weights1=new_original_weights, weights2=numpy.ones(len(target), dtype=float)))
        
def draw_distributions2(original, target, new_original_weights, branches):
    plt.figure(figsize=[7.5, 3.5])
    for id, column in enumerate(branches, 1):
        xlim = numpy.percentile(numpy.hstack([target[column]]), [0.01, 99.99])
        plt.hist(original[column], weights=new_original_weights, range=xlim, **hist_settings,label="MC")
        plt.hist(target[column], range=xlim, **hist_settings,label="Data")
        plt.title(column)
        plt.legend(loc="upper left")
        print('KS over ', column, ' = ', ks_2samp_weighted(original[column], target[column], 
                                         weights1=new_original_weights, weights2=numpy.ones(len(target), dtype=float)))

def main(args):
    global branches
    branches = ['strip_El_P','strip_El_Theta','strip_Ph1_P','strip_Ph1_Theta','strip_Ph2_P','strip_Ph2_Theta','mm2_egg']
    global branches0
    branches0 = ['strip_El_P','strip_El_Theta','strip_Ph1_P','strip_Ph1_Theta','strip_Ph2_P','strip_Ph2_Theta']
    MC_File = args[0]
    Data_File = args[1]
    #MC_File = 'Tested_BM_Sim.root'
    #Data_File = 'Data_NP_Theta_g_5.root'

    ####################################
    #Get info from MC
    ####################################
    print("Getting info from MC...")
    
    file = uproot.open(MC_File)
    original = file['eppi0']  
    data_dict = {branch: [] for branch in branches}
    
    for arrays in original.iterate(branches):
        for branch in branches:
            # Extract the specific component from each vector in the branch
            data_dict[branch].extend(arrays[branch].tolist())
            
    indices = []
    for arrays in original.iterate('bestCandidateFlag'):
        vector = arrays['bestCandidateFlag']
        for i, vec in enumerate(vector):
            # Find indices where the value is 1
            indices.extend(np.where(vec == 1)[0])

    #print(len(indices),len(data_dict['strip_El_P']))

    df = pd.DataFrame({branch: [entry[indices[count]] for count, entry in enumerate(data_dict[branch])] for branch in branches})
    original0=df
    original=original0.drop(['mm2_egg'],axis=1)


    ####################################
    #Get info from Data
    ####################################
    print("Getting info from Data...")

    file = uproot.open(Data_File)
    target = file['eppi0']  
    data_dict = {branch: [] for branch in branches}
    
    for arrays in target.iterate(branches):
        for branch in branches:
            # Extract the specific component from each vector in the branch
            data_dict[branch].extend(arrays[branch].tolist())

    indices = []
    for arrays in target.iterate('bestCandidateFlag'):
        vector = arrays['bestCandidateFlag']
        for i, vec in enumerate(vector):
            # Find indices where the value is 1
            indices.extend(np.where(vec == 1)[0])

    #print(len(indices),len(data_dict['strip_El_P']))

    df = pd.DataFrame({branch: [entry[indices[count]] for count, entry in enumerate(data_dict[branch])] for branch in branches})
    target0=df
    target=df.drop(['mm2_egg'],axis=1)
    
    ####################################
    #Prepare train and test samples
    ####################################
    print("Training...")
    
    # divide original samples into training ant test parts
    original_train, original_test = train_test_split(original)
    # divide target samples into training ant test parts
    target_train, target_test = train_test_split(target)
    
    #Remove mm2_egg for the training
    branches = branches0
    original_weights = numpy.ones(len(original))
    original_weights_train = numpy.ones(len(original_train))
    original_weights_test = numpy.ones(len(original_test))
    
    ####################################
    #Plot some distribution for control
    ####################################

    #Original distributions
    draw_distributions(original, target, original_weights)
    plt.savefig('Reweighting_Plots/dists_orig.pdf')
    #Train sample distributions
    draw_distributions(original_train, target_train, original_weights_train)
    plt.savefig('Reweighting_Plots/dists_train.pdf')
    #Test sample distributions
    draw_distributions(original_test, target_test, original_weights_test)
    plt.savefig('Reweighting_Plots/dists_test.pdf')
    
    ####################################
    #Gradient Boosted Reweighter
    ####################################

    reweighter = reweight.GBReweighter(n_estimators=50, learning_rate=0.25, max_depth=5, min_samples_leaf=1000, 
                                       gb_args={'subsample': 0.4})
    reweighter.fit(original_train, target_train)    
    gb_weights_test = reweighter.predict_weights(original_test)
    # validate reweighting rule on the test part comparing 1d projections
    draw_distributions(original_test, target_test, gb_weights_test)
    plt.savefig('Reweighting_Plots/dists_reweighted.pdf')

    #Plot mm2_egg
    branches3 = ['mm2_egg']
    weights_all = numpy.ones(len(original))
    draw_distributions2(original0, target0, weights_all, branches3)
    plt.savefig('Reweighting_Plots/mm2_egg_before.pdf')
    weights_all = reweighter.predict_weights(original)
    draw_distributions2(original0, target0, weights_all, branches3)
    plt.savefig('Reweighting_Plots/mm2_egg_reweighted.pdf')
    
    ####################################
    #Folding reweighter:
    # With `FoldingReweighter` one can simpler
    # do cross-validation and in the end obtain
    # unbiased weights for the whole original
    # sample
    ####################################
    print("Folding reweighting...")
    
    # define base reweighter
    reweighter_base = reweight.GBReweighter(n_estimators=50, 
                                            learning_rate=0.25, max_depth=5, min_samples_leaf=100, 
                                            gb_args={'subsample': 0.4})
    reweighter = reweight.FoldingReweighter(reweighter_base, n_folds=2)
    # it is not needed divide data into train/test parts; rewighter can be train on the whole samples
    reweighter.fit(original, target)
    
    # predict method provides unbiased weights prediction for the whole sample
    # folding reweighter contains two reweighters, each is trained on one half of samples
    # during predictions each reweighter predicts another half of samples not used in training
    folding_weights = reweighter.predict_weights(original)
    
    draw_distributions(original, target, folding_weights)
    
    #Plot weighted mm2_egg
    weights_all = reweighter.predict_weights(original)
    draw_distributions2(original0, target0, weights_all, branches3) 
    plt.savefig('Reweighting_Plots/mm2_egg_Folding_reweighted.pdf')
    
    ####################################
    #Create Output
    ####################################
    print("Creating output...")

    # Define the output text file path
    output_text_file = 'rWeights.dat'
    
    # Write the list to the text file
    with open(output_text_file, 'w') as file:
        for item in weights_all:
            file.write(f"{item}\n")

    print(f"List written to: {output_text_file}")
    print("Done.")

if __name__ == "__main__":
    args = [str(x) for x in sys.argv[1:]]     
    x =main(args)
    print(x)
