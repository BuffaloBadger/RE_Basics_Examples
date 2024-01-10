import pandas as pd
import numpy as np

# filepaths
filepath_to_data = './reb_19_1/Data/'
filepath_to_results = './reb_19_1/Results/'
filepath_to_figures = '../RE_Basics/Graphics/'

# read the data from the .csv file
df = pd.read_csv(filepath_to_data + 'reb_19_1_data.csv')
        # columns: Experiment, T, CA0, t, CA

# get the temperatures used in the experiments
block_temperatures = df['T'].unique() + 273.15

# create a dataframe to store the data file filenames
filenames = pd.DataFrame(columns=['Data Filename'])

# create and save the data blocks
for T_K in block_temperatures:
    # create the block
    block = df[df['T'] == T_K - 273.15]

    print(block)

    # save the block
    T_as_text = format(T_K-273.15,'.0f')
    filename = 'reb_19_1_data_' + T_as_text + '.csv'
    block.to_csv(filepath_to_data + filename,index=False)

    # save the filename
    filenames.loc[len(filenames.index)] = [filename]

# save the data file filenames
filenames.to_csv(filepath_to_data + 'reb_19_1_data_filenames.csv',index=False)