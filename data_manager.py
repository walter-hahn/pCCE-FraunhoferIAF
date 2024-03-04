import time
import pickle
import os
from constants import max_size

def save_position_data(data):
    timestr = time.strftime("%Y%m%d")
    HE_data = open("positions_140_spins_2D{}.pkl".format(timestr), "wb")

    pickle.dump(data, HE_data)

    HE_data.close()


def save_result_data_old(data, concentration, thickness):
    timestr = time.strftime("%Y%m%d")
    HE_data = open("P1_{}_{}_{}.pkl".format(concentration, thickness, timestr), "wb")

    pickle.dump(data, HE_data)

    HE_data.close()


def save_result_data(data, concentration, thickness, folder_path='pickleFiles'):
    """
    Saves the given data to a pickle file in a specified folder.

    Parameters:
    data: Data to be saved.
    concentration (float): Concentration parameter used in the file name.
    thickness (float): Thickness parameter used in the file name.
    folder_path (str): Path to the folder where the file will be saved. Defaults to 'results'.
    """
    # Ensure the folder exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Create the file name with time stamp
    timestr = time.strftime("%Y%m%d")
    file_name = "P1_{}_{}_{}.pkl".format(concentration, thickness, timestr)

    # Combine the folder path and file name
    file_path = os.path.join(folder_path, file_name)

    # Save the data
    with open(file_path, "wb") as HE_data:
        pickle.dump(data, HE_data)
def save_avg_result_data_old(data):
    timestr = time.strftime("%Y%m%d")
    HE_data = open("averaged_CCE2+{}_2D_{}_2spins.pkl".format(max_size, timestr), "wb")

    pickle.dump(data, HE_data)

    HE_data.close()

def save_avg_result_data(data, folder_path='pickleFiles'):
    """
    Saves the given averaged data to a pickle file in a specified folder.

    Parameters:
    data: Data to be saved.
    folder_path (str): Path to the folder where the file will be saved. Defaults to 'averaged_results'.
    """
    # Ensure the folder exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Create the file name with time stamp
    timestr = time.strftime("%Y%m%d")
    file_name = "averaged_CCE2+{}_2D_{}_2spins.pkl".format(max_size, timestr)

    # Combine the folder path and file name
    file_path = os.path.join(folder_path, file_name)

    # Save the data
    with open(file_path, "wb") as HE_data:
        pickle.dump(data, HE_data)
def load_data(file):
    a_file = open(file, "rb")

    output = pickle.load(a_file)

    return output
