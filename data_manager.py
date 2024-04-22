import time
import pickle
import os


def save_position_data(data):
    timestr = time.strftime("%Y%m%d")
    HE_data = open("positions_140_spins_2D{}.pkl".format(timestr), "wb")

    pickle.dump(data, HE_data)

    HE_data.close()


def save_result_data(data, concentration, thickness, folder_path='pickleFiles', mode="regular"):
    """
    Saves the given data to a pickle file in a specified folder.

    Parameters:
    data: Data to be saved.
    concentration (float): Concentration parameter used in the file name.
    thickness (float): Thickness parameter used in the file name.
    folder_path (str): Path to the folder where the file will be saved. Defaults to 'pickleFiles'.
    mode (str): Specifies whether it's for regular data or averaged data. Defaults to "regular".
    """
    # Ensure the folder exists
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)

    # Create the file name with time stamp
    timestr = time.strftime("%Y%m%d")
    if mode == "regular":
        file_name = "all_results_{}_{}_{}.pkl".format(concentration, thickness, timestr)
    elif mode == "avg":
        file_name = "averaged_{}_{}_{}.pkl".format(concentration, thickness, timestr)
    else:
        raise ValueError("Invalid mode. Use 'regular' or 'avg'.")

    # Combine the folder path and file name
    file_path = os.path.join(folder_path, file_name)

    # Save the data
    with open(file_path, "wb") as HE_data:
        pickle.dump(data, HE_data)


def load_data(file):
    a_file = open(file, "rb")

    output = pickle.load(a_file)

    return output
