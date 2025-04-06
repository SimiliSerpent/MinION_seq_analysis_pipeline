import time
import sys
import os
import json
import pathlib


def check_dir(path_to_dir):
    """Check that specified directory actually exists.

    Positionnal arguments:
    path_to_dir (str) - path to the specified directory to check
    """
    if not os.path.isdir(path_to_dir):
        raise ValueError('The specified directory does not exist:', path_to_dir)

    return


def check_file(path_to_file):
    """Check that specified file actually exists.

    Positionnal arguments:
    path_to_file (str) - path to the specified file to check
    """
    if not os.path.isfile(path_to_file):
        raise ValueError('The specified file does not exist:', path_to_file)

    return


def print_json(json_object):
    """Pretty print json object.
    """
    print(json.dumps(json_object, indent=4))

    return 0


def save_json(path_to_file, json_object):
    """Save dictionnary to desired location.

    Positionnal arguments:
    path_to_file (str) - path to the storage location of the json object
    """
    with open(path_to_file, 'w') as json_file:
        json_file.write(json.dumps(json_object))

    return


def load_json(path_to_file):
    """Load dictionnary from specified location.

    Positionnal arguments:
    path_to_file (str) - path to the storage location of the json object

    Return:
    json_object (dict) - dictionnary / json loaded from the file
    """
    check_file(path_to_file)
    with open(path_to_file) as json_file:
        json_object = json.loads(json_file.read())

    return json_object


def send_text(text, verbosity, verbosity_threshold=0, indentation=0):
    """Print message to standard output.

    Positionnal arguments:
    text (str) - the message to send to std output
    verbosity (int) - the requested level of verbosity

    Keywords arguments:
    verbosity_threshold (int) - the minimum level of verbosity
                                required to display the message (default 0)
    indentation (int) - the desired indentation level (default 0)
    """
    if verbosity >= verbosity_threshold:
        print(time.strftime("%Hh%M\'%Ss ", time.localtime())
              + indentation*'\t'
              + '- '
              + text)

    return


def get_unblocked(path_to_unblocked_file, v=0):
    """Retrieve the ids of unblocked reads.

    Positionnal arguments:
    path_to_unblocked_file (str) - the path to the file containing the list of
                                   unblocked read ids

    Keywords arguments:
    v (int) - the requested level of verbosity (default 0)

    Return:
    unblocked_ids (list of str) - list of ids of unblocked reads
    """
    # read all lines in the txt file
    send_text('Loading ids', v, 2, 1)
    with open(path_to_unblocked_file) as unblock_file:
        all_lines = unblock_file.readlines()

    # retrieve the ids
    send_text('Stripping lines', v, 2, 1)
    unblocked_ids = [line.strip() for line in all_lines]

    return unblocked_ids


def get_chan_ids(path_to_chan_toml, v=0):
    """Retrieve the ids of channels and whether they perform AS or not.

    Positionnal arguments:
    path_to_chan_toml (str) - the path to the toml file with the ids of the
                              channels for each region

    Keywords arguments:
    v (int) - the requested level of verbosity (default 0)

    Return:
    region_ids (dict) - dictionnary containing, for each region, the list of
                        channels ids
    """
    # initialize regions ids dictionnary
    chan_ids = {'AS': [], 'control': []}
    # read all lines in the txt file
    send_text('Loading regions', v, 2, 1)
    with open(path_to_chan_toml) as regions_file:
        all_lines = regions_file.readlines()

    # retrieve the ids
    send_text('Parsing regions', v, 2, 1)
    get_ids = False
    for line in all_lines:
        if line[:4] == 'name':
            region_name = line.split('\"')[1]
            if 'control' in region_name:
                region_tag = 'control'
            else:
                region_tag = 'AS'
            get_ids = True
        elif get_ids:
            chan_ids[region_tag] += line.strip()[12:-1].split(', ')
            get_ids = False

    return chan_ids


def find_file(path_to_dir, file_subname):
    """Retrieve the names of files containing the given string in the given dir.

    Retrieve the name of all the files in the specified directory
    (recursively) and having the specified sub-string in their file name.

    Positionnal arguments:
    path_to_dir (str) - the path to the directory where to look for the files
    file_subname (str) - the file name sub-string to look for

    Return:
    file_names (list of str) - list of the names of the found files
    """
    found_paths = [
        r_path for r_path in pathlib.Path(path_to_dir).rglob('*' + file_subname + '*')
    ]
    file_names = [
        r_path.name for r_path in found_paths if os.path.isfile(r_path)
    ]
    return file_names

def find_path(path_to_dir, file_subname):
    """Retrieve the paths to files containing the given string in the given dir.

    Retrieve the relative paths to all the files in the specified directory
    (recursively) and having the specified sub-string in their file name.

    Positionnal arguments:
    path_to_dir (str) - the path to the directory where to look for the files
    file_subname (str) - the file name sub-string to look for

    Return:
    file_paths (list of str) - list of the relative Posix paths to the found
        files
    """
    found_paths = [
        r_path for r_path in pathlib.Path(path_to_dir).rglob('*' + file_subname + '*')
    ]
    file_paths = [
        r_path.name for r_path in found_paths if os.path.isfile(r_path)
    ]
    return file_paths


def get_channels_number(path_to_icarust_ini):
    """Retrieve the total number of channels.

    Positionnal arguments:
    path_to_icarust_ini (str) - the path to the Icarust toml .ini config file

    Return:
    channel_number (int) - the total number of channel for the simulation
    """
    with open(path_to_icarust_ini) as ini_file:
        ini_lines = ini_file.readlines()
    channel_number = ''
    for ini_line in ini_lines:
        if ini_line.startswith('channels'):
            channel_number = int(ini_line.split('=')[1].strip())
            break
    if channel_number == '':
        raise ValueError('Failed to find channels number in',
                         path_to_icarust_ini)
    return channel_number


if __name__ == "__main__":
    raise ValueError('You are not supposed to use this as main!')
    sys.exit(1)
