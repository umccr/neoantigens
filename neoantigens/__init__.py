from os.path import isfile, join, dirname, abspath


def package_path():
    return dirname(abspath(__file__))