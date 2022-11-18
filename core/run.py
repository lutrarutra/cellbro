import glob, argparse, sys, time, os, multiprocessing

import App

import matplotlib.pyplot as plt

def plot():
    plt.plot([1,2,1])
    plt.show()

def main(args):
    app = App.App()

    app.loop()
    for process_key in list(app.processes.keys()):
        app.processes[process_key].join()


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("-s", "--sc", type=str, help="Path to scRNA-seq data file (.csv/.tsv)", required=True)
    # args = parser.parse_args()

    args = parser.parse_args()

    main(args)