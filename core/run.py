import glob, argparse, sys, time, os, multiprocessing

import App

import matplotlib.pyplot as plt

def main(args):
    app = App.App()
    app.loop()

    for process_key in list(app.processes.keys()):
        app.processes[process_key].join()


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    main(args)