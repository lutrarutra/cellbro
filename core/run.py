import glob, argparse, sys, time, os

import App

def main(args):
    app = App.App()
    app.loop()


if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("-s", "--sc", type=str, help="Path to scRNA-seq data file (.csv/.tsv)", required=True)
    # args = parser.parse_args()

    args = parser.parse_args()

    main(args)