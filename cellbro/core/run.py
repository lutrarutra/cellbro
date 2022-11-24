import argparse

import App as App

def main(args):
    app = App.App()
    app.run()

if __name__ ==  "__main__":
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    main(args)