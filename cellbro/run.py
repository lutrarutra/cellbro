from cellbro.core.app import App
import argparse
# from cellbro.util.Dataset import Dataset

# dataset = Dataset("/home/lutrarutra/Documents/dev/bioinfo/cellbrowser/data/vas.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("path", help="Path to the dataset")
    parser.add_argument("--debug", help="Run in debug mode", action="store_true")
    parser.add_argument("--host", help="Host address", default="127.0.0.1")
    parser.add_argument("--port", help="Port", default=8050)

    args = parser.parse_args()

    app = App(args.path)
    app.run(debug=args.debug, host=args.host, port=args.port)
