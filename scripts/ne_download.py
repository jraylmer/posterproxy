"""Download Natural Earth data and extract the files into the appropriate
subdirectories to be understood by the naturalearth.py module. This script
requires the python requests package (pip install requests).

Natural Earth datasets are available at three different 'scale's (1:110m,
1:50m, and 1:10m, in order of coarsest to finest) for different 'feature's
(land, coastline, etc.).

Usage: python ne_download.py [-s scale] [-f feature] [-d directory]

where the optional flags -s and -f are the scales and features as described
above (can be more than one). Currently only "land" and "coastline" are
possible (as other features are not currently processed in the naturalearth.py
module).

The optional flag -d sets the main output directory. If -d is not specified,
the default directory specified in the configuration is used. If the data is
downloaded and extracted into a different directory, the (user) configuration
file must be updated accordingly.

The flag -a may also be specified to download all scales/features.

"""

from argparse import ArgumentParser
import io
from pathlib import Path
import zipfile

import requests

from posterproxy import config


def main():

    prsr = ArgumentParser(usage="Download Natural Earth datasets")
    prsr.add_argument("-s", "--scale", type=str, default=["110m"],
                      nargs="*", choices=["110m", "50m", "10m"])
    prsr.add_argument("-f", "--feature", type=str, default=["land"],
                      nargs="*", choices=["land", "coastline"])
    prsr.add_argument("-a", "--all", action="store_true",
                      help="Download all datasets")
    prsr.add_argument("-d", "--directory", type=str, default="",
                      help="Save to this directory rather than config/default")
    cmd = prsr.parse_args()

    # All data to be saved in subdirectories of this directory:
    if cmd.directory == "":
        root_dir_out = Path(config.ne_data_dir)
    else:
        root_dir_out = Path(cmd.directory).resolve()
        x = input(f"Saving to {str(root_dir_out)}\nContinue? [y/n] ")
        if "n" in x.lower():
            return

    # URL follows a standard format (replacements are scale, scale, feature):
    url = "https://naciscdn.org/naturalearth/{}/physical/ne_{}_{}.zip"

    # Determine which to download from command-line arguments:
    if cmd.all:
        scales   = ["110m", "50m", "10m"]
        features = ["land", "coastline"]
    else:
        scales   = cmd.scale
        features = cmd.feature

    for scale in scales:
        for feature in features:

            # Make output subdirectory:
            dir_out = Path(root_dir_out, f"ne_{scale}_{feature}").resolve()
            dir_out.mkdir(parents=True, exist_ok=True)

            # Access the data at the URL, download, and extract to subdirectory:
            r = requests.get(url.format(scale, scale, feature))
            z = zipfile.ZipFile(io.BytesIO(r.content))
            z.extractall(dir_out)

            print(f"Saved {scale} {feature} data to: {str(dir_out)}")

if __name__ == "__main__":
    main()

