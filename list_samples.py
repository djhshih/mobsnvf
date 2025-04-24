# %%
import glob
import os
import argparse

# %%
## Parsing the command line arguments
parser = argparse.ArgumentParser(
    description="Strip an optional prefix/suffix from bam file names"
)
parser.add_argument(
    "-p", "--prefix",
    default="",
    help="prefix to remove"
)
parser.add_argument(
    "-s", "--suffix",
    default="",
    help="suffix to remove (without the .bam extension)"
)

args = parser.parse_args()

# %%
## Listing the bam files available in the bam directory
sample_paths = glob.glob("bam/*.bam")

## Removing the prefix and suffix from the sample names
sample_names = [
    os.path.basename(p)
      .removesuffix(".bam")
      .removesuffix(args.suffix)
      .removeprefix(args.prefix)
    for p in sample_paths
]

# %%
## Writing the sample names to a file
with open("samples.txt", "w") as sample:
    for name in sample_names:
        sample.write(f"{name}\n")
