import glob
import os

# Use automatically the directory where the script is launched
path_to_target = os.getcwd()

# Find mutpos files in the current directory
path_to_file_list = glob.glob(os.path.join(path_to_target, "*_final.bam.pileup.mutpos"))

for feed_path in path_to_file_list:
    feed = os.path.basename(feed_path)

    # Output file name
    out = feed.replace("_final.bam.pileup.mutpos", "_mutposreports")

    # Sample name
    sample = feed.replace(".sscs1_final.bam.pileup.mutpos", "")
    sample = sample.replace("_final.bam.pileup.mutpos", "")

    with open(feed_path, "r") as file:
        with open(out, "w") as outfile:
            for i, line in enumerate(file):
                line = line.rstrip("\n")

                if i == 0:
                    # Add proper header name
                    outfile.write(line + "\tsample\n")
                else:
                    # Add sample name to data rows
                    outfile.write(line + f"\t{sample}\n")
