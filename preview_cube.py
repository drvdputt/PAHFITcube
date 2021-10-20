"""Simple script for a quick look at the cube merge results"""

import plotting
import argparse
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("merged_cube_file")
    args = ap.parse_args()
    plotting.plot_cube(args.merged_cube_file, args.merged_cube_file)
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
