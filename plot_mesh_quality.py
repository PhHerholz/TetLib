import os
from PIL import Image
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("folder" )
args = parser.parse_args()

folder = args.folder
plot_folder = os.path.join(folder, 'plots')
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)

runs = []
for name in os.listdir(folder):
    splt = name.split(".")
    if splt[-1] == 'csv':
        runs.append(name.split('metrics.csv')[0])

for i, run in enumerate(runs):
    print("{}/{}".format(i+1, len(runs)))
    path_base = os.path.join(folder, run)

    df = pd.read_csv(path_base + "metrics.csv")

    opmoved=False
    # parse config
    splt = run.split('_')[:-1]
    n_orbitpoints = int(splt[0])
    cellSize = float(splt[1])
    re_ratio = float(splt[2])
    n_flips  = int(splt[3])
    maxpointmove = float(splt[4])
    use_lloyd = bool(int(splt[5]))
    use_perturb = bool(int(splt[6]))
    use_exude = bool(int(splt[7]))

    if len(splt) > 8:
        n_orbitpoints = int(splt[8])
        opmoved=True

    title_str = "N:{} {},  Cellsize {}, RE Ratio {},   {}{}{} {}Flips".format(
        n_orbitpoints,
        "(moved by max {})".format(maxpointmove) if opmoved else "",
        cellSize,
        re_ratio,
        "+lloyd" if
        use_lloyd
        else "",
        "+perturb" if
        use_perturb
        else "",
        "+exude" if
        use_exude
        else "", n_flips)


    rows = len(df.keys())
    cols = 2

    fig = plt.figure(figsize=(10, 10))
    plt.suptitle(title_str)
    for i, (k,v) in enumerate(df.items()):
        plt.subplot(rows, cols, 1+i*cols)
        plt.title(k)
        plt.axis('off')
        plt.imshow(Image.open(path_base + "{}out.png".format(k)))
        plt.subplot(rows, cols, 2+i*cols)
        df[k].hist(bins=200)
    plt.savefig(os.path.join(plot_folder,
                             "{}_metric_comparison.png".format(run)), dpi=400)
