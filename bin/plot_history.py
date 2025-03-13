#!/usr/bin/env python3

import argparse
import seaborn as sns
import pickle

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Plot training history.')
parser.add_argument('--history_files', type=str, help='Space-separated list of history pickle files.')
args = parser.parse_args()

pp = PdfPages("reconstruction_loss.pdf")

for i in sorted(args.history_files.split()):
    with open(i, 'rb') as f:
        history = pickle.load(f)
    temp = history['reconstruction_loss_train'].merge(history['reconstruction_loss_validation'], 
                                                      how='outer', left_index=True, right_index=True)
    plt.clf()
    sns.lineplot(temp).set_title(i)
    pp.savefig(plt.gcf())
pp.close()

