import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

def Save_Map(outfile, arr, labels_x, labels_y, t_hold):
    plt.figure(figsize=(10, 8))

    labels_x1 = [i[1] for i in labels_x]
    labels_y1 = [i[1] for i in labels_y]

    df_cm = pd.DataFrame(arr, index = [i for i in labels_y], columns = [i for i in labels_x])
    plt.title('Contact Map\n\nThreshold = ' + str(t_hold) + '\n')

    myColors = ((0.0, 0.0, 0.0, 0.0), (0.0, 0.8, 0.0, 1.0), (0.8, 0.0, 0.0, 1.0), (0.0, 0.0, 0.8, 1.0), (1.0, 1.0, 0.0, 1.0), (0.5, 0.5, 0.5, 1.0), (0.0, 0.0, 0.0, 1.0), (0.63, 0.13, 0.94, 1.0))
    cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

    ax1 = sns.heatmap(df_cm, cmap=cmap, linewidths=.5, linecolor='lightgray', square=True, annot=False, cbar_kws={'label': ""},vmax=7)

    # Manually specify colorbar labelling after it's been generated
    colorbar = ax1.collections[0].colorbar
    colorbar.set_ticks([0, 1, 2, 3, 4, 5, 6, 7])
    colorbar.set_ticklabels(["no contact", "Phosphate", "Ribose", "Watson-Crick", "Sugar", "Hoogstean", "Face", "Multiclass"])

    plt.xlabel('\nResidue Protein')
    plt.ylabel('\nResidue NA')

    plt.yticks(np.arange(0.5,len(labels_y),1), labels=labels_y1, rotation=0)
    plt.xticks(np.arange(0.5,len(labels_x),1), labels=labels_x1, rotation=0)

    plt.savefig(outfile + "_map_fig_" + str(t_hold) + ".png", dpi=300)