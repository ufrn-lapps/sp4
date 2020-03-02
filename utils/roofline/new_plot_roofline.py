from json import load
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from sys import argv
from os.path import isfile
import numpy as np
from math import log2, log10
import seaborn as sns


def read_input(filename):

    if not isfile(filename):
        print('File %s does not exist.' % filename)

    with open(filename) as json_file:
        data = load(json_file)
        return data


def get_peak_performance(memory_bandwidth, compute_max, oi):
    """
    For a given operational intensity returns the peak performance associated.

    memory_bandwidth: float
        Given the memory bandwidth in 
    """
    peak_performance = 0

    # Convert OI to float
    oi = float(oi)

    peak_performance = min(oi*memory_bandwidth, compute_max)

    return peak_performance


def configure_plot(data):

    # Paramters
    xmin = data.get("xmin", 0.5)
    xmax = data.get("xmax", 64)
    ymin = data.get("ymin", 2)
    ymax = data.get("ymax", 100000)

    sns.set(style="whitegrid", color_codes=True)

    fig = plt.figure(1, figsize=(10.67, 6.6))
    plt.clf()
    # plt.title(data['title'])
    ax = fig.gca()
    ax.set_xscale('log', basex=2)
    ax.set_yscale('log', basey=2)
    ax.set_xlabel('Operational Intensity (FLOPs/Byte)', fontsize=20)
    ax.set_ylabel('Performance (GFLOPs/sec)', fontsize=20)

    # Set grid
    plt.grid(b=True, which='major')

    # Or if you want different settings for the grids:
    # ax.grid(which='minor', alpha=0.2)
    # ax.grid(which='major', color='#004C9901', linestyle='--')

    ax.xaxis.set_minor_formatter(mticker.ScalarFormatter())
    ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
    # print([2**x for x in range(1+int(log2(xmax)))])
    plt.xticks([2**x for x in range(1+int(log2(xmax)))], fontsize=18)
    ax.set_xlim(xmin, xmax)

    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%d'))
    # print([2**x for x in range(int(log2(ymax)) + 1)])
    plt.yticks([2**x for x in range(int(log2(ymax)) + 1)], fontsize=18)
    ax.set_ylim(ymin, ymax)

    # Plot roofs
    x_values = np.arange(0, xmax, 0.1)
    for roof in data['roofs']:
        print(roof['memory'], roof['compute'])
        y_values = [min(roof['compute'], roof['memory']*x)
                    for x in x_values]
        ax.plot(x_values, y_values, c='k', ls='-', lw='2')

        # Plot compute text above peak line
        if not roof['compute_text'].isspace():
            ax.text(x_values[-1], roof['compute'], roof['compute_text'],
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    fontsize=14)

        # Plot memory text above peak line
        if not roof['memory_text'].isspace():
            alpha = np.arctan(np.log10(ax.get_xlim()[1]/ax.get_xlim()[0]) / np.log10(ax.get_ylim()[1]/ax.get_ylim()[0])
                              * fig.get_size_inches()[1]/fig.get_size_inches()[0])
            start_point = xmin*1.02

            ax.text(xmin*1.02, xmin*roof['memory']*1.05,
                    roof['memory_text'],
                    horizontalalignment='left',
                    verticalalignment='bottom',
                    rotation=180*alpha/np.pi,
                    fontsize=20)

        space = 0.5
        # Plot points
        for n, point in enumerate(roof['points']):
            # Draw the point
            ax.plot(point['ai'], point['gflops'],
                    c=point.get("color", "r"),
                    marker=point.get("marker", "o"),
                    linestyle='-',
                    ms=10)

            # Draw vertical line
            plt.stem([point['ai']],
                     [min(point['ai']*roof['memory'], roof['compute'])],
                     # color=point.get("color", "r"),
                     linefmt='%s--' % point.get("color", "r"),
                     markerfmt='g,')

            # Draw point label
            ax.text(point['ai']*1.01, ymin*1.05 + n*space,
                    point['label'],
                    verticalalignment='bottom',
                    color=point.get("color", "r"),
                    rotation=90,
                    fontsize=14)

            # Draw percent label
            if not point.get('time', ' ').isspace():

                if point['label'] == 'SO=12':
                    ax.text(point['ai']*1.01, point['gflops']*1.3,
                            point['time'],
                            color=point.get("color", "r"),
                            fontsize=14)

                    ax.text(point['ai']*1,
                            min(point['ai'] * roof['memory'],
                                roof['compute'])*0.7,
                            "%.2f%%" % (100 *
                                        float(point['gflops'] /
                                              min(point['ai'] * roof['memory'], roof['compute']))),
                            color=point.get("color", "r"),
                            verticalalignment='bottom',
                            fontsize=14)

                else:
                    ax.text(point['ai']*1.01, point['gflops']*1.1,
                            point['time'],
                            color=point.get("color", "r"),
                            fontsize=14)

                    ax.text(point['ai']*0.95,
                            min(point['ai'] * roof['memory'],
                                roof['compute'])*1.2,
                            "%.2f%%" % (100 *
                                        float(point['gflops'] /
                                              min(point['ai'] * roof['memory'], roof['compute']))),
                            color=point.get("color", "r"),
                            verticalalignment='bottom',
                            fontsize=14)


if __name__ == "__main__":

    # For V100 Single Precision
    # memory = 900
    # compute = 14000
    print(get_peak_performance(900.0, 14000.0, argv[1]))

    # data = read_input(argv[1])
    # configure_plot(data)

    # plt.savefig(filename+'.png')
    # plt.savefig(filename+'.eps')
    # plt.show()
    # plt.savefig(data['title'] + '.eps', format='eps', dpi=1000)
