import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def read_file(filename, size):
    data = np.loadtxt(filename)
    # import pdb; pdb.set_trace()
    data = data.reshape(size)
    return data


def scatter_plot(filename, shape):

    size = np.prod(shape)

    data = read_file(filename, size)

    X, Y, Z = np.mgrid[-1:1:(shape[0] * 1j), -1:1:(shape[1] * 1j), -1:1:(shape[2] * 1j)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Title and axis labels
    plt.title('Pressure Field')
    ax.set_xlabel('X-axis', fontweight='bold')
    ax.set_ylabel('Y-axis', fontweight='bold')
    ax.set_zlabel('Z-axis', fontweight='bold')

    # Configure colormap.
    # Many options avaiable. Check all options in: 
    # https://matplotlib.org/tutorials/colors/colormaps.html

    my_cmap = plt.get_cmap('Spectral')
    scat = ax.scatter3D(X, Y, Z, c=data, cmap=my_cmap)

    # Add colorbar to the chart
    fig.colorbar(scat, shrink=0.5, aspect=5)

    # Maximize window
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()

    plt.show()


scatter_plot('input.txt', (7, 7, 7))
