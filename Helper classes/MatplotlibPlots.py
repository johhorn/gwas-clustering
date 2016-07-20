
class MatplotlibPlots(object):
    """
    Class for basic plotting using matplotlib. For custom figures, an instance of the class Figure (as created by
    matplotlib or pyplot functions) can be passed as argument to saveFigure() along with an instance of
    GalaxyRunSpecficFile. After this call, the GalaxyRunSpecificFile have access to the stored image.

    The functions available in this file serve as examples of implementation, and can be used directly.
    They take as argument data related to the figure, as well as an instance of GalaxyRunSpecificFile,
    so that the figure can be saved directly.

    Example usage from a tool:

    dendrogramFile = GalaxyRunSpecificFile(['Image', 'deondrogram.pdf'], galaxyFn)
    MatplotlibPlots.dendrogramClusteringPlot(linkageMatrix, labels.values(), dendrogramFile)
    print dendrogramFile.getLink('Dendrogram here')

    For plots that should be used in articles, pdf format is recommended. It scales very well, as they do not compress
    or blur the figure at different sizes.
    """

    @classmethod
    def saveFigure(cls, fig, fileLocation):
        from matplotlib import pyplot as plt
        fig.savefig(fileLocation.getDiskPath(ensurePath=True))
        plt.close(fig)

    @classmethod
    def dendrogramClusteringPlot(cls, linkageMatrix, labels, fileLocation):
        from matplotlib import pyplot as plt
        import scipy.cluster.hierarchy as sch
        import numpy as np

        sch.set_link_color_palette(['black'])

        fig, axes = plt.subplots()
        fig.subplots_adjust(bottom=0.4)
        plt.title('Hierarchical Clustering Dendrogram')
        plt.ylabel('Distance')

        sch.dendrogram(
            linkageMatrix,
            labels=labels,
            leaf_rotation=270,   # rotates the x axis labels
            color_threshold=np.inf
        )
        plt.plot()
        cls.saveFigure(fig, fileLocation)

    @classmethod
    def seabornHeatmapPlot(cls, data, labels, fileName, max=1, min=0, cmap='Reds'):
        from matplotlib import pyplot as plt
        import seaborn as sns

        a4_dims = (11.7, 8.27)
        heatmap, axes = plt.subplots(figsize=a4_dims)
        sns.heatmap(data, vmin=min, vmax=max, xticklabels=labels, yticklabels=labels,
                    linewidths=0, square=True, cmap=cmap)
        plt.subplots_adjust(left=0.4, bottom=0.4)
        plt.xticks(rotation=270)
        plt.yticks(rotation=0)
        axes.invert_xaxis()
        cls.saveFigure(heatmap, fileName)

    @classmethod
    def pointGraphY(cls, y, fileLocation, xlabel='', ylabel='', xticks=None):
        import matplotlib.pyplot as plt
        import seaborn as sns

        sns.set_style("darkgrid")

        points, axes = plt.subplots(figsize=(8, 7))
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        if xticks:
            plt.xticks(range(0, len(xticks)), xticks, rotation=270)
            plt.subplots_adjust(bottom=0.3)

        plt.plot(y)
        cls.saveFigure(points, fileLocation)

    @classmethod
    def pointGraph(cls, x, y, fileLocation, xlabel='', ylabel=''):
        import matplotlib.pyplot as plt
        import seaborn as sns

        sns.set_style("darkgrid")

        points, axes = plt.subplots()
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)
        plt.plot(x, y)
        cls.saveFigure(points, fileLocation)

    @classmethod
    def multipleLineGraph(cls, xdata, ydataList, ydataLabels, fileLocation, xlabel='', ylabel=''):
        import matplotlib.pyplot as plt
        import seaborn as sns

        lineCount = len(ydataList)
        labelCount = len(ydataLabels)

        sns.set_style("darkgrid")
        sns.set_palette(sns.color_palette("Blues_d", len(ydataLabels)))

        points, axes = plt.subplots()
        axes.set_xlabel(xlabel)
        axes.set_ylabel(ylabel)

        for i in range(0, lineCount):
            if lineCount == labelCount:
                plt.plot(xdata, ydataList[i], label=ydataLabels[i])
            else:
                plt.plot(xdata, ydataList[i])

        if lineCount == labelCount:
            handles, labels = axes.get_legend_handles_labels()
            axes.legend(handles[::-1], labels[::-1], prop={'size': 11}, loc=1)

        cls.saveFigure(points, fileLocation)

    @classmethod
    def histogramRugPlot(cls, data, bins, fileLocation, label=''):
        from matplotlib import pyplot as plt
        import seaborn as sns

        histogram, axes = plt.subplots()
        sns.distplot(data, bins=bins, kde=False, rug=True, axlabel=label)
        cls.saveFigure(histogram, fileLocation)

    @classmethod
    def histogramPlot(cls, data, bins, fileLocation, label=''):
        from matplotlib import pyplot as plt
        import seaborn as sns

        histogram, axes = plt.subplots()
        sns.distplot(data, bins=bins, kde=False, rug=False, axlabel=label)
        cls.saveFigure(histogram, fileLocation)
