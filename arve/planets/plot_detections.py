from   matplotlib.colors  import LinearSegmentedColormap
import matplotlib.pyplot  as     plt
from   matplotlib.ticker  import AutoMinorLocator
import numpy              as     np

class plot_detections:

    def plot_detections(self, vmin:float=0, vmax:float=2) -> plt.Figure:
        """Plot detections.

        :param vmin: minimum value in color map, defaults to 0
        :type vmin: float, optional
        :param vmax: maximum value in color map, defaults to 2
        :type vmax: float, optional
        :return: figure with injection-recovery test results
        :rtype: plt.Figure
        """

        # read detections
        detections = self.detections
        x_var, y_var, x_arr, y_arr, x_map, y_map = [detections[key] for key in ["x_var", "y_var", "x_arr", "y_arr", "x_map", "y_map"]]
        detection_arr, detection_map, scale      = [detections[key] for key in ["detection_arr", "detection_map", "scale"]]

        # include array and map
        include_arr = detection_arr is not None
        include_map = detection_map is not None

        # color map
        colors = ["black", "lightskyblue", "firebrick"]
        values = [0.0    , 0.5           , 1.0        ]
        linear = list(zip(values,colors))
        cmap   = LinearSegmentedColormap.from_list("rg", linear, N=256)
        cmap.set_bad("white", alpha=1)

        # initiate figure
        fig = plt.figure()

        # if map is included
        if include_map:

            # x and y arrays for plotting
            if scale == "linear":
                x_map_plot = x_map
                y_map_plot = y_map
            if scale == "log":
                x_map_plot = np.log10(x_map)
                y_map_plot = np.log10(y_map)

            # plot 2D detection map
            dx = np.nanmedian(np.diff(x_map_plot))
            dy = np.nanmedian(np.diff(y_map_plot))
            extent = [x_map_plot[0]-dx/2, x_map_plot[-1]+dx/2, y_map_plot[0]-dy/2, y_map_plot[-1]+dy/2]
            plt.imshow(detection_map.T, cmap=cmap, origin="lower", vmin=vmin, vmax=vmax, extent=extent, aspect="auto", interpolation="none")
        
        # if array is included
        if include_arr:

            # x and y arrays for plotting
            if scale == "linear":
                x_arr_plot = x_arr
                y_arr_plot = y_arr
            if scale == "log":
                x_arr_plot = np.log10(x_arr)
                y_arr_plot = np.log10(y_arr)
            
            # plot 1D detection array
            idx_nan = np.isnan(detection_arr)
            plt.scatter(x_arr_plot[idx_nan], y_arr_plot[idx_nan], c="w"          , cmap=cmap, vmin=vmin, vmax=vmax, marker="o", edgecolors="k")
            plt.scatter(x_arr_plot         , y_arr_plot         , c=detection_arr, cmap=cmap, vmin=vmin, vmax=vmax, marker="o", edgecolors="k")

        # plot labels
        if x_var == "P": xlabel = "$P$ [d]"
        if x_var == "a": xlabel = "$a$ [AU]"
        if y_var == "K": ylabel = "$K$ [km/s]"
        if y_var == "m": ylabel = "$m$ [$M_{\mathrm{Earth}}$]"
        if scale == "log":
            xlabel = "$\log_{10}\,$" + xlabel
            ylabel = "$\log_{10}\,$" + ylabel
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.gca().tick_params(axis="both", which="both", direction="in", top=True, right=True)
        plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
        plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
        
        # color bar
        cb = plt.colorbar(pad=0.01)
        cb.ax.tick_params(axis="both", which="both", direction="in")
        cb.ax.yaxis.set_minor_locator(AutoMinorLocator())
        if y_var == "K": cb.set_label("$K_{\mathrm{rec}}/K_{\mathrm{inj}}$")
        if y_var == "m": cb.set_label("$m_{\mathrm{rec}}/m_{\mathrm{inj}}$")

        # return figure
        return fig