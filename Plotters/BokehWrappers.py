#! /usr/bin/env python

from bokeh.plotting import figure, output_file, show
from bokeh.models import LabelSet, Label, Title

def plot_single_drift_simulation(xdata, ydata, output_filename):
    output_file(output_filename)

    p = figure(
        plot_width=600, plot_height=600,
        x_axis_label='Generations', y_axis_label='Frequency',
        y_range=(0,1)
    )

    # add a line renderer
    p.line(xdata, ydata, line_width=1)

    title = Title(text="The Effect of Genetic Drift", align="center")
    p.add_layout(title, "above")

    show(p)


if __name__ == '__main__':
    main()
