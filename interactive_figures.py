#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

class InteractiveFigure:
    
    def __init__(self, title, calc):
        self.title = title
        self.calc = calc

    def __enter__(self):
        self.ax = plt.figure().gca()
        self.ax.set_title(self.title)
        self.ax.set_xlim(0, 1)
        return self.ax

    def __exit__(self, *args):
        self.ax.figure.canvas.mpl_connect('motion_notify_event', self.event_handler)
        self.ax.figure.canvas.mpl_connect('button_press_event', self.event_handler)
        self.ax._handler_reference = self

    def event_handler(self, event):
        from matplotlib.backend_bases import MouseButton
        if event.inaxes and event.button == MouseButton.LEFT:
            self.update(event.xdata, event.ydata)
            self.ax.figure.canvas.draw_idle()

class InteractiveFigurePos(InteractiveFigure):

    def __exit__(self, *args):
        self.mu_line, = self.ax.plot([], [], ".-")
        self.set_line = self.ax.axvline(-1, c="k", ls="--")
        self.update(0.5, 0)
        super().__exit__(*args)

    def update(self, x, y):
        self.set_line.set_xdata([x,x])
        xdata, ydata = self.calc(x)
        self.mu_line.set_data(xdata, ydata)

class InteractiveFigureSlope(InteractiveFigure):

    def __exit__(self, *args):
        self.mu_line, = self.ax.plot([], [], ".-")
        self.set_line, = self.ax.plot([], [], "k.--")
        self.update(0, 0)
        super().__exit__(*args)


    def update(self, x, y):
        ylim = self.ax.get_ylim()
        yc = (ylim[0] + ylim[1])/2
        slope = (y-yc)/(x-0.5)
        self.set_line.set_data([0, 0.5, 1], [yc-0.5*slope, yc, yc+0.5*slope ])
        xdata, ydata = self.calc(slope)
        self.mu_line.set_data(xdata, ydata)

