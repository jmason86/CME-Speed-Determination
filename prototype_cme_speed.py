# -*- coding: utf-8 -*-
"""
========================================
Prototyping CME speed determination tool
========================================
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button

import astropy.units as u
from astropy.coordinates import SkyCoord

import sunpy.map
import sunpy.coordinates.wcs_utils
from sunpy.net import Fido, attrs as a

from matplotlib.lines import Line2D

import warnings
warnings.filterwarnings('ignore')

# Find some data to download
stereo = (a.vso.Source('STEREO_B') &
          a.Instrument('EUVI') &
          a.Time('2011-01-01T00:14:00', '2011-01-01T03:00:00'))

aia = (a.Instrument('SWAP') &
       a.vso.Sample(24 * u.hour) &
       a.Time('2011-01-01T00:14:00', '2011-01-01T03:00:00'))

wave = a.Wavelength(17 * u.nm, 18 * u.nm)
result = Fido.search(wave, aia | stereo)

# Let's inspect the result
print(result)

# and download the files
downloaded_files = Fido.fetch(result)
print(downloaded_files)

# Let's create a dictionary with the two maps, which we crop to full disk.
#maps = {m.name.split(' ', 1)[0]: m for m in sunpy.map.Map(downloaded_files)}
maps = {m.name.split(' ', 1)[1]: m for m in sunpy.map.Map(downloaded_files)}  # TODO: load the two instruments into different lists or dictionaries or maps (if those can be concatenated) instead of into a single dictionary

# Prepare button
icon_next = plt.imread('https://i.imgur.com/zk94EAK.png')
icon_done = plt.imread("https://i.imgur.com/jHGPyFy.png")

##############################################################################
# Now, let's plot both maps

fig = plt.figure(figsize=(10, 4))
ax1 = fig.add_subplot(1, 2, 1, projection=maps['SWAP'])
maps['SWAP'].plot(axes=ax1)
ax1_map_name = ax1.axes.title.get_text().split(' ', 1)[0]

ax2 = fig.add_subplot(1, 2, 2, projection=maps['EUVI-B'])
maps['EUVI-B'].plot(axes=ax2)

# Setup initial interaction parameters
is_last_map = False
line_of_sight_is_defined = False


def onclick(event):
    global clicked_map, other_map, line_of_sight_is_defined
    clicked_map = which_map_clicked(event)
    other_map = which_is_other_map()
    clicked_skycoord = get_clicked_skycoord(event)
    draw_clicked_circle(clicked_skycoord)

    if not line_of_sight_is_defined:
        translate_skycoord_to_other_map(clicked_skycoord)
        draw_translated_line()
        line_of_sight_is_defined = True

    closeout_clicks(event)

    return True


def which_map_clicked(event):
    instrument_name = event.inaxes.title.get_text().split(' ', 1)[0]
    return instrument_name


def which_is_other_map():
    return np.setdiff1d(list(maps.keys()), [clicked_map])[0]


def get_clicked_skycoord(event):
    ix, iy = event.xdata, event.ydata
    clicked_skycoord = maps[clicked_map].pixel_to_world(ix * u.pix, iy * u.pix)
    return clicked_skycoord


def translate_skycoord_to_other_map(clicked_skycoord):
    global line_coords
    point_to_line = clicked_skycoord.realize_frame(clicked_skycoord.spherical * np.linspace(0.9, 1.1, 1e6) * u.AU)
    line_coords = point_to_line.transform_to(maps[other_map].coordinate_frame)


def draw_clicked_circle(clicked_skycoord):
    if ax1_map_name == clicked_map:
        ax1.plot_coord(clicked_skycoord, color='g', marker='o', fillstyle='none')
    else:
        ax2.plot_coord(clicked_skycoord, color='g', marker='o', fillstyle='none')


def draw_translated_line():
    if ax1_map_name == other_map:
        ax_lim = ax1.axis()
        ax1.plot_coord(line_coords, color='g', picker=5)
        ax1.axis(ax_lim)
    else:
        ax_lim = ax2.axis()
        ax2.plot_coord(line_coords, color='g', picker=5)
        ax2.axis(ax_lim)
    plt.draw()


def closeout_clicks(event):
    if line_of_sight_is_defined:
        fig.canvas.mpl_disconnect(cid1)


def pick_los_point(event):
    if isinstance(event.artist, Line2D):
        index = int(np.median(event.ind))
        skycoord_3d = line_coords[index]
        draw_3d_points(skycoord_3d)

        print(skycoord_3d)  # TODO: need to pass this to another function that'll do something with it
    return True


def draw_3d_points(skycoord_3d):
    skycoord_3d_in_other_map = skycoord_3d.transform_to(maps[other_map].coordinate_frame)

    if ax1_map_name == other_map:
        ax2.plot_coord(skycoord_3d, color='blue', marker='o')
        ax1.plot_coord(skycoord_3d_in_other_map, color='blue', marker='o')
    else:
        ax1.plot_coord(skycoord_3d, color='blue', marker='o')
        ax2.plot_coord(skycoord_3d_in_other_map, color='blue', marker='o')

    plt.draw()


def next_map_clicked(event):
    if is_last_map:
        button_axes.images[0].set_data(icon_done)
        fig.canvas.draw_idle()
    else:
        load_new_maps()


def load_new_maps():
    pass


cid1 = fig.canvas.mpl_connect('button_press_event', onclick)
fig.canvas.mpl_connect('pick_event', pick_los_point)

button_axes = plt.axes([0.83, 0.04, 0.22, 0.22])
start_button = Button(button_axes, '', image=icon_next)
start_button.on_clicked(next_map)

plt.show()
