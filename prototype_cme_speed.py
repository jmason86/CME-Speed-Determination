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

# Define Sun-Earth Line (SEL)

# Find some data to download
stereo = (a.vso.Source('STEREO_B') &
          a.Instrument('EUVI') &
          a.Time('2011-01-01T00:14:00', '2011-01-01T03:00:00'))

aia = (a.Instrument('SWAP') &
       a.vso.Sample(24 * u.hour) &
       a.Time('2011-01-01T00:14:00', '2011-01-01T03:00:00'))

wave = a.Wavelength(17 * u.nm, 18 * u.nm)
result_left = Fido.search(wave, aia)
result_right = Fido.search(wave, stereo)

# Let's inspect the result
print(result_left)
print(result_right)

# and download the files
downloaded_files_left = Fido.fetch(result_left, progress=False)
downloaded_files_right = Fido.fetch(result_right, progress=False)

# Let's create a dictionary with the two maps, which we crop to full disk.
#maps = {m.name.split(' ', 1)[0]: m for m in sunpy.map.Map(downloaded_files)}
#maps = {m.name.split(' ', 1)[1]: m for m in sunpy.map.Map(downloaded_files)}  # TODO: load the two instruments into different lists or dictionaries or maps (if those can be concatenated) instead of into a single dictionary
maps_left = sunpy.map.Map(downloaded_files_left, sequence=True)
maps_right = sunpy.map.Map(downloaded_files_right, sequence=True)


# Prepare button
icon_next = plt.imread('https://i.imgur.com/zk94EAK.png')
icon_done = plt.imread("https://i.imgur.com/jHGPyFy.png")

##############################################################################
# Now, let's plot both maps

fig = plt.figure(figsize=(10, 4))
ax_left = fig.add_subplot(1, 2, 1, projection=maps_left[0])
maps_left[0].plot(axes=ax_left)

ax_right = fig.add_subplot(1, 2, 2, projection=maps_right[0])
maps_right[0].plot(axes=ax_right)

# Setup initial interaction parameters
is_last_map = False
map_sequence_index = 0
line_of_sight_is_defined = False


def onclick(event):
    global line_of_sight_is_defined
    which_map_clicked(event)
    clicked_skycoord = get_clicked_skycoord(event)
    draw_clicked_circle(clicked_skycoord)

    if not line_of_sight_is_defined:
        translate_skycoord_to_other_map(clicked_skycoord)
        draw_translated_line()
        line_of_sight_is_defined = True

    closeout_clicks(event)

    return True


def which_map_clicked(event):
    global clicked_map, other_map
    if event.inaxes.colNum == 0:
        clicked_map = maps_left[map_sequence_index]
        other_map = maps_right[map_sequence_index]
    else:
        clicked_map = maps_right[map_sequence_index]
        other_map = maps_left[map_sequence_index]


def get_clicked_skycoord(event):
    ix, iy = event.xdata, event.ydata
    clicked_skycoord = clicked_map.pixel_to_world(ix * u.pix, iy * u.pix)
    return clicked_skycoord


def translate_skycoord_to_other_map(clicked_skycoord):
    global line_coords
    point_to_line = clicked_skycoord.realize_frame(clicked_skycoord.spherical * np.linspace(0.9, 1.1, 1e6) * u.AU)
    line_coords = point_to_line.transform_to(other_map.coordinate_frame)


def draw_clicked_circle(clicked_skycoord):
    if clicked_map == maps_left[map_sequence_index]:
        ax_left.plot_coord(clicked_skycoord, color='g', marker='o', fillstyle='none')
    else:
        ax_right.plot_coord(clicked_skycoord, color='g', marker='o', fillstyle='none')


def draw_translated_line():
    if clicked_map == maps_right[map_sequence_index]:
        ax_lim = ax_left.axis()
        ax_left.plot_coord(line_coords, color='g', picker=5)
        ax_left.axis(ax_lim)
    else:
        ax_lim = ax_right.axis()
        ax_right.plot_coord(line_coords, color='g', picker=5)
        ax_right.axis(ax_lim)
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
    skycoord_3d_in_other_map = skycoord_3d.transform_to(other_map.coordinate_frame)

    if clicked_map == maps_right[map_sequence_index]:
        ax_right.plot_coord(skycoord_3d, color='blue', marker='o')
        ax_left.plot_coord(skycoord_3d_in_other_map, color='blue', marker='o')
    else:
        ax_left.plot_coord(skycoord_3d, color='blue', marker='o')
        ax_right.plot_coord(skycoord_3d_in_other_map, color='blue', marker='o')

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
button = Button(button_axes, '', image=icon_next)
button.on_clicked(next_map_clicked)

plt.show()
