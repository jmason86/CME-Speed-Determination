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
result_left = Fido.search(wave, aia)
result_right = Fido.search(wave, stereo)

# and download the files
downloaded_files_left = Fido.fetch(result_left, progress=False)
downloaded_files_right = Fido.fetch(result_right, progress=False)

# Load the maps into two arrays that'll be plotted side by side, with left/right corresponding to spacecraft position


def which_map_on_left():
    global maps_left, maps_right
    if maps_left[0].wcs.heliographic_observer.lon.min() > maps_right[0].wcs.heliographic_observer.lon.min():
        maps_left, maps_right = maps_right, maps_left


maps_left = sunpy.map.Map(downloaded_files_left, sequence=True)
maps_right = sunpy.map.Map(downloaded_files_right, sequence=True)
which_map_on_left()

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

# Main return value
skycoord_3d_array = []  # or np.array?


def onclick(event):
    global line_of_sight_is_defined
    which_map_clicked(event)
    clicked_skycoord = get_clicked_skycoord(event)

    if not line_of_sight_is_defined:
        draw_clicked_circle(clicked_skycoord)
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
    point_to_line = clicked_skycoord.realize_frame(clicked_skycoord.spherical * np.linspace(0, 1000, 1e6) * u.solRad)
    line_coords = point_to_line.transform_to(other_map.coordinate_frame)


def draw_clicked_circle(clicked_skycoord):
    if clicked_map == maps_left[map_sequence_index]:
        ax_left.plot_coord(clicked_skycoord, color='g', marker='o', markersize=8, fillstyle='none')
    else:
        ax_right.plot_coord(clicked_skycoord, color='g', marker='o', markersize=8, fillstyle='none')
    plt.draw()  # FIXME: Figure out why this only works for the first map in the sequence


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
    pass


def pick_los_point(event):
    if isinstance(event.artist, Line2D):
        index = int(np.median(event.ind))
        skycoord_3d = SkyCoord(line_coords[index])
        skycoord_3d_array.append(skycoord_3d.transform_to(sunpy.coordinates.frames.HeliographicStonyhurst))
        print(skycoord_3d_array[-1])
        draw_3d_points(skycoord_3d)
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
    global map_sequence_index, is_last_map, line_of_sight_is_defined

    if not is_last_map:
        map_sequence_index += 1
        if map_sequence_index == min(len(maps_left), len(maps_right)) - 1:
            is_last_map = True

        load_new_maps()
        clear_clicked_annotations()
        line_of_sight_is_defined = False


def load_new_maps():
    maps_left[map_sequence_index].plot(axes=ax_left)
    maps_right[map_sequence_index].plot(axes=ax_right)


def clear_clicked_annotations():
    del ax_left.lines[:]
    del ax_right.lines[:]
    plt.draw()
    pass


def done_clicked(event):
    fig.canvas.mpl_disconnect(cid1)
    fig.canvas.mpl_disconnect(cid2)
    print(skycoord_3d_array)  # eventually want to return this as the main return of this program?


# Set up user click interactions
cid1 = fig.canvas.mpl_connect('button_press_event', onclick)
cid2 = fig.canvas.mpl_connect('pick_event', pick_los_point)

# Set up buttons
icon_next = plt.imread('https://i.imgur.com/4bu7tvv.png')
ax_button_next = plt.axes([0.89, 0.10, 0.10, 0.10])
button_next = Button(ax_button_next, '', image=icon_next)
button_next.on_clicked(next_map_clicked)

icon_done = plt.imread("https://i.imgur.com/JBazCVv.png")
ax_button_done = plt.axes([0.89, 0.01, 0.10, 0.10])
button_done = Button(ax_button_done, '', image=icon_done)
button_done.on_clicked(done_clicked)

plt.show()
