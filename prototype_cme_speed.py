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
from astropy import coordinates

import sunpy.map
import sunpy.coordinates.wcs_utils
from sunpy.net import Fido, attrs as a

from matplotlib.lines import Line2D

import warnings
warnings.filterwarnings('ignore')

start_time = '2012-07-23T02:30:00'
end_time = '2012-07-23T05:30:00'

# Find some data to download
off_sun_earth_line_imager = (a.vso.Source('STEREO_A') &
                             a.Instrument('EUVI') &
                             a.Time(start_time, end_time))

sun_earth_line_imager = (a.Instrument('SWAP') &
                         #a.Detector('C2') &
                         a.Time(start_time, end_time))

wave = a.Wavelength(17 * u.nm, 18 * u.nm)
result_left = Fido.search(wave, sun_earth_line_imager)
result_right = Fido.search(wave, off_sun_earth_line_imager)

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


def temporally_align_map_sequences():
    global maps_left, maps_right
    maps_driver, maps_to_subsample = which_maps_drive_vs_drop()
    maps_subsampled = drop_extraneous_maps(maps_driver, maps_to_subsample)
    if maps_driver == maps_left:
        maps_right = maps_subsampled
    else:
        maps_left = maps_subsampled


def which_maps_drive_vs_drop():
    maps_driver = min(maps_left, maps_right, key=len)
    maps_to_subsample = max(maps_left, maps_right, key=len)
    return maps_driver, maps_to_subsample


def drop_extraneous_maps(maps_driver, maps_to_subsample):
    indices_to_keep = []
    for map_driver in maps_driver:
        delta_ts = []
        for i in range(len(maps_to_subsample)):
            delta_ts.append(map_driver.date - maps_to_subsample[i].date)
        indices_to_keep.append(delta_ts.index(min(delta_ts, key=abs)))
    tmp = [maps_to_subsample[i] for i in indices_to_keep]
    return sunpy.map.MapSequence(tmp)

temporally_align_map_sequences()

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
    if is_not_plot_clicked(event):
        return False

    which_map_clicked(event)
    clicked_skycoord = get_clicked_skycoord(event)

    if not line_of_sight_is_defined:
        draw_clicked_circle(clicked_skycoord)
        translate_skycoord_to_other_map(clicked_skycoord)
        draw_translated_line()
        line_of_sight_is_defined = True

    closeout_clicks(event)

    return True


def is_not_plot_clicked(event):
    return not hasattr(event.inaxes, 'colNum')


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


def done_clicked(event):
    fig.canvas.mpl_disconnect(cid1)
    fig.canvas.mpl_disconnect(cid2)
    print(skycoord_3d_array)  # eventually want to return this as the main return of this program?

    compute_kinematics()
    plot_kinematics()


def compute_kinematics():
    global distances, speeds, accelerations, delta_t_sec
    distances = compute_distances()
    delta_t_sec = compute_delta_time()
    speeds = compute_speeds(delta_t_sec)
    accelerations = compute_accelerations(delta_t_sec)
    print(distances)
    print(speeds)
    print(accelerations)


def compute_distances():
    sun_coord = coordinates.get_sun(maps_left[0].date)
    return [skycoord_3d.separation_3d(sun_coord) for skycoord_3d in skycoord_3d_array]


def compute_delta_time():
    obstime = [skycoord.obstime for skycoord in skycoord_3d_array]
    delta_t = [t - obstime[0] for t in obstime]
    return [dt.sec for dt in delta_t]


def compute_speeds(delta_t_sec):
    distance_values = [(distance.to(u.km)).value for distance in distances]
    return np.gradient(distance_values, delta_t_sec) * (u.km / u.second)


def compute_accelerations(delta_t_sec):
    return np.gradient(speeds.value, delta_t_sec) * (u.km / u.second / u.second)


def plot_kinematics():
    fig2 = plt.figure(figsize=(10, 10))
    ax_d = fig2.add_subplot(3, 1, 1)
    plt.plot(delta_t_sec, [d.value for d in distances], axes=ax_d)
    ax_s = fig2.add_subplot(3, 1, 2)
    plt.plot(delta_t_sec, speeds.value, axes=ax_s)
    ax_a = fig2.add_subplot(3, 1, 3)
    plt.plot(delta_t_sec, accelerations.value, axes=ax_a)

    ax_d.set_ylabel('height [R$_{\odot}$)]')
    ax_s.set_ylabel('speed [km s$^{-1}$]')
    ax_a.set_ylabel('acceleration [km s$^{-2}$]')
    ax_a.set_xlabel('time [seconds since start]')



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
