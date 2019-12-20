# -*- coding: utf-8 -*-
"""
========================================
Prototyping CME speed determination tool
========================================
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox,  CheckButtons
import matplotlib.colors as colors

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
end_time = '2012-07-23T03:00:00'
#end_time = '2012-07-23T05:30:00'

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

# Plot both maps
fig = plt.figure(figsize=(10, 4))
ax_left = fig.add_subplot(1, 2, 1, projection=maps_left[0])
maps_left[0].plot_settings['cmap'] = plt.get_cmap('Greys_r')
maps_left[0].plot(axes=ax_left)

ax_right = fig.add_subplot(1, 2, 2, projection=maps_right[0])
maps_right[0].plot_settings['cmap'] = plt.get_cmap('Greys_r')
maps_right[0].plot(axes=ax_right)

# Indicate how many maps there are and which this is (e.g., 1/12)
text_n_maps = plt.text(0.94, 0.2, '1/{}'.format(len(maps_left)), horizontalalignment='center', transform=fig.transFigure)

# Prepare kinematic profile plot
fig2, (ax_d, ax_s, ax_a) = plt.subplots(3, figsize=(10, 10))
fig2.suptitle('Radial Kinematic Profiles')
line_d, = ax_d.plot([], [], 'o-')
ax_d.set_ylabel('height [R$_{\odot}$)]')
line_s, = ax_s.plot([], [], 'o-')
ax_s.set_ylabel('speed [km s$^{-1}$]')
text_mean_speed = plt.text(0.02, 0.9, 'mean speed = [] km/s', transform=ax_s.transAxes)
line_a, = ax_a.plot([], [], 'o-')
ax_a.set_ylabel('acceleration [km s$^{-2}$]')
ax_a.set_xlabel('time [seconds since start]')

# Setup initial interaction parameters
is_last_map = False
map_sequence_index = 0
line_of_sight_is_defined = False
is_checked_base_difference = False
is_checked_running_difference = False

# Main return value
skycoord_3d_array = []


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

        compute_kinematics()
        plot_kinematics()
        put_maps_figure_back_in_focus()
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
        update_map_counter()

        line_of_sight_is_defined = False


def load_new_maps():
    if is_checked_base_difference:
        difference_map(base=True)
    elif is_checked_running_difference:
        difference_map()
    else:
        maps_left[map_sequence_index].plot_settings['cmap'] = plt.get_cmap('Greys_r')
        maps_right[map_sequence_index].plot_settings['cmap'] = plt.get_cmap('Greys_r')
        maps_left[map_sequence_index].plot(axes=ax_left)
        maps_right[map_sequence_index].plot(axes=ax_right)
        plt.draw()


def clear_clicked_annotations():
    del ax_left.lines[:]
    del ax_right.lines[:]
    plt.draw()


def update_map_counter():
    text_n_maps.set_text('{0}/{1}'.format(map_sequence_index + 1, len(maps_left)))


def difference_clicked(label):
    global is_checked_base_difference, is_checked_running_difference
    if label == 'Base':
        is_checked_base_difference = not is_checked_base_difference
        if is_checked_base_difference:
            difference_map(base=True)
        else:
            load_new_maps()
    elif label == 'Running':
        is_checked_running_difference = not is_checked_running_difference
        if is_checked_running_difference:
            difference_map()
        else:
            load_new_maps()


def difference_map(base=False):
    if map_sequence_index < 1:
        print('Click Next Map first. Difference is not defined at first time.')
        return

    if base:  # Base difference
        index = map_sequence_index
    else:  # Running difference
        index = 1

    diff_left = maps_left[map_sequence_index].data - maps_left[map_sequence_index - index].data
    diff_right = maps_right[map_sequence_index].data - maps_right[map_sequence_index - index].data
    header_left = maps_left[map_sequence_index - index].fits_header
    header_right = maps_right[map_sequence_index - index].fits_header

    diff_left_map = sunpy.map.Map(diff_left, header_left)
    diff_right_map = sunpy.map.Map(diff_right, header_right)

    diff_left_map.plot(axes=ax_left,
                       cmap=plt.get_cmap('Greys'),
                       norm=colors.Normalize(vmin=-50, vmax=50))
    diff_right_map.plot(axes=ax_right,
                        cmap=plt.get_cmap('Greys'),
                        norm=colors.Normalize(vmin=-50, vmax=50))
    plt.draw()


def power_clicked(text):
    scaling = eval(text)
    maps_left[map_sequence_index].plot_settings['norm'] = colors.PowerNorm(gamma=scaling)
    maps_right[map_sequence_index].plot_settings['norm'] = colors.PowerNorm(gamma=scaling)

    maps_left[map_sequence_index].plot(axes=ax_left)
    maps_right[map_sequence_index].plot(axes=ax_right)
    plt.draw()


def log_clicked(text):
    scaling = eval(text)
    if scaling >= maps_left[map_sequence_index].max():
        print('The input minimum for log scaling is too large for the left map. Try a smaller number.')
        return
    if scaling >= maps_right[map_sequence_index].max():
        print('The input minimum for log scaling is too large for the right map. Try a smaller number.')
        return

    maps_left[map_sequence_index].plot_settings['norm'] = colors.LogNorm(scaling, maps_left[map_sequence_index].max())
    maps_right[map_sequence_index].plot_settings['norm'] = colors.LogNorm(scaling, maps_right[map_sequence_index].max())

    maps_left[map_sequence_index].plot(axes=ax_left)
    maps_right[map_sequence_index].plot(axes=ax_right)
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
    speeds = compute_speeds()
    accelerations = compute_accelerations()
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


def compute_speeds():
    if map_sequence_index >= 1:
        distance_values = [(distance.to(u.km)).value for distance in distances]
        return np.gradient(distance_values, delta_t_sec) * (u.km / u.second)
    else:
        return None


def compute_accelerations():
    if map_sequence_index >= 2:
        return np.gradient(speeds.value, delta_t_sec) * (u.km / u.second / u.second)
    else:
        return None


def plot_kinematics():
    if distances is not None:
        line_d.set_data(delta_t_sec, [d.value for d in distances])
        ax_d.relim()
        ax_d.autoscale_view()
    if speeds is not None:
        line_s.set_data(delta_t_sec, speeds.value)
        ax_s.relim()
        ax_s.autoscale_view()
        text_mean_speed.set_text('mean speed = {0:.0f} km/s'.format(np.mean(speeds.value)))
    if accelerations is not None:
        line_a.set_data(delta_t_sec, accelerations.value)
        ax_a.relim()
        ax_a.autoscale_view()

    plt.show(block=False)
    fig2.canvas.draw_idle()  # misnomer function -- actually triggers draw (.draw() doesn't work)


def put_maps_figure_back_in_focus():
    plt.figure(fig.number)
    plt.get_current_fig_manager().show()


# Set up user click interactions
cid1 = fig.canvas.mpl_connect('button_press_event', onclick)
cid2 = fig.canvas.mpl_connect('pick_event', pick_los_point)

# Set up text inputs for intensity scaling
plt.figure(fig.number)  # Ensures the right figure will receive the textboxes/buttons
text_scaling = plt.text(0.94, 0.52, 'Intensity Scaling', horizontalalignment='center', transform=fig.transFigure)
ax_textbox_power = plt.axes([0.96, 0.40, 0.03, 0.10])
textbox_power = TextBox(ax_textbox_power, 'x$^n$', initial='1')
textbox_power.on_submit(power_clicked)

ax_textbox_log = plt.axes([0.96, 0.30, 0.03, 0.10])
textbox_log = TextBox(ax_textbox_log, 'log(n, max)', initial='N/A')
textbox_log.on_submit(log_clicked)


# Set up buttons
text_difference = plt.text(0.94, 0.78, 'Diff Image', horizontalalignment='center', transform=fig.transFigure)
ax_button_diff = plt.axes([0.89, 0.60, 0.10, 0.15])
diff_labels = ['Base', 'Running']
checkbox_diff = CheckButtons(ax_button_diff, diff_labels)
checkbox_diff.on_clicked(difference_clicked)

icon_next = plt.imread('https://i.imgur.com/4bu7tvv.png')
ax_button_next = plt.axes([0.89, 0.10, 0.10, 0.10])
button_next = Button(ax_button_next, '', image=icon_next)
button_next.on_clicked(next_map_clicked)

icon_done = plt.imread("https://i.imgur.com/JBazCVv.png")
ax_button_done = plt.axes([0.89, 0.01, 0.10, 0.10])
button_done = Button(ax_button_done, '', image=icon_done)
button_done.on_clicked(done_clicked)

plt.show()
