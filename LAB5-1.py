import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons
import numpy as np
from scipy import signal

global_noise = None

def harmonic_with_noise(t, amplitude, frequency, phase, noise_mean, noise_std, show_noise):
    global global_noise
    harmonic_signal = amplitude * np.sin(2 * np.pi * frequency * t + phase)
    if show_noise:
        if global_noise is None:
            global_noise = np.random.normal(noise_mean, noise_std, size=len(t))
        return harmonic_signal + global_noise
    else:
        return harmonic_signal

def apply_filter(signal_data, cutoff_freq, fs, filter_type='lowpass', filter_order=4):
    nyquist = 0.5 * fs
    Wn = cutoff_freq / nyquist
    b, a = signal.butter(filter_order, Wn, btype=filter_type)
    filtered_signal = signal.filtfilt(b, a, signal_data)
    return filtered_signal

initial_amplitude = 1.0
initial_frequency = 1.0
initial_phase = 0.0
initial_noise_mean = 0.0
initial_noise_std = 0.1
initial_cutoff_frequency = 5.0

show_noise = False
t = np.linspace(0, 10, 1000)
fs = 1 / (t[1] - t[0])

fig, ax = plt.subplots()
ax.grid(True)
plt.subplots_adjust(left=0.1, bottom=0.5)

noisy_signal = harmonic_with_noise(t, initial_amplitude, initial_frequency, initial_phase, initial_noise_mean, initial_noise_std, show_noise)

line_noisy, = ax.plot(t, noisy_signal, lw=2, color='purple', label='Noisy Signal')
line_filtered, = ax.plot(t, [np.nan]*len(t), lw=2, color='orange', label='Filtered Signal')
line_original, = ax.plot(t, initial_amplitude * np.sin(2 * np.pi * initial_frequency * t + initial_phase), lw=2, color='blue', label='Original Harmonic')

axcolor = 'lightgoldenrodyellow'
ax_amplitude = plt.axes([0.1, 0.4, 0.65, 0.03], facecolor=axcolor)
ax_frequency = plt.axes([0.1, 0.35, 0.65, 0.03], facecolor=axcolor)
ax_phase = plt.axes([0.1, 0.3, 0.65, 0.03], facecolor=axcolor)
ax_noise_mean = plt.axes([0.1, 0.25, 0.65, 0.03], facecolor=axcolor)
ax_noise_std = plt.axes([0.1, 0.2, 0.65, 0.03], facecolor=axcolor)
ax_cutoff_freq = plt.axes([0.1, 0.15, 0.65, 0.03], facecolor=axcolor)

s_amplitude = Slider(ax_amplitude, 'Amplitude', 0.1, 10.0, valinit=initial_amplitude)
s_frequency = Slider(ax_frequency, 'Frequency', 0.1, 10.0, valinit=initial_frequency)
s_phase = Slider(ax_phase, 'Phase', 0.0, 2*np.pi, valinit=initial_phase)
s_noise_mean = Slider(ax_noise_mean, 'Noise Mean', -1.0, 1.0, valinit=initial_noise_mean)
s_noise_std = Slider(ax_noise_std, 'Noise STD', 0.01, 1.0, valinit=initial_noise_std)
s_cutoff_freq = Slider(ax_cutoff_freq, 'Cutoff Frequency', 0.1, 10.0, valinit=initial_cutoff_frequency)

def update(val):
    global global_noise
    amplitude = s_amplitude.val
    frequency = s_frequency.val
    phase = s_phase.val
    noise_mean = s_noise_mean.val
    noise_std = s_noise_std.val
    cutoff_freq = s_cutoff_freq.val
    show_noise = check.get_status()[0]

    if show_noise:
        global_noise = np.random.normal(noise_mean, noise_std, size=len(t))
    else:
        global_noise = None

    noisy_signal = harmonic_with_noise(t, amplitude, frequency, phase, noise_mean, noise_std, show_noise)
    line_noisy.set_ydata(noisy_signal)
    line_original.set_ydata(amplitude * np.sin(2 * np.pi * frequency * t + phase))

    if show_noise:
        filtered_signal = apply_filter(noisy_signal, cutoff_freq, fs)
        line_filtered.set_ydata(filtered_signal)
    else:
        line_filtered.set_ydata([np.nan]*len(t))

    fig.canvas.draw_idle()

for slider in [s_amplitude, s_frequency, s_phase, s_noise_mean, s_noise_std, s_cutoff_freq]:
    slider.on_changed(update)

rax = plt.axes([0.025, 0.7, 0.15, 0.1], facecolor=axcolor)
check = CheckButtons(rax, ['Show Noise'], [show_noise])
check.on_clicked(update)

resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
button_reset = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

def reset(event):
    for slider in [s_amplitude, s_frequency, s_phase, s_noise_mean, s_noise_std, s_cutoff_freq]:
        slider.reset()
    if check.get_status()[0]:
        check.set_active(0)

button_reset.on_clicked(reset)

plt.legend()
plt.show()
