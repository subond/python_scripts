from scipy.signal import butter, lfilter, freqz
import matplotlib.pyplot as plt
import numpy as np

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq

    if np.logical_or(low == 1.0 , high == 1.0):
	print 'your nyquist frequency is equal to one of your band-pass filter bounds. This will cause significant errors. low/nyq , high/nyq :', low, high

    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5, axis_to_filter=0):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data, axis=axis_to_filter)
    return y

def convert_frequencies_to_hertz(delta_time, units ):
    
    if (units == 'days'):
	freq = 1./(delta_time*24.*60.*60.)

    if (units == 'years'):
	freq = 1./(delta_time*365.25*24.*60.*60.)

    if (units == 'hours'):
	freq = 1./(delta_time*60.*60.)
    
    return freq

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import freqz

    # Sample rate and desired cutoff frequencies (in Hz).
    fs = convert_frequencies_to_hertz( 1, 'days')
    lowcut = convert_frequencies_to_hertz( 12, 'days')
    highcut = convert_frequencies_to_hertz( 2., 'days')

#    fs     = 5000.
#    lowcut =  500.
#    highcut = 1250.

    print fs, lowcut, highcut


    plt.figure(1)
    plt.clf()
    for order in [3, 6, 9]:
        b, a = butter_bandpass(lowcut, highcut, fs, order=order)
        w, h = freqz(b, a, worN=2000)
        plt.plot((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)],
             '--', label='sqrt(0.5)')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Gain')
    plt.grid(True)
    plt.legend(loc='best')



#    # Filter a noisy signal.
#    T = 0.05
#    nsamples = T * fs
#    t = np.linspace(0, T, nsamples, endpoint=False)
#    a = 0.02
#    f0 = 600.0
#    x = 0.1 * np.sin(2 * np.pi * 1.2 * np.sqrt(t))
#    x += 0.01 * np.cos(2 * np.pi * 312 * t + 0.1)
#    x += a * np.cos(2 * np.pi * f0 * t + .11)
#    x += 0.03 * np.cos(2 * np.pi * 2000 * t)
#    plt.figure(2)
#    plt.clf()
#    plt.plot(t, x, label='Noisy signal')
#
#    y = butter_bandpass_filter(x, lowcut, highcut, fs, order=6)
#    plt.plot(t, y, label='Filtered signal (%g Hz)' % f0)
#    plt.xlabel('time (seconds)')
#    plt.hlines([-a, a], 0, T, linestyles='--')
#    plt.grid(True)
#    plt.axis('tight')
#    plt.legend(loc='upper left')

    plt.show()
