# contains de functions presentend in the jupyter notebook
# Simone Mencarelli Sep 22.

import matplotlib.pyplot as plt
from design_functions import *

## private functions
def nadir_return_plotter_rg(ax, h, prf_min, prf_max, c=299792458, prf_resolution=1):
    """
    plots the nadir returns for the timing diagram on the passed axis
    :param ax: matplotlib axis object
    :param h: platform height
    :param prf_min: minimum prf
    :param prf_max: maximum prf
    :param c: c: optional default speed of light
    :param prf_resolution: optional default 1 HZ, PRF axis resolution
    :return:
    """
    # horizon slant range
    rmax = np.floor(np.sqrt((h + 6371e3) ** 2 - 6371e3 ** 2))
    # step 1: find n min and n max
    n_min = np.floor(2 * h * prf_min / c)
    n_max = np.floor(2 * rmax * prf_max / c)
    nn = np.arange(n_min, n_max)
    #nn = np.arange(1, 20)

    # step 2: find m min and m max (how many times can H be divided by c/2prf)
    m_min = np.floor(prf_min * 2 * h / c)
    m_max = np.ceil(prf_max * 2 * h / c)
    mm = np.arange(m_min, m_max)
    # for every m
    for m in mm:
        # step 3 find the prf axis of validity
        ppmin = m * c / (2 * h) + prf_resolution
        ppmax = (m + 1) * c / (2 * h) - prf_resolution
        if ppmax > ppmin:
            prff = np.linspace(ppmin, ppmax, int(1 + np.abs((ppmax - ppmin)) / prf_resolution))
        else:
            prff = np.zeros(1)
        # for every n
        for n in nn:
            # step 4: find the slant range of the nadir return
            R_nad = h % (c / (2 * prff)) + n * (c / (prff * 2))
            R_nad = np.where(R_nad > rmax, rmax, R_nad)
            # step 5: find the ground range of the nadir return
            rg_nad, th = range_slant_to_ground(R_nad, h)
            # step 6: plot the damn thing
            ax.plot(prff, rg_nad / 1000, 'k')


def range_transmit_event(n, dutycycle, prf, h=500e3, c=299792458):
    """
    finds a list of n curves for the enveloped slant ranges associated to transmit events (bot and eot)
    :param n: maximum order of the return
    :param dutycycle: duty cycle
    :param prf: prf axis
    :param c: optional default speed of light
    :return:
    """
    rmax = np.floor(np.sqrt((h + 6371e3) ** 2 - 6371e3 ** 2))
    transmits = []
    for nn in n:
        bot = nn * (c / (2 * prf)) - (c * dutycycle / (2 * prf))
        eot = nn * (c / (2 * prf)) + (c * dutycycle / (2 * prf))
        bot = np.where(bot < rmax, bot, rmax)
        eot = np.where(eot < rmax, eot, rmax)
        tev = (bot, eot)
        transmits.append(tev)
    return transmits

## interface
def time_diagram_plotter(ax, prf, dutycycle, h=500e3, re=6371e3, c=299792458):
    """
    plots a timing diagram over the given axis
    :param ax: matplotlib axis
    :param prf: pulse repetition frequency
    :param dutycycle: duty cycle of the radar
    :param h: satellite height
    :param re: earth radius
    :param c: optional default speed of light
    :return:
    """
    nadir_return_plotter_rg(ax, h, prf[0], prf[-1], prf_resolution=5)
    n = np.arange(0, 200)
    transmit_events = range_transmit_event(n, dutycycle, prf)
    # conversion to ground range
    transmit_events_rg = []
    transmit_events_theta = []
    for nn in range(len(n)):
        transmit_event_rg_bot, transmit_event_theta_bot = range_slant_to_ground(transmit_events[nn][0])
        transmit_event_rg_eot, transmit_event_theta_eot = range_slant_to_ground(transmit_events[nn][1])
        transmit_events_rg.append((transmit_event_rg_bot, transmit_event_rg_eot))
        transmit_events_theta.append((transmit_event_theta_bot, transmit_event_theta_eot))

    transmit_events_rg = np.where(np.isnan(transmit_events_rg), 0, transmit_events_rg)
    for nn in range(len(n)):
        ax.fill_between(prf, transmit_events_rg[nn][1] / 1000, transmit_events_rg[nn][0] / 1000, color='orange')
    ax.set_xlabel('PRF [Hz]')