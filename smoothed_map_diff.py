#!/usr/bin/env python

"""
Read a set of Healpix maps and perform the following steps:

1. Combine them linearly using user-provided weights
2. Apply a smoothing filter to the result
3. Degrade the map down to some user-specified NSIDE
4. Save the result
"""

# for interactive use maptplotlib should be already 
# inizialized, so this does not have any effect,
# as intended.
import matplotlib
matplotlib.use('AGG')

import numpy as np
import healpy as hp
import sys
import logging as log
import math
import itertools
from optparse import OptionParser
from cached_map_reader import CachedMapReader

#####################################################################


def parse_command_line():
    "Interpret the command line parameters using an OptionParser object"
    parser = OptionParser(usage="%prog [OPTS] MAP1 COEFF1 [MAP2 COEFF2] ...")

    parser.add_option('--degraded-nside', metavar='NSIDE',
                      default=32, type='int', dest='degraded_nside',
                      help='Value of NSIDE for the output map '
                      '(default: %default). Use -1 if you want to '
                      'skip this operation. ')
    parser.add_option('--smoothing-angle', metavar='DEGREES',
                      default=2.0, type='float',
                      dest='smoothing_angle',
                      help='Value of the FHWM (in degrees) to be used '
                      'by the smoother (default: %default). Set it to 0.0'
                      'to skip this operation.')
    parser.add_option('--components', '-c', default='I,Q,U', metavar='STRING',
                      dest='stokes_components',
                      help='Comma-separated list of Stokes components to'
                      'analyze (default: "%default")')
    parser.add_option('--output', '-o', metavar="FILE",
                      dest='output_file_name', default='output.fits',
                      help='Name of the output file that will contain the map'
                      ' (default: "%default")')

    return parser.parse_args()

#####################################################################


def validate_parameters(options, args):
    """Check the validity of the parameters in 'options'.

    The 'options' and "args" variables are returned by parse_command_line."""

    if options.degraded_nside <= 0:
        options.degraded_nside = -1
    else:
        # Remember that NSIDE must always be some power of 2
        nside_log2 = math.log(options.degraded_nside) / math.log(2)
        if 2 ** int(nside_log2) != options.degraded_nside:
            log.fatal("NSIDE = %d is not an integer power of 2",
                      options.degraded_nside)
            sys.exit(1)

    if options.smoothing_angle < 0.0:
        log.fatal("Invalid value for the smoothing angle (%f)",
                  options.smoothing_angle)
        sys.exit(1)

    # We want to convert a string of Stokes components like "I,q"
    # into something recognized by hp.read_map, i.e. (0, 1).
    # The use of "frozenset" removes duplicates.
    comp = frozenset([x
                      for x in options.stokes_components.upper().split(",")
                      if x != ''])
    log.info('Selected Stokes components are %s',
             ', '.join(comp))
    component_map = {'I': 0, 'Q': 1, 'U': 2}
    try:
        stokes_components = tuple([component_map[i]
                                   for i in comp])
    except KeyError:
        log.fatal('Unknown Stokes component %s in string "%s" '
                  '(available choices are %s)',
                  sys.exc_value,
                  options.stokes_components,
                  ', '.join(['"%s"' % x for x in component_map.keys()]))
        sys.exit(1)

    # Now overwrite options.stokes_components: we do not need the
    # user-provided string any longer
    options.stokes_components = stokes_components

    if len(args) == 0:
        log.fatal('You must specify at least one file name (together with its '
                  'weight) on the command line')
        sys.exit(1)

    if len(args) % 2 != 0:
        log.fatal('Wrong number of command-line parameters: each map '
                  'file name must be associated with a number (the weight)')
        sys.exit(1)

#####################################################################


def smooth_and_combine_maps(maps_and_weights,
                            output_file_name,
                            reader_function,
                            stokes_components=(0, 1, 2),
                            smoothing_angle=2.0,
                            degraded_nside=32):
    '''Combine a set of FITS maps into another.

    The maps are specified by `maps_and_weights', a list of tuples
    of the form (PATH, WEIGHT), where PATH is a string containing the
    path of the FITS file, and '''

    combined_map = []
    for map_file_name, weight in maps_and_weights:
        log.info("Reading map %s (weight %f)", map_file_name, weight)
        pixels = reader_function(map_file_name)
        if len(combined_map) == 0:
            combined_map = [x * weight for x in pixels]
        else:
            for idx in xrange(len(pixels)):
                combined_map[idx] = combined_map[idx] + weight * pixels[idx]

        # Maps with NSIDE=1024 are large, so it's better to free memory we
        # are not going to use any longer
        del pixels

    component_name = {0: 'I', 1: 'Q', 2: 'U'}
    smoothed_map = []
    if smoothing_angle > 0.0:
        log.info('Applying smoothing filter to %d Stokes maps',
                 len(combined_map))
        sm_angle = np.deg2rad(smoothing_angle)
        for idx, diff_component in enumerate(combined_map):
            nan_mask = np.isnan(diff_component)
            diff_component[nan_mask] = 0.0

            log.info('Applying a smoothing filter to %s map...',
                     component_name[stokes_components[idx]])
            smoothed_component = hp.smoothing(diff_component,
                                                  sm_angle)
            smoothed_component[nan_mask] = np.NaN

            smoothed_map.append(smoothed_component)
    else:
        smoothed_map = combined_map

    if degraded_nside > 0:
        log.info('Degrading the map to NSIDE = %d', degraded_nside)
        degraded_map = hp.ud_grade(smoothed_map, degraded_nside)

    log.info('Saving the map into %s', output_file_name)
    hp.write_map(output_file_name, degraded_map)

#####################################################################


if __name__ == "__main__":
    log.basicConfig(format='[%(asctime)s %(levelname)s '
                    + '%(filename)s:%(lineno)d] %(message)s',
                    level=log.DEBUG)

    OPTIONS, ARGS = parse_command_line()
    validate_parameters(OPTIONS, ARGS)

    # The following instruction splits ARGS into a list of
    # 2-tuples, e.g.
    # ['a', 0.1, 'b', 0.2] => [('a', 0.1), ('b', 0.2)]
    MAPS_AND_WEIGHTS = [(x[0], float(x[1]))
                        for x in zip(*[itertools.islice(ARGS, i, None, 2)
                                       for i in (0, 1)])]

    READER = CachedMapReader(fields=OPTIONS.stokes_components)
    smooth_and_combine_maps(MAPS_AND_WEIGHTS,
                            OPTIONS.output_file,
                            READER,
                            OPTIONS.stokes_components,
                            OPTIONS.smoothing_angle,
                            OPTIONS.degraded_nside)
