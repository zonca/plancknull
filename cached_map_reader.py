#!/usr/bin/env python

import healpy
from collections import namedtuple, deque
import numpy as np
import os.path

MapCacheEntry = namedtuple('MapCacheEntry', ('path', 'maps'))

class CachedMapReader:
    '''This class behaves like a function and reads Healpix maps.

    The class can be used as a substitute of healpy.read_map. It
    implements a caching mechanism, so that if the same map is asked
    twice, the FITS file will be loaded only once.

    The return value is either a NumPy array of pixels (if only one
    field among I, Q, U is chosen) or a tuple of NumPy arrays.

    Unlike healpy.read_map, masked pixels are NaN.'''

    def __init__(self, max_num_of_maps = 10, fields=(0, 1, 2)):
        '''Initialize the cache.

        The value of `max_num_of_maps' specifies how many maps to
        retain in the cache, and it should be chosen according to the
        available memory. The value of `fields' is analogue to the
        field accepted by healpy.read_map.'''

        self.cache = deque(maxlen=max_num_of_maps)
        self.fields = fields

    def __call__(self, path):
        '''Read a map and return the list of pixels.'''

        # Using `realpath' we can match 'foo.fits' with './foo.fits'
        real_path = os.path.realpath(path)
        matches = filter(lambda e: e.path == real_path, self.cache)
        if matches:
            # Note that we might want to move this map to the
            # beginning of the deque. For the moment, we keep
            # things as simple as possible.
            print "Cache hit"
            return matches[0].maps
        else:
            print "Cache miss"
            maps = healpy.read_map(path, field=self.fields)

            for pixels in maps:
                pixels[pixels < -1.6e+30] = np.NAN

            # Since we specified `maxlen' in the call to `deque'
            # this call to `appendleft' will automatically
            # trigger a call to `pop' if the list gets too large.
            self.cache.appendleft(MapCacheEntry(path=real_path,
                                                maps=maps))

            return maps
