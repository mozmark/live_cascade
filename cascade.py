"""This module implements a bloom filter cascade data structure.

Requires the bitarray library: http://pypi.python.org/pypi/bitarray/
"""
from __future__ import absolute_import
import math
import hashlib
import copy
from struct import unpack, pack, calcsize

try:
    import bitarray
except ImportError:
    raise ImportError('cascade requires bitarray')

def make_hashfuncs(num_slices, num_bits):
    if num_bits >= (1 << 31):
        fmt_code, chunk_size = 'Q', 8
    elif num_bits >= (1 << 15):
        fmt_code, chunk_size = 'I', 4
    else:
        fmt_code, chunk_size = 'H', 2
    total_hash_bits = 8 * num_slices * chunk_size
    if total_hash_bits > 384:
        hashfn = hashlib.sha512
    elif total_hash_bits > 256:
        hashfn = hashlib.sha384
    elif total_hash_bits > 160:
        hashfn = hashlib.sha256
    elif total_hash_bits > 128:
        hashfn = hashlib.sha1
    else:
        hashfn = hashlib.md5

    fmt = fmt_code * (hashfn().digest_size // chunk_size)
    num_salts, extra = divmod(num_slices, len(fmt))
    if extra:
        num_salts += 1
    salts = tuple(hashfn(hashfn(pack('I', i)).digest()) for i in range(0, num_salts))

    def _hash_maker(key):
        if isinstance(key, str):
            key = key.encode('utf-8')
        else:
            key = str(key).encode('utf-8')
        i = 0
        for salt in salts:
            h = salt.copy()
            h.update(key)
            d = h.digest()
            for uint in unpack(fmt, d):
                try:
                    yield uint % num_bits
                    i += 1
                    if i >= num_slices:
                        return
                except Exception as e:
                    print("something went wrong!", e)

    return _hash_maker, hashfn

class Bucket:
    def __init__(self):
        self.includes = set()
        self.excludes = set()

    def include(self, item):
        related_exclusions = set()
        if len(self.includes) == 0:
            related_exclusions.update(self.excludes)
        self.includes.add(item)
        return related_exclusions

    def exclude(self, item):
        related_inclusions = set()
        if len(self.excludes) == 0:
            related_inclusions.update(self.includes)
        self.excludes.add(item)
        return related_inclusions

    def get_true_negatives(self):
        if not len(self.includes) > 0:
            return self.excludes
        return []

    def __repr__(self):
        return self.__str__()

    def __str__(self):
        return "includes: %s, excludes %s" % (self.includes, self.excludes)

class CascadeLayer:
    def __init__(self, capacity, error_rate=0.1, depth=0):
        if not (0 < error_rate < 1):
            raise ValueError("Error_Rate must be between 0 and 1.")
        if not capacity > 0:
            raise ValueError("Capacity must be > 0")
        # given M = num_bits, k = num_slices, P = error_rate, n = capacity
        #       k = log2(1/P)
        # solving for m = bits_per_slice
        # n ~= M * ((ln(2) ** 2) / abs(ln(P)))
        # n ~= (k * m) * ((ln(2) ** 2) / abs(ln(P)))
        # m ~= n * abs(ln(P)) / (k * (ln(2) ** 2))
        num_slices = int(math.ceil(math.log(1.0 / error_rate, 2)))
        bits_per_slice = int(math.ceil(
            (capacity * abs(math.log(error_rate))) /
            (num_slices * (math.log(2) ** 2))))
        self._setup(error_rate, num_slices, bits_per_slice, capacity)
        self.buckets = [None for x in range(self.num_bits)]
        self.exclusions = []
        self.childLayer = None
        self.depth = depth
        # set the "salt" for this layer. Note: Keeler makes the point that so
        # long as each layer size is relatively prime, we shouldn't get collisions.
        # If this is the case, we can be smarter about layer bitarray sizes and maybe
        # do away with per-layer salting.
        self.salt = "a" * depth

    def _setup(self, error_rate, num_slices, bits_per_slice, capacity):
        self.error_rate = error_rate
        self.num_slices = num_slices
        self.bits_per_slice = bits_per_slice
        self.capacity = capacity
        self.num_bits = num_slices * bits_per_slice
        self.make_hashes, self.hashfn = make_hashfuncs(self.num_slices, self.bits_per_slice)

    def get_layer_true_negatives(self):
        excludes = set()
        for bucket in self.buckets:
            if None != bucket:
                excludes.update(bucket.get_true_negatives())
        return excludes

    def get_layer_excludes(self):
        excludes = set()
        for bucket in self.buckets:
            if None != bucket:
                excludes.update(bucket.excludes)
        return excludes

    def make_layer_bit_array(self):
        ba = bitarray.bitarray(self.num_bits, endian='little')
        ba.setall(False)
        for i in range(self.num_bits):
            bucket = self.buckets[i]
            if None != bucket:
                if len(bucket.includes) > 0:
                    ba[i] = True
        return ba

    """Add a batch of entries to the filter cascade."""
    def add_items(self, inclusions, exclusions, propagate = False):
        # loop over the elements that should *not* be there. Create a new layer
        # that *includes* the false positives and *excludes* the true positives
        for elem in exclusions:
            self.record_for_layer(elem, False)

        possible_false_positives = set()

        # loop over the elements that should be there. Add them to the filter.
        for elem in inclusions:
            possible_false_positives.update(self.record_for_layer(elem, True))

        possible_false_positives.update(exclusions)

        if propagate:
            false_positives = possible_false_positives.difference(self.get_layer_true_negatives())
            if len(false_positives) > 0 or None != self.childLayer:
                entries = set()
                entries.update(inclusions)
                if None == self.childLayer:
                    self.childLayer = CascadeLayer(
                                # TODO: MDG - the size of the child layer should be calculated in
                                # a smarter way than this.
                                len(false_positives),
                                # 0.5, # <- suggested p for subsequent layers
                                self.error_rate,
                                self.depth + 1
                              )
                self.childLayer.add_items(false_positives, entries, True)
        # Propagation works like this:
        # Add the included entries, there will be new possible false positives
        # to check on account of the fact that some of these will hit buckets
        # with no previous excludes. Stick these in the set of possible
        # false positives along with the new excludes for this update. Check the
        # possible false positives to get a set of *actual* false positives.
        # Do the add_items thing for the next layer, with these new false positives
        # as the *includes* and any new additions as the *excludes*.

    """Remove a batch of entries from the filter cascade."""
    def remove_items(self, items):
        for bucket in self.buckets:
            if None != bucket:
                bucket.includes.difference_update(items)
                bucket.excludes.difference_update(items)
        if None != self.childLayer:
            self.childLayer.remove_items(items)

    def __contains__(self, elem):
        if self.filter_contains(elem):
            if self.childLayer is None:
                return True
            else:
                return not elem in self.childLayer

    def filter_contains(self, key):
        """Tests a key's membership in this bloom filter.
        """
        bits_per_slice = self.bits_per_slice
        hashes = self.make_layer_hashes(key)
        offset = 0
        for k in hashes:
            if len(self.buckets[offset + k].includes) == 0:
                return False
            offset += bits_per_slice
        return True

    def make_layer_hashes(self, key):
        if isinstance(key, str):
            return self.make_hashes(self.salt + key)
        return self.make_hashes(self.salt + str(key))

    """returns affected *excludes* for new include items, affected *includes* for
    new exclude items"""
    def record_for_layer(self, key, add):
        """Adds a record to this filter cascade."""
        bits_per_slice = self.bits_per_slice
        hashes = self.make_layer_hashes(key)
        offset = 0
        affected_items = set()
        for k in hashes:
            bucket = self.buckets[offset + k]
            if None == bucket:
                bucket = Bucket()
                self.buckets[offset + k] = bucket
            offset += bits_per_slice
            if add:
                affected_items.update(bucket.include(key))
            else:
                affected_items.update(bucket.exclude(key))
        return affected_items

    def bit_count(self):
        if self.childLayer is None:
            return len(self.make_layer_bit_array())
        return len(self.make_layer_bit_array()) + self.childLayer.bit_count()

    def layer_count(self):
        if self.childLayer is None:
            return 1
        else:
            return self.childLayer.layer_count() + 1