# live_cascade
A sketch (in Python) of a technique for live updates of multi-level bloom filters

# An explanation
(from an explanation I sent a co-worker)

I think I have a way to significantly reduce the cost of calculating filter updates. Have a read and let me know what you think.

Consider the set of certificates, C. There are two subsets; A and B. A represents the sets of revoked certs. B, the set of non-revoked certs.

Let's start off by thinking about how a filter cascade represents these data. The cascade consists of layers. Each layer is a bit array, each position representing a bit that is set in a hash of some element of a set. Which set the element belongs to depends on the layer we're looking at (0th and even layers; set A, otherwise set B).

Imagine, instead of a bit array of m elements, we have m buckets. Each bucket contains two lists. The first contains data pertaining elements of set A which have that bit set that, the other, data pertaining to elements of B.

To build the filter, we take each element of the set C. We calculate all hashes needed for a filter cascade of depth d.

Next, we populate the lists in each bucket in the first layer with these sets of d hashes. List X gets the data for elements of A for which any of the 0th set of hashes have this bit set. List Y gets the data for elements of B.

(I'm being lazy from this point, A and B now refer to sets of hashes for the elements of the original A and B)

We can then iterate over the buckets and remove elements from B which exist in buckets with empty X lists. Once this is complete, we can start again on the next layer (incrementing the index of the set of hashes). For each layer, we reverse the sources of X and Y. We repeat until A and B are empty.

We can use these buckets produce the bit arrays for each layer. To do this, we iterate over the buckets and set only bits with non-empty X lists.

For updates, we add new sets of A and B. If an entry results in an addition to a previously empty X list, we must consider the Y elements when processing the next layer.

We only need to recompute from scratch if d changes (or we need to resize).

I'm going to sketch this out in code later on.
