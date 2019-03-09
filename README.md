# ibim

This is the beginnings of a fairly general-purpose implementation of Walter Burkhard's Interpolation-based Index
Management hash technique. It builds on his work with linear hashing, which requires a hash function that maps
(as uniformly as possible) a key-space into the floating point [0.0 1.0), with the requirement that for any keys
K1 and K2 and hash function h, K1 < K2 => h(K1) < h(K2). This is clearly not as generally applicable as some
hash schemes where the mapping function is less constrained and can be designed to introduce a great deal of
entropy to help make more even use of the hash-table elements. With linear hashing, some care must be taken when
choosing a hash function to meet the above constraints while also trying to achieve uniform occupation of the
hash buckets. So, whereas a typical hash implementation may have a default hash function for strings which may, for
example, columnate the characters into an integral type of the size suitable to index into the hash table and then
continuously overlay the bit patters of the ASCII character values, xoring over the bytes below. This works well in
general, because minor difference in text are likely to produce vastly, and pseudo-randomly different hash function
outputs - h("john") will be spacially largely unrelated to h("joan"). With linear hashing, this can require a little
more conscious understanding of ones input data space. Imagince a key set of case-insensitive, alphabetic first names,
that happen to be very heavy with j-names, such "john", "joe", "joan", "joel", etc. A common linear-hashing-friendly
function that meets the basic requirement is to interpret the string as a base-26 fraction. This maps the entire
possible key space into [0.0 1.0) and has the desired quality that if a name N1 lexigraphically sorts before name N2,
then h(N1) < h(N2). But because the first character will be so dominant in determining the hash value, and because
we have noted that our key set is skewed toward "j" names, if we use this basic technique, our table will be over
utilized in a region around h(N) ~= 0.34 ("j" = 9 over 26), leading to either run-time behavior
approximating O(n) for a fixed table or wasted memory in an expanding table. Thus, if there is reason to believe
that the key set has statistical non-uniformities in it, the designer of the table may want to apply differing weights
to each key element or otherwise apply smoothing to the hash function. That's the bad news - linear hashing requires
more care to use as a general-purpose data-structure than some other approaches. The good news is that linear hashing
makes it easy to do range-query look-ups because of the order-preserving nature of the hash function. Bowles then extends
the idea to allow for multi-dimensional key-spaces with a generalized way of using a single table to store data with
multiple independent keys, along with techniques for identifying sets of hash buckets corresponding multi-dimensional
key regions. It's not quite as trivial as creating hyper-cubes of floating-point key-spaces, and it is best
described in Bowles' original paper on the subject, which I will try to find a reference to and place here in this
readme (though copyright concerns may forbid).

This implementation is also a playground for sharpening my variadic template skills, since that is how the data
structure is described. To begin with, it is in an XCode project, but the code is generic C++14 (or so) - whatever
language level Apple Clang supports in Xcode 10+. I think it may be C++11 compatible, the question being around some
of the variadic template support functions involving std::index_sequence - I think that only comes in at the C++14. If
this project has any broad appeal, though, it will probably just be as a starting-point for some implementation ideas
around IBIM or as an example (for good or ill) of variadic template use.
