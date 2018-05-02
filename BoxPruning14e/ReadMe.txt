Box pruning example
-------------------

This is the simplest box pruning implementation I've ever coded,
yet it sometimes performs surprisingly well compared to more clever
schemes. See for yourself.

"Bipartite" pruning finds intersections of one set of boxes with
another set of boxes.

"Complete" pruning finds intersections within a single set of boxes.

The nice thing is that it doesn't rely on coherence to do well,
doesn't use extra precomputed data structures, doesn't use extra
ram. First shot is roughly as fast as next ones, which is nice for
one-shot queries.


Pierre Terdiman
February, 10, 2002

www.codercorner.com