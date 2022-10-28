.. Copyright (c) 2020, Finn Smith

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

Konieczny
=========

``libsemigroups`` contains a preliminary implementation of Konieczny's
generalisation of the Lallement-McFadden algorithm for computing the structure
of finite semigroups that act on sets; see the documentation below for further
details. It is preliminary in the sense that the interface exposed by the
relevant classes is minimal, and that certain optimisations remain to be
performed.

The implementation of this algorithm is generic, and can be easily adapted
for user-defined types. In other words, it is possible to directly apply this
algorithm to any types that satisfy some minimal prerequisites, via traits
classes and/or the adapters described here.

The implementation of the Konieczny and Lallement-McFadden algorithm can be
used "interactively", in the sense that they can be run for a particular amount
of time, or until some condition is met; for further details see, for example,
:cpp:any:`void libsemigroups::Runner::run_for(std::chrono::nanoseconds)` and
:cpp:func:`libsemigroups::Runner::run_until`.

The classes implementing the Konieczny and Lallement-McFadden algorithm are:

.. toctree::
   :maxdepth: 1
  
   _generated/libsemigroups__konieczny
   _generated/libsemigroups__koniecznytraits
   _generated/libsemigroups__konieczny__dclass
