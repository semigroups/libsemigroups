.. Copyright (c) 2019, J. D. Mitchell

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

.. _Adapters:

Adapters
========

Declared in ``adapters.hpp``.

This file contains class templates for adapting a user-defined type for use
with the algorithms in ``libsemigroups``. There are explicit implementations for
relevant element types from ``libsemigroups`` and ``HPCombi``, details of which
can be found elsewhere in this documentation.

No default implementation
~~~~~~~~~~~~~~~~~~~~~~~~~

The following adapters have no default implementation for arbitrary types:

.. toctree::
   :maxdepth: 1

   _generated/libsemigroups__complexity.rst
   _generated/libsemigroups__degree.rst
   _generated/libsemigroups__imageleftaction.rst
   _generated/libsemigroups__imagerightaction.rst
   _generated/libsemigroups__increasedegree.rst
   _generated/libsemigroups__inverse.rst
   _generated/libsemigroups__lambda.rst
   _generated/libsemigroups__lambdavalue.rst
   _generated/libsemigroups__one.rst
   _generated/libsemigroups__product.rst
   _generated/libsemigroups__rank.rst
   _generated/libsemigroups__rho.rst
   _generated/libsemigroups__rhovalue.rst

Default implementation
~~~~~~~~~~~~~~~~~~~~~~

The following adapters do have a default implementation for arbitrary types:

.. The following are not in the yml setup because JDM couldn't figure out how
   to generate the rst page for EqualTo::operator().

.. toctree::
   :maxdepth: 1

   api/equalto.rst
   api/hash.rst
   api/less.rst
   api/ontuples.rst
   api/onsets.rst
   _generated/libsemigroups__rankstate.rst
   api/swap.rst
