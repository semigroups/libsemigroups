.. Copyright (c) 2020, J. D. Mitchell

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

Matrix
======

There are several different classes in ``libsemigroups`` for representing
matrices over various semirings. There are up to three different
representations for every type of matrix:

1. those whose dimension and arithmetic operations can be defined at compile
   time: :cpp:any:`StaticMatrix`.
2. those whose arithmetic operation can be defined at compile time but whose
   dimesions can be set at run time: :cpp:any:`DynamicMatrix` 
3. those whose arithmetic operation and dimensions can be set at run time:
   :cpp:any:`DynamicMatrix`

It's unlikely that you will want to use the classes in described on this page
directly, but rather through the aliases described on the other matrix pages
(such as, for example, :cpp:any:`BMat`). 

Over specific semirings
-----------------------

The following matrix classes are provided which define matrices over some
specific semirings:

.. toctree::
   :maxdepth: 1

   bmat/index
   intmat
   maxplusmat
   minplusmat
   maxplustruncmat/maxplustruncmat
   minplustruncmat
   ntpmat 
   projmaxplus 


Over arbitrary semirings
------------------------

The following general matrix classes are provided which can be used to define
matrices over arbitrary semirings:

.. toctree::
   :maxdepth: 1

   dynamicmatrix-compile
   dynamicmatrix-run
   staticmatrix
   variable-templates
   validation


Row views
---------

A row view is a lightweight representations of a row of matrix.  The following
row view classes are provided:

.. toctree::
   :maxdepth: 1

   dynamicrowview
   staticrowview


Helpers
-------

The namespace :cpp:any:`matrix_helpers` contains a number of helper functions
for certain types of matrices:

.. toctree::
   :maxdepth: 1

   helpers/pow
   helpers/rows
   helpers/row_basis
   helpers/row_space_size
   helpers/bitset_rows
   helpers/bitset_row_basis

Adapters
---------------

There are various specialisations of the adapters described on :ref:`Adapters`
for the matrices desribed on this page:

.. toctree::
   :maxdepth: 1

   adapters/complexity
   adapters/degree
   adapters/increasedegreeby
   adapters/one
   adapters/product
