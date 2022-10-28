.. Copyright (c) 2021, Finn Smith

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

Rho
===

This page contains details of the specialisations of the
:ref:`adapters<Adapters>` :cpp:any:`Rho` and :cpp:any:`RhoValue` for truncated
max-plus matrices. 

.. cpp:struct:: template <typename Mat>  \
                RhoValue<Mat>

   Specialization of the adapter :cpp:any:`RhoValue` for instances of
   :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>`.
 
   :tparam Mat: 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>` for some threshold ``T``, row
      dimension ``R``, and column dimension ``C``.

   .. cpp:type:: type = detail::StaticVector1<StaticMaxPlusTruncMat<T, C, \
                  R>::Row, C>

      For :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>` objects, ``type`` is
      :cpp:any:`StaticVector1\<StaticMaxPlusTruncMat\<T, C, R>::Row, C>`. This
      represents a column space basis of a 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>`.


.. cpp:struct::  template <typename Mat> \
                 Rho<Mat, typename RhoValue<Mat>>

   Specialization of the adapter :cpp:any:`Rho` for instances of
   :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>`.
      
   :tparam Mat: 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>` for some threshold ``T`` row
      dimension ``R``, and column dimension ``C``.

   .. cpp:function:: void operator()(typename RhoValue<Mat>& res, Mat \
                                     const& x) const
      
      Modifies ``res`` to contain a column space basis of ``x``.

      :param res: the container to hold the result
      :param x:   the matrix

      :returns: (None).

      :complexity: 
         :math:`O(C ^ 2 R)`.
