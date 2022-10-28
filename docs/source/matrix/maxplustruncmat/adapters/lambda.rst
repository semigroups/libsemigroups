.. Copyright (c) 2021, Finn Smith

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

Lambda
======

This page contains details of the specialisations of the
:ref:`adapters<Adapters>` :cpp:any:`Lambda` and :cpp:any:`LambdaValue` for
truncated max-plus matrices. 

.. cpp:struct:: template <typename Mat>  \
                LambdaValue<Mat>

   Specialization of the adapter :cpp:any:`LambdaValue` for instances of
   :cpp:any:`StaticMaxPlusTruncMat`.
 
   :tparam Mat: 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>` for some threshold ``T``, row
      dimension ``R``, and column dimension ``C``.

   .. cpp:type:: type = detail::StaticVector1<StaticMaxPlusTruncMat<T, R, \
                  C>::Row, R>

      For :cpp:any:`StaticMaxPlusTruncMat\<T, N>` objects, ``type`` is
      :cpp:any:`StaticVector1\<StaticMaxPlusTruncMat\<T, R, C>::Row, R>`. This
      represents a row space basis of a 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>`.


.. cpp:struct::  template <typename Mat> \
                 Lambda<Mat, typename LambdaValue<Mat>>
   
   Specialization of the adapter :cpp:any:`Lambda` for instances of
   :cpp:any:`StaticMaxPlusTruncMat`.
      
   :tparam Mat: 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>` for some threshold ``T`` row
      dimension ``R``, and column dimension ``C``.

   .. cpp:function:: void operator()(typename LambdaValue<Mat>& res, Mat \
                                     const& x) const
      
      Modifies ``res`` to contain the row space basis of ``x``.

      :param res: the container to hold the result
      :param x:   the matrix

      :returns: (None).

      :complexity: 
         :math:`O(R ^ 2 C)`.
