.. Copyright (c) 2021, Finn Smith

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

Rank
====

This page contains details of the specialisations of the :ref:`adapter<Adapters>`
:cpp:any:`Rank` for truncated max-plus matrices. 


.. cpp:struct:: template <typename Mat> \
                Rank<Mat, RankState<Mat>>

   Specialization of the adapter :cpp:any:`Rank` for instances of
   :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>`.

   :tparam Mat: 
      :cpp:any:`StaticMaxPlusTruncMat\<T, R, C>` for some threshold ``T``, row
      dimension ``R``, and column dimension ``C``.
  
   .. cpp:function:: size_t operator()(RankState<Mat> const& state, Mat const& x) const

      Returns the rank of ``x``.

      :param state: a const reference to the `RankState<Mat>`

      :param x: the matrix

      :returns: 
        the rank of the parameter ``x``, which is a value of type ``size_t``.

      :complexity:
        :math:`O(TRCn)` where :math:`n` is the rank of ``x``, i.e. the size of
        the row space of ``x``.
