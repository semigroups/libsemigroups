.. Copyright (c) 2021, Finn Smith

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.


ImageRightAction
================

.. cpp:namespace:: libsemigroups

.. cpp:struct:: template <typename Mat> \
                ImageRightAction<Mat, detail::StaticVector1<Mat::Row, \
                Mat::nr_rows>>
  
   This struct is a specialisation of :cpp:any:`ImageRightAction` for square
   truncated max-plus matrices and :cpp:any:`detail::StaticVector1` of rows. 

   :tparam Mat: 
      :cpp:any:`StaticMaxPlusTruncMat\<T, N>` for some threshold ``T`` and
      dimension ``N``.

   .. cpp:type:: result_type

      This is an alias for :cpp:any:`detail::StaticVector1\<Mat::Row, Mat::nr_rows>`.

   .. cpp:function:: void operator()(result_type& res,      \
                                     result_type const& pt, \ 
                                     Mat const& x) const
      
      Given a row basis ``pt`` and a square truncated max-plus matrix ``x``
      such that ``x`` has the same number of rows as each element of ``pt``,
      this function computes the row basis of the space obtained by right
      multiplication of the rows in ``pt`` by the matrix ``x``, and modifies
      ``res`` to contain the resulting row basis.

      :param res: container to store the result.

      :param pt:  the container on which to act.

      :param x:   the matrix which is acting. 

      :returns: (None). 

      :complexity: 
         :math:`O(a ^ 2 r)` where :math:`a` is the size of ``pt`` and :math:`r`
         is the number of rows in ``x``. 
