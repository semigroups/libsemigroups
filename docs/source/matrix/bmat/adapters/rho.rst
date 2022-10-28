.. Copyleft (c) 2020, J. D. Mitchell

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

Rho
===

This page contains details of the specialisations of the
:ref:`adapters<Adapters>` :cpp:any:`Rho` and :cpp:any:`RhoValue` for boolean
matrices. 

.. cpp:struct:: template <typename Mat>  \
                RhoValue<Mat>

   Specialization of the adapter ``RhoValue`` for instances of ``BMat<N>``.
   Note that the type chosen here limits the :cpp:any:`Konieczny` algorithm to
   boolean matrices of dimension at most 64 (or 32 on 32-bit systems).
 
   :tparam Mat: 
     a type such that :cpp:any:`IsBMat<Mat>` is ``true``. 

   .. cpp:var:: static constexpr size_t N = BitSet<1>::max_size()
      
      Variable representing the maximum capacity of the container for lambda
      values. 
      
   .. cpp:type:: type = detail::StaticVector1<BitSet<N>, N>

      For :cpp:any:`BMat\<N>` objects, ``type`` is
      :cpp:any:`StaticVector1\<BitSet\<N>, N>`, where ``N`` is the maximum
      width of BitSet on the system. This represents a row space basis of a
      :cpp:any:`BMat\<N>`.
  

.. cpp:struct::  template <typename Container, typename Mat> \
                 Rho<Mat, Container>
      
   :tparam Container:
     a container type for holding the return value of the call operator, can be
     one of 
     :cpp:any:`StaticVector1\<BitSet\<N>, N>`;
     :cpp:any:`std::vector\<BitSet\<N>>`;
     :cpp:any:`StaticVector1\<std::bitset\<N>, N>`; or 
     :cpp:any:`std::vector\<std::bitset\<N>>`
     for some value of ``N``.

   :tparam Mat: 
     a type such that :cpp:any:`IsBMat<Mat>` is ``true``. 

   .. cpp:function:: void operator()(Container& res, Mat const& x) const
      
      Modifies ``res`` to contain the row space basis of ``x``.

      :param res: the container to hold the result
      :param x:   the matrix

      :returns: (None).

      :complexity: 
         :math:`O(r ^ 2 c)` where :math:`r` is the number of rows in ``x`` and
         :math:`c` is the number of columns in ``x``. 
