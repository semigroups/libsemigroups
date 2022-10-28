.. Copyright (c) 2020, J. D. Mitchell

   Distributed under the terms of the GPL license version 3.

   The full license is in the file LICENSE, distributed with this software.

.. cpp:namespace:: libsemigroups

Truncated max-plus matrices
===========================

Defined in ``matrix.hpp``.

This page describes the functionality for :math:`n \times n`  matrices over the
finite quotient of the max-plus semiring by the congruence :math:`t = t + 1`
for arbitrary :math:`n` and :math:`t`. The value :math:`t` is referred to as
the *threshold*. 

There are three types of such matrices where:

1. the dimension is known at compile-time;
2. the dimension is to be defined a run time but the arithmetic operations
   are known at compile-time (i.e. the value of :math:`t` is known at compile
   time)
3. both the dimension and the arithmetic operations (i.e. :math:`t`) are to be 
   defined a run time.

All three of these types can be accessed via the alias template
:cpp:any:`MaxPlusTruncMat\<T, R, C, Scalar>`: if ``T`` has value ``0``, then the
threshold can be set at run time, and if ``R`` or ``C`` is ``0``, then the
dimension can be set at run time.  The default value of ``T`` is ``0``, ``R``
is ``0``, and of ``C`` is ``R``. 

The alias :cpp:any:`MaxPlusTruncMat\<T, R, C, Scalar>` is either
:cpp:any:`StaticMatrix` or :cpp:any:`DynamicMatrix`, please refer to the
documentation of these class templates for more details. The only substantial
difference in the interface of :cpp:any:`StaticMatrix` and
:cpp:any:`DynamicMatrix` is that the former can be default constructed and the
latter should be constructed using the dimensions. 

**Example**

.. code-block:: cpp

   MaxPlusTruncMat<11, 3> m;       // default construct an uninitialized 3 x 3 static matrix with threshold 11
   MaxPlusTruncMat<11> m(4, 4);    // construct an uninitialized 4 x 4 dynamic matrix with threshold 11
   MaxPlusTruncSemiring sr(11);    // construct a truncated max-plus semiring with threshold 11 
   MaxPlusTruncMat<>  m(sr, 5, 5); // construct an uninitialized 5 x 5 dynamic matrix with threshold 11 (defined at run time)


.. cpp:struct:: template <size_t T, typename Scalar> \
                MaxPlusTruncProd

    This is a stateless struct with a single call operator of signature:
    ``Scalar operator()(Scalar const x, Scalar const y) const noexcept``
    that returns :math:`x \otimes y` which is defined by 

    .. math::

       x\otimes y = 
       \begin{cases}
       \min\{x + y, T\}   & \text{if } x \neq -\infty\text{ and }y \neq -\infty \\
       -\infty & \text{if } x = -\infty \text{ or }y = -\infty; \\
       \end{cases}

    representing multiplication in the quotient of the max-plus semiring by the
    congruence :math:`T = T + 1`.


.. cpp:class:: template <typename Scalar = int> \
               MaxPlusTruncSemiring final 
 
   This class represents the **max-plus truncated semiring** consists of the
   integers :math:`\{0, \ldots , t\}` for some value :math:`t` (called the
   **threshold** of the semiring) and :math:`-\infty`. Instances of this class
   can be used to define the value of the threshold :math:`t` at run time. 

   :tparam Scalar: 
      the type of the elements of the semiring. This must be a signed integer
      type. 

   .. cpp:function:: MaxPlusTruncSemiring() = delete

      Deleted default constructor. 
   
   .. cpp:function:: MaxPlusTruncSemiring(MaxPlusTruncSemiring const&) = default

      Default copy constructor.
   .. cpp:function:: MaxPlusTruncSemiring(MaxPlusTruncSemiring&&) = default

      Default move constructor.

   .. cpp:function:: MaxPlusTruncSemiring& operator=(MaxPlusTruncSemiring const&) \
                 = default

      Default copy assignment operator.

   .. cpp:function:: MaxPlusTruncSemiring& operator=(MaxPlusTruncSemiring&&) \
                 = default
      
      Default move assignment operator.
    
   .. cpp:function:: explicit MaxPlusTruncSemiring(Scalar const threshold) 
      
      Construct from threshold.

      :param threshold: the threshold.
     
      :throws:
        :cpp:any:`LibsemigroupsException` if ``threshold`` is less than zero.
      
      :complexity: Constant.
    
   .. cpp:function:: Scalar zero() const noexcept 
      
      Returns :math:`-\infty`; representing the additive identity of the
      quotient of the max-plus semiring.
      
      :parameters: (None)

      :returns: A value of type ``Scalar``.

      :exceptions: 
        This function is ``noexcept`` and is guaranteed never to throw.
      
      :complexity: Constant.
   
   .. cpp:function:: Scalar one() const noexcept 
       
      Returns :math:`0`; representing the multiplicative identity of the
      quotient of the max-plus semiring.

      :parameters: (None)

      :returns: A value of type ``Scalar``.

      :exceptions: 
        This function is ``noexcept`` and is guaranteed never to throw.
      
      :complexity: Constant.

   .. cpp:function:: Scalar plus(Scalar const x, Scalar const y) const noexcept
       
      Returns :math:`x \oplus y` which is defined by 

      .. math::

        x\oplus y = 
        \begin{cases}
        \max\{x, y\}   & \text{if } x \neq -\infty\text{ and }y \neq -\infty \\
        -\infty & \text{if } x = -\infty \text{ or }y = -\infty; \\
        \end{cases}

     representing addition in the max-plus semiring (and its quotient).
       
     :param x: scalar
     :param y: scalar

     :returns: A value of type ``Scalar``.

     :exceptions: 
       This function is ``noexcept`` and is guaranteed never to throw.
     
     :complexity: Constant.

   .. cpp:function:: Scalar prod(Scalar const x, Scalar const y) const noexcept
       
      Returns :math:`x \otimes y` which is defined by 

     .. math::

        x\otimes y = 
        \begin{cases}
        \min\{x + y, t\}   & \text{if } x \neq -\infty\text{ and }y \neq -\infty \\
        -\infty & \text{if } x = -\infty \text{ or }y = -\infty; \\
        \end{cases}

     where :math:`t` is the threshold; representing multiplication in the
     quotient of the max-plus semiring.
     
     :param x: scalar
     :param y: scalar
     
     :returns: A value of type ``Scalar``.

     :exceptions: 
       This function is ``noexcept`` and is guaranteed never to throw.
     
     :complexity: Constant.
     
   .. cpp:function:: Scalar threshold() const noexcept
      
      Returns the threshold value used to construct ``this``. 

      :parameters: (None)

      :returns: A value of type ``Scalar``.

      :exceptions: 
        This function is ``noexcept`` and is guaranteed never to throw.
      
      :complexity: Constant.


.. cpp:type:: template <typename Scalar> \
              DynamicMaxPlusTruncMatSR = \
              DynamicMatrix<MaxPlusTruncSemiring<Scalar>, Scalar>
    
   Alias for truncated max-plus matrices with dimensions and threshold defined
   at runtime. 

   :tparam Scalar: 
     The type of the entries in the matrix. 


.. cpp:type:: template <typename T, typename Scalar>       \
              DynamicMaxPlusTruncMat                       \
              = DynamicMatrix<MaxPlusPlus<Scalar>,         \
                              MaxPlusTruncProd<T, Scalar>, \
                              MaxPlusZero<Scalar>,         \
                              IntegerZero<Scalar>,         \ 
                              Scalar>

   Alias for the type of dynamic max-plus matrices where the dimension
   is defined at run time, but the threshold is defined at compile-time.
   
   :tparam T: 
     the threshold. 
   
   :tparam Scalar: 
     the type of the entries in the matrix. 


.. cpp:type::  template <size_t T,                                               \
                         size_t R,                                               \
                         size_t C,                                               \
                         typename Scalar>                                        \
               StaticMaxPlusTruncMat = StaticMatrix<MaxPlusPlus<Scalar>,         \
                                                    MaxPlusTruncProd<T, Scalar>, \
                                                    MaxPlusZero<Scalar>,         \
                                                    IntegerZero<Scalar>,         \
                                                    R,                           \
                                                    C,                           \
                                                    Scalar>

   Alias for static max-plus truncated matrices where the threshold and
   dimensions are defined at compile-time. 

   :tparam T: the threshold. 

   :tparam R: the number of rows.  

   :tparam C: the number of columns.  

   :tparam Scalar: 
     The type of the entries in the matrix. 


.. cpp:type:: template <size_t T = 0, \
                        size_t R = 0, \
                        size_t C = R, \ 
                        typename Scalar = int> \
              MaxPlusTruncMat = std::conditional_t< \
                                  R == 0 || C == 0, \
                                  std::conditional_t<T == 0, \
                                    DynamicMaxPlusTruncMatSR<Scalar>,  \
                                    DynamicMaxPlusTruncMat<T, Scalar>>, \
                                  StaticMaxPlusTruncMat<T, R, C, Scalar>>
   
   :tparam T: 
     the threshold. A value of ``0`` indicates that the value will be set at
     run time (default: ``0``).

   :tparam R: 
     the number of rows.  A value of ``0`` indicates that the value will be set
     at run time (default: ``0``).

   :tparam C: 
     the number of columns.  A value of ``0`` indicates that the value will be
     set at run time (default: ``R``).

   :tparam Scalar: 
     The type of the entries in the matrix (default: ``int``). 


.. cpp:var:: template <typename T> \
             static constexpr bool IsMaxPlusTruncMat

   This variable has value ``true`` if the template parameter ``T`` is the same
   as :cpp:any:`MaxPlusTruncMat\<T, R, C, Scalar>` for some values of ``T``,
   ``R``, ``C``, and ``Scalar``; and ``false`` if it is not. 

Truncated max-plus matrix adapters
----------------------------------

The following adapters are implemented to so that
:cpp:any:`StaticMaxPlusTruncMat\<T, R, C, Scalar>` can be used with
:cpp:any:`Konieczny`. 

.. toctree::
   :maxdepth: 1
   
   adapters/imagerightaction
   adapters/imageleftaction
   adapters/lambda
   adapters/rho
   adapters/rank
