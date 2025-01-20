Information for developers
==========================

Building the documentation from source
--------------------------------------

If you are checking out the libsemigroups git repo, then the following steps
are required to build the documentation:

1. Install doxygen (using homebrew or apt-get)
2. Check out the doxygen-awesome-css/ submodule using:

   git submodule update --init --recursive

3. Type `make doc`

See the `manual <https://libsemigroups.readthedocs.io/en/latest/install.html>`_
for details of how to do this.

Naming conventions
------------------

====================  ============================
Object                Convention
--------------------  ----------------------------
class                 NamedLikeThis
struct                NamedLikeThis
function              named\_like\_this
enum class            namedlikethis
enum class member     named_like_this
stateful type alias   named\_like\_this\_type
stateless type alias  LikeThis
iterator              iterator or const\_iterator
global variable       NAMED\_LIKE\_THIS
====================  ============================

Try to ensure that any new functionality added to ``libsemigroups`` satisfies
the strong exception guarantee, and if not that this is documented.

Some more conventions:

* If a function does not check its arguments (including `*this` for member
  functions), then it must have the suffix `_no_checks` if there is any
  possibility that the function will fail. This indicates that anything may
  happen if the function is called on a non-valid object or for non-valid
  arguments.

  In the python bindings ``libsemigroups_pybind11`` it should not be possible
  to construct or modify an object so that it is invalid. For example, it
  shouldn't be possible to create an invalid `Bipartition` or to modify a valid
  `Bipartition` so that it is invalid. So, it is safe to call the
  `Bipartition::rank_no_checks` member function as if it was the
  `Bipartition::rank` function. Hence in the python bindings we bind a
  functions called `Bipartition.rank` to the C++ function
  `Bipartition::rank_no_checks`.

* All class helper functions (i.e. those free functions taking a class as an
  argument, and using public member functions of that class) should be in a
  helper namespace with the same name as the class (but in lower case). I.e.
  `bipartition::one`, or `ptransf::one`. Exceptions are:

  - operators such as `operator*` etc. This is because these are more or less
    difficult to use if they are in a nested namespace.
  - functions that are (or should be) implemented for all of the public facing
    classes in ``libsemigroups``, such as, for example
    `to_human_readable_repr`

* Functions that construct an object must be of one of three types;
  constructors, `make` functions, or `to` functions:

  * A constructor should be used to construct an object without checking its
    arguments.
  * A `make` function should be a free function template in the libsemigroups
    namespace that is used to construct an object from non-libsemigroups
    objects such as containers, integers and strings. It should check its
    arguments.
  * A `to` function should a free function template in the libsemigroups
    namespace that is used to convert from one `libsemigroups` type to another.
    A typical signature might look something like
    `to<ToddCoxeter>(Presentation<word_type> p)`. It should check its
    arguments.

* The documentation for each function may contain the following section
  indicators in the following order:
  #. ``\tparam``
  #. ``\param``
  #. ``\return``
  #. ``\exceptions`` or ``\throws``
  #. ``\complexity``
  #. ``\note``
  #. ``\warning``
  #. ``\sa``
  This can be checked by running ``etc/check_doc_order.py``.

Debugging and valgrinding
-------------------------

To run `lldb`:

.. code-block:: bash

    ./configure --enable-debug && make test_all
    libtool --mode=execute lldb test_all

`test_all` is the name of the check program produced by `make check`. Similarly
to run valgrind you have to do:

.. code-block:: bash

    ./configure --enable-debug --disable-hpcombi && make test_all
    libtool --mode=execute valgrind --leak-check=full ./test_all [quick] 2>&1 | tee --append valgrind.txt

Adding new test cases
---------------------

Any new tests should be tagged with one of the following:

========  =======
Tag       Runtime
--------  -------
quick     < 200ms
standard  < 3s
extreme   > 3s
========  =======

They should be declared using

.. code-block:: cpp

    LIBSEMIGROUPS_TEST_CASE(classname, number, message, tags)

"tags" should include "[FilePrefix]" where "FilePrefix" would be
"cong-pair" in the file "tests/cong-pair.test.cc", if the file prefix is the
same as "classname", then it should not be included. Tags are case insensitive.

Making a release
----------------

A ***bugfix release*** is one of the form `x.y.z -> x.y.z+1`, and
a ***non-bugfix release*** is one of the form `x.y.z -> x+1.y.z` or `x.y+1.z`.

Use the script `etc/release-libsemigroups.py`.

References
----------

- <https://conda.io/docs/build_tutorials.html>
- <https://github.com/conda/conda-recipes>
