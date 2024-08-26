//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-2024 Finn Smith
// Copyright (C) 2019-2024 James D. Mitchell
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

// This file contains declarations related to word graphs (which are basically
// deterministic automata without initial or accept states).

#ifndef LIBSEMIGROUPS_WORD_GRAPH_HPP_
#define LIBSEMIGROUPS_WORD_GRAPH_HPP_

#include <algorithm>  // for uniform_int_distribution
#include <cstddef>    // for size_t
#include <cstdint>
#include <ostream>        // for operator<<
#include <queue>          // for queue
#include <random>         // for mt19937
#include <stack>          // for stack
#include <string>         // for to_string
#include <tuple>          // for tie
#include <type_traits>    // for is_integral, is_unsigned
#include <unordered_set>  // for unordered_set
#include <utility>        // for pair
#include <vector>         // for vector

#include "config.hpp"     // for LIBSEMIGROUPS_EIGEN_EN...
#include "constants.hpp"  // for UNDEFINED
#include "debug.hpp"      // for LIBSEMIGROUPS_ASSERT
#include "dot.hpp"        // for Dot
#include "exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "forest.hpp"     // for Forest
#include "order.hpp"      // for Order
#include "ranges.hpp"     // for ??
#include "types.hpp"      // for word_type

#include "detail/containers.hpp"  // for DynamicArray2
#include "detail/int-range.hpp"   // for IntRange
#include "detail/stl.hpp"         // for IsIterator
#include "detail/uf.hpp"          // for Duf

#ifdef LIBSEMIGROUPS_EIGEN_ENABLED
#include "detail/eigen.hpp"
#else
#include "matrix.hpp"
#endif

namespace libsemigroups {

  // TODO(doc)
  struct WordGraphBase {};

  //! \defgroup word_graph_group Word graphs and related functionality
  //!
  //! This page collects various classes and functions related to word graphs
  //! in ``libsemigroups``.

  //! \ingroup word_graph_group
  //!
  //! Defined in ``word-graph.hpp``.
  //!
  //! This class represents a word graph. If the word graph has \c n nodes,
  //! they are represented by the numbers \f$\{0, ..., n - 1\}\f$, and every
  //! node has the same number \c m of out-edges (edges with source that node
  //! and target any other node). The number \c m is referred to as the
  //! *out-degree* of the word graph, or any of its nodes.
  //!
  //! \tparam Node the type of the nodes in the word graph, must be an unsigned
  //! integer type.
  template <typename Node>
  class WordGraph : private WordGraphBase {
    static_assert(std::is_integral<Node>(),
                  "the template parameter Node must be an integral type!");
    static_assert(
        std::is_unsigned<Node>(),
        "the template parameter Node must be an unsigned integral type!");

    template <typename N>
    friend class WordGraph;

   public:
    ////////////////////////////////////////////////////////////////////////
    // WordGraph - typedefs - public
    ////////////////////////////////////////////////////////////////////////

    //! The type of nodes in a word graph.
    using node_type = Node;

    //! The type of edge labels in a word graph.
    using label_type = Node;

    //! Unsigned integer type.
    using size_type = std::size_t;

#ifdef LIBSEMIGROUPS_EIGEN_ENABLED
    //! The type of the adjacency matrix.
    using adjacency_matrix_type
        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
#else
    //! The type of the adjacency matrix.
    using adjacency_matrix_type = IntMat<0, 0, int64_t>;
#endif

    //! The type of an iterator pointing to the nodes of a word graph.
    using const_iterator_nodes =
        typename detail::IntRange<Node>::const_iterator;

    //! The type of a reverse iterator pointing to the nodes of a word graph.
    using const_reverse_iterator_nodes =
        typename detail::IntRange<Node>::const_reverse_iterator;

    //! The type of an iterator pointing to the out-edges of a node in a
    //! word graph.
    using const_iterator_targets =
        typename detail::DynamicArray2<Node>::const_iterator;

    ////////////////////////////////////////////////////////////////////////
    // WordGraph - constructors + destructor - public
    ////////////////////////////////////////////////////////////////////////

    //! Construct from number of nodes and out degree.
    //!
    //! \param m the number of nodes in the word graph (default: 0).
    //! \param n the out-degree of every node (default: 0).
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes, and \p n is
    //! the out-degree of the word graph.
    // Not noexcept
    explicit WordGraph(size_type m = 0, size_type n = 0);

    //! Re-initialize the word graph to have \p m nodes and out-degree \p n
    //!
    //! This functions puts a word graph into the state that it would have been
    //! in if it had just been newly constructed with the same parameters \p m
    //! and \p n.
    //!
    //! \param m the number of nodes in the word graph.
    //! \param n the out-degree of every node.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes, and \p n is the
    //! out-degree of the word graph.
    WordGraph& init(size_type m = 0, size_type n = 0);

    //! Default copy constructor
    WordGraph(WordGraph const&);

    // TODO(1) other versions from OtherNode
    template <typename OtherNode>
    WordGraph(WordGraph<OtherNode> const&);

    template <typename OtherNode>
    WordGraph& init(WordGraph<OtherNode> const&);

    //! Default move constructor
    WordGraph(WordGraph&&);

    //! Default copy assignment constructor
    WordGraph& operator=(WordGraph const&);

    //! Default move assignment constructor
    WordGraph& operator=(WordGraph&&);

    ~WordGraph();

    //! Construct a random word graph from number of nodes and out-degree.
    //!
    //! \param number_of_nodes the number of nodes
    //! \param out_degree the out-degree of every node
    //! \param mt a std::mt19937 used as a random source (defaults to:
    //! std::mt19937(std::random_device()()))
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p number_of_nodes is less than \c 2
    //! * \p out_degree is less than \c 2
    //!
    //! \par Complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes, and \p n is
    //! the out-degree of the word graph.
    static WordGraph random(size_type    number_of_nodes,
                            size_type    out_degree,
                            std::mt19937 mt
                            = std::mt19937(std::random_device()()));

    //! Construct a random word graph from number of nodes, edges, and
    //! out-degree.
    //!
    //! \param number_of_nodes the number of nodes
    //! \param out_degree the out-degree of every node
    //! \param number_of_edges ??
    //! \param mt a std::mt19937 used as a random source (defaults to:
    //! std::mt19937(std::random_device()()))
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p number_of_nodes is less than \c 2
    //! * \p out_degree is less than \c 2
    //! * \p number_of_edges exceeds the product of \p number_of_nodes and \p
    //! out_degree
    //!
    //! \par Complexity
    //! At least \f$O(mn)\f$ where \p m is the number of nodes, and \p n is the
    //! out-degree of the word graph.
    static WordGraph random(size_type    number_of_nodes,
                            size_type    out_degree,
                            size_type    number_of_edges,
                            std::mt19937 mt
                            = std::mt19937(std::random_device()()));

    //! Construct a random acyclic word graph from number of nodes, edges, and
    //! out-degree.
    //!
    //! \param number_of_nodes the number of nodes
    //! \param out_degree the out-degree of every node
    //! \param number_of_edges TODO
    //! \param mt a std::mt19937 used as a random source (defaults to:
    //! std::mt19937(std::random_device()()))
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p number_of_nodes is less than \c 2
    //! * \p out_degree is less than \c 2
    //! * \p number_of_edges exceeds the product of \p number_of_nodes and \p
    //! out_degree
    //! * \p number_of_edges exceeds the product of \p number_of_nodes and \p
    //! number_of_nodes - 1 divided by 2.
    //!
    //! \par Complexity
    //! At least \f$O(mn)\f$ where \p m is the number of nodes, and \p n is the
    //! out-degree of the word graph.
    static WordGraph random_acyclic(size_type    number_of_nodes,
                                    size_type    out_degree,
                                    size_type    number_of_edges,
                                    std::mt19937 mt
                                    = std::mt19937(std::random_device()()));

    //! Ensures that the word graph has capacity for a given number of nodes,
    //! and out-degree.
    //!
    //! \note
    //! Does not modify number_of_nodes() or out_degree().
    //!
    //! \param m the number of nodes
    //! \param n the out-degree
    //!
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes and \p n is the
    //! out-degree.
    //!
    //! \iterator_validity
    //! \iterator_invalid
    // Not noexcept because DynamicArray2::add_cols isn't.
    WordGraph& reserve(size_type m, size_type n);

    ////////////////////////////////////////////////////////////////////////
    // WordGraph - modifiers - public
    ////////////////////////////////////////////////////////////////////////

    //! Adds a number of new nodes.
    //!
    //! \param nr the number of nodes to add.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //! \strong_guarantee
    //!
    //! \complexity
    //! Linear in `number_of_nodes() + nr`.
    //!
    //! \iterator_validity
    //! \iterator_invalid
    // Not noexcept because DynamicArray2::add_rows isn't.
    WordGraph& add_nodes(size_type nr);

    //! Adds to the out-degree.
    //!
    //! \param nr the number of new out-edges for every node.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //! \strong_guarantee
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is the number of nodes, and \c n is the new out
    //! degree of the word graph.
    //!
    //! \iterator_validity
    //! \iterator_invalid
    // Not noexcept because DynamicArray2::add_cols isn't.
    WordGraph& add_to_out_degree(size_type nr);

    //! Add an edge from one node to another with a given label.
    //!
    //! If \p i and \p j are nodes in \c this, and \p lbl is in the range `[0,
    //! out_degree())`, then this function adds an edge from \p i to \p j
    //! labelled \p lbl.
    //!
    //! \param m the source node
    //! \param lbl the label of the edge from \p m to \p n
    //! \param n the range node
    //!
    //! \throws LibsemigroupsException if \p i, \p j, or \p lbl is
    //! not valid.
    //! \strong_guarantee
    //!
    //! \complexity
    //! Constant.
    // Not noexcept because throw_if_node_out_of_bounds/label aren't
    // return *this by reference.
    WordGraph& target(node_type m, label_type lbl, node_type n);

    //! Add an edge from one node to another with a given label.
    //!
    //! \param i the source node
    //! \param j the range node
    //! \param lbl the label of the edge from \p i to \p j
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    inline WordGraph& target_no_checks(node_type  m,
                                       label_type lbl,
                                       node_type  n) {
      _dynamic_array_2.set(m, lbl, n);
      return *this;
    }

    //! Remove an edge from a node with a given label.
    //!
    //! \param i the source node
    //! \param lbl the label of the edge from \p i
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    inline WordGraph& remove_target_no_checks(node_type s, label_type a) {
      _dynamic_array_2.set(s, a, UNDEFINED);
      return *this;
    }

    // TODO doc
    WordGraph& remove_target(node_type s, label_type a);

    // TODO(doc)
    WordGraph& remove_label(label_type a);

    // TODO(doc)
    WordGraph& remove_label_no_checks(label_type a);

    //! Remove all of the edges in the word graph.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes and \p n is the
    //! out-degree.
    //!
    inline WordGraph& remove_all_targets() {
      std::fill(_dynamic_array_2.begin(), _dynamic_array_2.end(), UNDEFINED);
      return *this;
    }

    //! Swaps the edge with specified label from one node with another.
    //!
    //! This function swaps the target of the edge from the node \p u labelled
    //! \p a with the target of the edge from the node \p v labelled \p a.
    //!
    //! \param u the first node
    //! \param v the second node
    //! \param a the label
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    // swap u - a - > u' and v - a -> v'
    WordGraph& swap_targets_no_checks(node_type u, node_type v, label_type a) {
      _dynamic_array_2.swap(u, a, v, a);
      return *this;
    }

    // TODO(doc)
    WordGraph& swap_targets(node_type u, node_type v, label_type a);

    //! Check if two word graphs are equal.
    //!
    //! \param that the word graph for comparisonb
    //!
    //! \returns
    //! A `bool`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    // swap u - a - > u' and v - a -> v'
    [[nodiscard]] bool operator==(WordGraph const& that) const {
      return _dynamic_array_2 == that._dynamic_array_2;
    }

    [[nodiscard]] bool operator!=(WordGraph const& that) const {
      return !operator==(that);
    }

    [[nodiscard]] bool operator<(WordGraph const& that) const {
      return _dynamic_array_2 < that._dynamic_array_2;
    }

    [[nodiscard]] bool operator<=(WordGraph const& that) const {
      return _dynamic_array_2 <= that._dynamic_array_2;
    }

    [[nodiscard]] bool operator>(WordGraph const& that) const {
      return _dynamic_array_2 > that._dynamic_array_2;
    }

    [[nodiscard]] bool operator>=(WordGraph const& that) const {
      return _dynamic_array_2 > that._dynamic_array_2;
    }

    ////////////////////////////////////////////////////////////////////////
    // WordGraph - nodes, targets, etc - public
    ////////////////////////////////////////////////////////////////////////

    //! Get the range of the edge with given source node and label.
    //!
    //! \param v the node
    //! \param lbl the label
    //!
    //! \returns
    //! Returns the node adjacent to \p v via the edge labelled \p lbl, or
    //! libsemigroups::UNDEFINED; both are values of type \ref node_type.
    //!
    //! \complexity
    //! Constant.
    //!
    //! \exception LibsemigroupsException if \p v or \p lbl is not
    //! valid.
    // Not noexcept because throw_if_node_out_of_bounds/label aren't
    [[nodiscard]] node_type target(node_type v, label_type lbl) const;

    //! Get the range of the edge with given source node and label.
    //!
    //! \param v the node
    //! \param lbl the label
    //!
    //! \returns
    //! Returns the node adjacent to \p v via the edge labelled \p lbl, or
    //! libsemigroups::UNDEFINED; both are values of type \ref node_type.
    //!
    //! \complexity
    //! Constant.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning This function is unsafe because it does not verify \p v or \p
    //! lbl is valid.
    // Not noexcept because DynamicArray2::get is not
    [[nodiscard]] node_type inline target_no_checks(node_type  v,
                                                    label_type lbl) const {
      return _dynamic_array_2.get(v, lbl);
    }

    //! Get the next target of a node that doesn't equal
    //! libsemigroups::UNDEFINED.
    //!
    //! If `target(v, i)` equals libsemigroups::UNDEFINED for every value in
    //! the range \f$[i, n)\f$, where \f$n\f$ is the return value of
    //! out_degree() then \c x.first and \c x.second equal
    //! libsemigroups::UNDEFINED.
    //!
    //! \param v the node
    //! \param i the label
    //!
    //! \returns
    //! Returns a std::pair
    //! \c x where:
    //! 1. \c x.first is adjacent to \p v via an edge labelled
    //!    \c x.second;
    //! 2. \c x.second is the minimum value in the range \f$[i, n)\f$ such that
    //!    `target(v, x.second)` is not equal to libsemigroups::UNDEFINED
    //!    where \f$n\f$ is the return value of out_degree(); and
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ equals out_degree().
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \sa next_label_and_target.
    //!
    //! \warning This function is unsafe because it is not verified that the
    //! parameter \p v represents a node of \c this.
    // Not noexcept because DynamicArray2::get is not
    [[nodiscard]] std::pair<label_type, node_type>
        next_label_and_target_no_checks(node_type, label_type) const;

    //! Get the next target of a node that doesn't equal
    //! libsemigroups::UNDEFINED.
    //!
    //! If `target(v, i)` equals libsemigroups::UNDEFINED for every value in
    //! the range \f$[i, n)\f$, where \f$n\f$ is the return value of
    //! out_degree() then \c x.first and \c x.second equal
    //! libsemigroups::UNDEFINED.
    //!
    //! \param v the node
    //! \param i the label
    //!
    //! \returns
    //! Returns a std::pair
    //! \c x where:
    //! 1. \c x.first is adjacent to \p v via an edge labelled
    //!    \c x.second; and
    //! 2. \c x.second is the minimum value in the range \f$[i, n)\f$ such that
    //!    `target(v, x.second)` is not equal to libsemigroups::UNDEFINED
    //!    where \f$n\f$ is the return value of out_degree();
    //!
    //! If no such value exists, then `{UNDEFINED, UNDEFINED}` is returned.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ equals out_degree().
    //!
    //! \throws LibsemigroupsException if \p v does not represent a node in \c
    //! this.
    //!
    //! \sa next_label_and_target_no_checks.
    // Not noexcept because next_label_and_target_no_checks is not
    [[nodiscard]] std::pair<label_type, node_type>
    next_label_and_target(node_type v, label_type i) const;

    //! Returns the number of nodes.
    //!
    //! \returns
    //! The number of nodes, a value of type \c Node.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] size_type inline number_of_nodes() const noexcept {
      return _nr_nodes;
    }

#ifndef PARSED_BY_DOXYGEN
    WordGraph& number_of_active_nodes(size_type val) {
      _num_active_nodes = val;
      return *this;
    }

    [[nodiscard]] size_type inline number_of_active_nodes() const noexcept {
      return _num_active_nodes;
    }
#endif

    //! Returns the number of edges.
    //!
    //! \returns
    //! The total number of edges, a value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is number_of_nodes() and \c n is out_degree().
    //!
    // Not noexcept because std::count isn't
    [[nodiscard]] size_type number_of_edges() const;

    //! Returns the number of edges with given source node.
    //!
    //! \param n the node.
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \throws LibsemigroupsException if \p n is not a node.
    //!
    //! \complexity
    //! \f$O(n)\f$ where \c n is out_degree().
    [[nodiscard]] size_type number_of_edges(node_type n) const;

    //! TODO(doc)
    [[nodiscard]] size_type number_of_edges_no_checks(node_type n) const;

    //! Returns the out-degree.
    //!
    //! \returns
    //! The number of out-edges of every node, a value of type \c Node.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] size_type out_degree() const noexcept {
      return _degree;
    }

    //! Returns a random access iterator pointing at the first node of the
    //! word graph.
    //!
    //! \returns
    //! An \ref const_iterator_nodes.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    //!
    const_iterator_nodes cbegin_nodes() const noexcept {
      return detail::IntRange<node_type>(0, number_of_nodes()).cbegin();
    }

    // no except correct?
    [[nodiscard]] auto nodes() const noexcept {
      return rx::seq<node_type>() | rx::take(number_of_nodes());
    }

    // no except correct?
    [[nodiscard]] auto labels() const noexcept {
      return rx::seq<label_type>() | rx::take(out_degree());
    }

    //! Returns a random access iterator pointing at the last node of the
    //! word graph.
    //!
    //! \returns
    //! An \ref const_reverse_iterator_nodes.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_reverse_iterator_nodes crbegin_nodes() const noexcept {
      return detail::IntRange<node_type>(0, number_of_nodes()).crbegin();
    }

    //! Returns a random access iterator pointing one-past-the-first node of
    //! the word graph.
    //!
    //! \returns
    //! An \ref const_reverse_iterator_nodes.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    //!
    [[nodiscard]] const_reverse_iterator_nodes crend_nodes() const noexcept {
      return detail::IntRange<node_type>(0, number_of_nodes()).crend();
    }

    //! Returns a random access iterator pointing one-past-the-last node of the
    //! word graph.
    //!
    //! \returns
    //! An \ref const_iterator_nodes.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    //!
    [[nodiscard]] const_iterator_nodes cend_nodes() const noexcept {
      return detail::IntRange<node_type>(0, number_of_nodes()).cend();
    }

    //! Returns a random access iterator pointing at the first target of a
    //! node.
    //!
    //! \param i a node in the word graph.
    //!
    //! \returns
    //! A \ref const_iterator_targets.
    //!
    //! \throws LibsemigroupsException if \p i is not valid.
    //!
    //! \complexity
    //! Constant.
    // Not noexcept because throw_if_node_out_of_bounds isn't
    [[nodiscard]] const_iterator_targets cbegin_targets(node_type i) const;

    //! Returns a random access iterator pointing at the first target of a
    //! node.
    //!
    //! \param i a node in the word graph.
    //!
    //! \returns
    //! A \ref const_iterator_targets.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    //!
    //! \sa
    //! \ref cbegin_targets.
    [[nodiscard]] const_iterator_targets
    cbegin_targets_no_checks(node_type i) const noexcept {
      return _dynamic_array_2.cbegin_row(i);
    }

    //! Returns a random access iterator pointing one-past-the-last target of
    //! a node.
    //!
    //! \param i a node in the word graph.
    //!
    //! \returns
    //! A \ref const_iterator_targets.
    //!
    //! \throws LibsemigroupsException if \p i is not valid.
    //!
    //! \complexity
    //! Constant.
    // Not noexcept because throw_if_node_out_of_bounds isn't
    [[nodiscard]] const_iterator_targets cend_targets(node_type i) const;

    //! Returns a random access iterator pointing one-past-the-last target of
    //! a node.
    //!
    //! \param i a node in the word graph.
    //!
    //! \returns
    //! A \ref const_iterator_targets.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    //!
    //! \sa
    //! \ref cend_targets.
    [[nodiscard]] const_iterator_targets
    cend_targets_no_checks(node_type i) const noexcept {
      return _dynamic_array_2.cbegin_row(i) + _degree;
    }

    // TODO(doc)
    [[nodiscard]] rx::iterator_range<const_iterator_targets>
    targets_no_checks(node_type n) const noexcept {
      return rx::iterator_range(cbegin_targets_no_checks(n),
                                cend_targets_no_checks(n));
    }

    // TODO(doc)
    [[nodiscard]] rx::iterator_range<const_iterator_targets>
    targets(node_type n) const;

    // TODO(doc)
    [[nodiscard]] auto
    labels_and_targets_no_checks(node_type n) const noexcept {
      return rx::enumerate(targets_no_checks(n));
    }

    // TODO(doc)
    [[nodiscard]] auto labels_and_targets(node_type n) const;

    //! Restrict the word graph to its first \p n nodes.
    //!
    //! \param n the number of nodes
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! This function performs no checks whatsoever and will result in a
    //! corrupted word graph if there are any edges from the nodes \f$0, \ldots,
    //! n - 1\f$ to nodes larger than \f$n - 1\f$.
    // TODO(0) update doc
    WordGraph& induced_subgraph_no_checks(node_type first, node_type last);

    // TODO(doc)
    WordGraph& induced_subgraph(node_type first, node_type last);

    // TODO(doc)
    template <typename Iterator,
              typename = std::enable_if_t<detail::IsIterator<Iterator>::value>>
    WordGraph& induced_subgraph_no_checks(Iterator first, Iterator last);

    // TODO(doc)
    //! \throws LibsemigroupsException if any value in the range specified by \p
    //! first and \p last is not a node of the word graph.
    //!
    //! \throws LibsemigroupsException if any target of any edge with
    //! source node in the range specified by \p first and \p last does not
    //! belong to the same range.
    template <typename Iterator,
              typename = std::enable_if_t<detail::IsIterator<Iterator>::value>>
    WordGraph& induced_subgraph(Iterator first, Iterator last);

    WordGraph& disjoint_union_inplace_no_checks(WordGraph<Node> const& that);

    WordGraph& disjoint_union_inplace(WordGraph<Node> const& that);

    // TODO(doc)
    // requires access to apply_row_permutation so can't be helper
    WordGraph& permute_nodes_no_checks(std::vector<node_type> const& p,
                                       std::vector<node_type> const& q,
                                       size_t                        m);

    // TODO(doc)
    WordGraph& permute_nodes_no_checks(std::vector<node_type> const& p,
                                       std::vector<node_type> const& q) {
      return permute_nodes_no_checks(p, q, p.size());
    }

   protected:
    // from WordGraphWithSources
    template <typename S>
    void apply_row_permutation(S const& p) {
      _dynamic_array_2.apply_row_permutation(p);
    }

   private:
    ////////////////////////////////////////////////////////////////////////
    // WordGraph - data members - private
    ////////////////////////////////////////////////////////////////////////

    size_type _degree;
    size_type _nr_nodes;
    // TODO(later) remove when WordGraphView is implemented
    size_type                           _num_active_nodes;
    mutable detail::DynamicArray2<Node> _dynamic_array_2;
  };

  //! \ingroup word_graph_group
  //!
  //! \brief Namespace containing helper functions for the \ref WordGraph
  //! class.
  //!
  //! Defined in ``word-graph.hpp``.
  //!
  //! \brief This namespace contains helper functions for the \ref WordGraph
  //! class.
  namespace word_graph {

    //////////////////////////////////////////////////////////////////////////
    // WordGraph - helper functions - in alphabetical order!!!
    //////////////////////////////////////////////////////////////////////////

    //! \brief Adds a cycle involving the specified range of nodes to a word
    //! graph.
    //!
    //! This function adds a cycle involving the specified range of nodes.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \tparam Iterator the type of an iterator pointing to nodes of
    //! a wordGraph
    //!
    //! \param wg the WordGraph object to add a cycle to.
    //! \param first an iterator to nodes of \p wg.
    //! \param last an iterator to nodes of \p wg.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m)\f$ where \f$m\f$ is the distance between \p first and \p last.
    //!
    //! \note
    //! The edges added by this function are all labelled \c 0.
    // TODO(1) add add_cycle with checks version.
    template <typename Node, typename Iterator>
    void add_cycle_no_checks(WordGraph<Node>& wg,
                             Iterator         first,
                             Iterator         last);

    //! \brief Adds a cycle consisting of \p N new nodes.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the WordGraph object to add a cycle to.
    //! \param N the length of the cycle and number of new nodes to add.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(N)\f$ where \f$N\f$ is the second parameter.
    //!
    //! \note
    //! The edges added by this function are all labelled \c 0.
    template <typename Node>
    void add_cycle(WordGraph<Node>& wg, size_t N) {
      size_t M = wg.number_of_nodes();
      wg.add_nodes(N);
      add_cycle_no_checks(wg, wg.cbegin_nodes() + M, wg.cend_nodes());
    }

    //! \brief Returns the adjacency matrix of a word graph.
    //!
    //! This function returns the adjacency matrix of the word graph \p wg. The
    //! type of the returned matrix depends on whether or not ``libsemigroups``
    //! is compiled with [eigen][] enabled. The returned matrix has the number
    //! of edges with source \c s and target \c t in the `(s, t)`-entry.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the word graph.
    //!
    //! \returns The adjacency matrix.
    //!
    //! [eigen]: http://eigen.tuxfamily.org/
    template <typename Node>
    [[nodiscard]] auto adjacency_matrix(WordGraph<Node> const& wg);

    //! \brief Returns a \ref Dot object representing a word graph.
    //!
    //! This function returns a \ref Dot object representing the word graph \p
    //! wg.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the word graph.
    //!
    //! \returns A \ref Dot object.
    template <typename Node>
    [[nodiscard]] Dot dot(WordGraph<Node> const& wg);

    //! \brief Compares two word graphs on a range of nodes.
    //!
    //! This function returns \c true if the word graphs \p x and \p y are
    //! equal on the range of nodes from \p first to \p last; and \c false
    //! otherwise.  The word graphs \p x and \p y are equal at a node \c s if:
    //! * the out-degrees of \p x and \p y coincide;
    //! * the edges with source \c s and label \c a have equal targets in \p x
    //! and \p y for every label \c a.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param x  the first word graph for comparison.
    //! \param y the second word graph for comparison.
    //! \param first the first node in the range.
    //! \param last the last node in the range plus \c 1.
    //!
    //! \returns Whether or not the word graphs are equal at the specified range
    //! of nodes.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning No checks are performed to ensure that the arguments
    //! are valid.
    //!
    //! \sa WordGraph::operator== for a comparison of two entire word graphs.
    template <typename Node>
    [[nodiscard]] bool equal_to_no_checks(WordGraph<Node> const& x,
                                          WordGraph<Node> const& y,
                                          Node                   first,
                                          Node                   last);

    //! \brief Compares two word graphs on a range of nodes.
    //!
    //! This function returns \c true if the word graphs \p x and \p y are
    //! equal on the range of nodes from \p first to \p last; and \c false
    //! otherwise.  The word graphs \p x and \p y are equal at a node \c s if:
    //! * the out-degrees of \p x and \p y coincide;
    //! * the edges with source \c s and label \c a have equal targets in \p x
    //! and \p y for every label \c a.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param x  the first word graph for comparison.
    //! \param y the second word graph for comparison.
    //! \param first the first node in the range.
    //! \param last the last node in the range plus \c 1.
    //!
    //! \returns Whether or not the word graphs are equal at the specified range
    //! of nodes.
    //!
    //! \throw LibsemigroupsException if \p first  is not a node
    //! in \p x or not a node in \p y; or if `last - 1` is not a node in \p or
    //! not a node in \p y.
    //!
    //! \sa WordGraph::operator== for a comparison of two entire word graphs.
    template <typename Node>
    [[nodiscard]] bool equal_to(WordGraph<Node> const& x,
                                WordGraph<Node> const& y,
                                Node                   first,
                                Node                   last);

    //! \brief Find the node that a path starting at a given node leads to (if
    //! any).
    //!
    //! This function attempts to follow the path in the word graph \p wg
    //! starting at the node \p from  labelled by the word defined by \p first
    //! and \p last. If this path exists, then the last node on that path is
    //! returned. If this path does not exist, then \ref UNDEFINED is returned.
    //!
    //! \tparam Node1
    //! the type used as the template parameter for the WordGraph \p wg.
    //!
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \tparam Iterator the type of \p first and \p last.
    //!
    //! \param wg a word graph.
    //! \param from the starting node.
    //! \param first an iterator point at the start of the word.
    //! \param last an iterator point one beyond the last letter of the word.
    //!
    //! \returns
    //! A value of type \p Node1. If one or more edges in \p path are not
    //! defined, then \ref UNDEFINED is returned.
    //!
    //! \throw LibsemigroupsException if \p from is not a node in the word
    //! graph or the word defined by \p first and \p last contains a value that
    //! is not an edge-label.
    //!
    //! \par Complexity
    //! Linear in the distance between \p first and \p last.
    template <typename Node1, typename Node2, typename Iterator>
    [[nodiscard]] Node1 follow_path(WordGraph<Node1> const& wg,
                                    Node2                   from,
                                    Iterator                first,
                                    Iterator                last);

    //! \brief Find the node that a path starting at a given node leads to (if
    //! any).
    //!
    //! This function attempts to follow the path in the word graph \p wg
    //! starting at the node \p from  labelled by the word \p path. If this path
    //! exists, then the last node on that path is returned. If this path does
    //! not exist, then \ref UNDEFINED is returned.
    //!
    //! \tparam Node1
    //! the type used as the template parameter for the WordGraph \p wg.
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \param wg a word graph.
    //! \param from the starting node.
    //! \param path the path to follow.
    //!
    //! \returns
    //! A value of type \p Node1. If one or more edges in \p path are not
    //! defined, then \ref UNDEFINED is returned.
    //!
    //! \throw LibsemigroupsException if \p from is not a node in the word
    //! graph or \p path contains a value that is not an edge-label.
    //!
    //! \par Complexity
    //! Linear in the length of \p path.
    // TODO(later) example
    // not noexcept because WordGraph::target isn't
    template <typename Node1, typename Node2>
    [[nodiscard]] Node1 follow_path(WordGraph<Node1> const& wg,
                                    Node2                   from,
                                    word_type const&        path) {
      static_assert(sizeof(Node2) <= sizeof(Node1));
      return follow_path(wg, from, path.cbegin(), path.cend());
    }

    //! \brief Follow the path from a specified node labelled by a word.
    //!
    //! This function returns the last node on the path in the word graph
    //! \p wg starting at the node \p from labelled by the word defined by \p
    //! first and \p last or \ref UNDEFINED.
    //!
    //! \tparam Node1
    //! the type used as the template parameter for the WordGraph \p wg.
    //!
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \param wg a word graph.
    //! \param from the source node.
    //! \param first iterator into a word.
    //! \param last iterator into a word.
    //!
    //! \exception
    //! \noexcept
    //!
    //! \returns A value of type \p Node1.
    //!
    //! \complexity
    //! At worst the distance from \p first to \p last.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed.
    template <typename Node1, typename Node2, typename Iterator>
    [[nodiscard]] Node1 follow_path_no_checks(WordGraph<Node1> const& wg,
                                              Node2                   from,
                                              Iterator                first,
                                              Iterator last) noexcept;

    //! \brief Follow the path from a specified node labelled by a word.
    //!
    //! This function returns the last node on the path in the word graph
    //! \p wg starting at the node \p from labelled by \p path or
    //! \ref UNDEFINED.
    //!
    //! \tparam Node1
    //! the type used as the template parameter for the WordGraph \p wg.
    //!
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \param wg a word graph.
    //! \param from the source node.
    //! \param path the word.
    //!
    //! \exception
    //! \noexcept
    //!
    //! \returns A value of type \p Node1.
    //!
    //! \complexity
    //! At worst the length of \p path.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed.
    template <typename Node1, typename Node2>
    [[nodiscard]] Node1 follow_path_no_checks(WordGraph<Node1> const& wg,
                                              Node2                   from,
                                              word_type const& path) noexcept {
      static_assert(sizeof(Node2) <= sizeof(Node1));
      return follow_path_no_checks(wg, from, path.cbegin(), path.cend());
    }

    //! \brief Check if a word graph is acyclic.
    //!
    //! This function returns \c true if the word graph \p wg is acyclic and \c
    //! false otherwise. A word graph is acyclic if every directed cycle in the
    //! word graph is trivial.
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the WordGraph object to check.
    //!
    //! \returns
    //! A value of type `bool`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the \ref WordGraph::out_degree.
    //!
    //!
    //! \par Example
    //! \code
    //! WordGraph<size_t> wg;
    //! wg.add_nodes(2);
    //! wg.add_to_out_degree(1);
    //! wg.target(0, 0, 1);
    //! wg.target(1, 0, 0);
    //! word_graph::is_acyclic(wg); // returns false
    //! \endcode
    // Not noexcept because detail::is_acyclic isn't
    template <typename Node>
    [[nodiscard]] bool is_acyclic(WordGraph<Node> const& wg);

    //! \brief Check if the word graph induced by the nodes reachable from a
    //! source node is acyclic.
    //!
    //! This function returns \c true if the word graph consisting of the nodes
    //! reachable from \p source in the word graph \p wg is acyclic and \c false
    //! if not. A word graph is acyclic if every directed cycle in the word
    //! graph is trivial.
    //!
    //! \tparam Node1
    //! the type used as the template parameter for the WordGraph \p wg.
    //!
    //! \tparam Node2 the type of the node \p source (must satisfy
    //! `sizeof(Node2) <= sizeof(Node1)`).
    //!
    //! \param wg the WordGraph object to check.
    //! \param source the source node.
    //!
    //! \returns
    //! A value of type `bool`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the \ref WordGraph::out_degree.
    //!
    //!
    //! \par Example
    //! \code
    //! WordGraph<size_t> wg;
    //! wg.add_nodes(4);
    //! wg.add_to_out_degree(1);
    //! wg.target(0, 0, 1);
    //! wg.target(1, 0, 0);
    //! wg.target(2, 0, 3);
    //! word_graph::is_acyclic(wg); // returns false
    //! word_graph::is_acyclic(wg, 0); // returns false
    //! word_graph::is_acyclic(wg, 1); // returns false
    //! word_graph::is_acyclic(wg, 2); // returns true
    //! word_graph::is_acyclic(wg, 3); // returns true
    //! \endcode
    // Not noexcept because detail::is_acyclic isn't
    template <typename Node1, typename Node2>
    [[nodiscard]] bool is_acyclic(WordGraph<Node1> const& wg, Node2 source);

    //! \brief Check if the word graph induced by the nodes reachable from a
    //! source node and from which a target node can be reached is acyclic.
    //!
    //! This function returns \c true if the word graph consisting of the nodes
    //! reachable from \p source and from which \p target is reachable, in the
    //! word graph \p wg, is acyclic; and \c false if not. A word graph is
    //! acyclic if every directed cycle of the word graph is trivial.
    //!
    //! \tparam Node1
    //! the type used as the template parameter for the WordGraph \p wg.
    //!
    //! \tparam Node2 the type of the nodes \p source and \p target (must
    //! satisfy `sizeof(Node2) <= sizeof(Node1)`).
    //!
    //! \param wg the WordGraph object to check.
    //! \param source the source node.
    //! \param target the target node.
    //!
    //! \returns
    //! A value of type `bool`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the \ref WordGraph::out_degree.
    //!
    //! \endcode
    // Not noexcept because detail::is_acyclic isn't
    template <typename Node1, typename Node2>
    [[nodiscard]] bool is_acyclic(WordGraph<Node1> const& wg,
                                  Node2                   source,
                                  Node2                   target);

    //! \brief Check if a word graph is compatible with some relations at a
    //! range of nodes.
    //!
    //! This function returns \c true if the word graph \p wg is compatible
    //! with the relations in the range \p first_rule to \p last_rule at every
    //! node in the range from \p first_node to \p last_node. This means that
    //! the paths with given sources that are labelled by one side of a relation
    //! leads to the same node as the path labelled by the other side of the
    //! relation.
    //!
    //! \tparam Node the type used as the template parameter for the WordGraph
    //! \p wg.
    //!
    //! \tparam Iterator1 the type of \p first_node.
    //!
    //! \tparam Iterator2 the type of \p last_node.
    //!
    //! \tparam Iterator3 the type of \p first_rule and \p last_rule.
    //!
    //! \param wg the word graph.
    //!
    //! \param first_node iterator pointing at the first node.
    //!
    //! \param last_node iterator pointing at one beyond the last node.
    //!
    //! \param first_rule iterator pointing to the first rule.
    //!
    //! \param last_rule iterator pointing one beyond the last rule.
    //!
    //! \return Whether or not the word graph is compatible with the given rules
    //! at each one of the given nodes.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    template <typename Node,
              typename Iterator1,
              typename Iterator2,
              typename Iterator3>
    [[nodiscard]] bool is_compatible_no_checks(WordGraph<Node> const& wg,
                                               Iterator1 first_node,
                                               Iterator2 last_node,
                                               Iterator3 first_rule,
                                               Iterator3 last_rule);

    //! \brief Check if a word graph is compatible with some relations at a
    //! range of nodes.
    //!
    //! This function returns \c true if the word graph \p wg is compatible
    //! with the relations in the range \p first_rule to \p last_rule at every
    //! node in the range from \p first_node to \p last_node. This means that
    //! the paths with given sources that are labelled by one side of a relation
    //! leads to the same node as the path labelled by the other side of the
    //! relation.
    //!
    //! \tparam Node the type used as the template parameter for the WordGraph
    //! \p wg.
    //!
    //! \tparam Iterator1 the type of \p first_node.
    //!
    //! \tparam Iterator2 the type of \p last_node.
    //!
    //! \tparam Iterator3 the type of \p first_rule and \p last_rule.
    //!
    //! \param wg the word graph.
    //!
    //! \param first_node iterator pointing at the first node.
    //!
    //! \param last_node iterator pointing at one beyond the last node.
    //!
    //! \param first_rule iterator pointing to the first rule.
    //!
    //! \param last_rule iterator pointing one beyond the last rule.
    //!
    //! \return Whether or not the word graph is compatible with the given rules
    //! at each one of the given nodes.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed.
    //!
    //! \throws LibsemigroupsException if any of the nodes in the range between
    //! \p first_node and \p last_node does not belong to \p wg (i.e. is greater
    //! than or equal to WordGraph::number_of_nodes).
    //!
    //! \throws LibsemigroupsException if any of the rules in the range between
    //! \p first_rule and \p last_rule contains an invalid label (i.e. one
    //! greater than or equal to WordGraph::out_degree).
    template <typename Node,
              typename Iterator1,
              typename Iterator2,
              typename Iterator3>
    [[nodiscard]] bool is_compatible(WordGraph<Node> const& wg,
                                     Iterator1              first_node,
                                     Iterator2              last_node,
                                     Iterator3              first_rule,
                                     Iterator3              last_rule);

    //! \brief Check if every node in a range has exactly \ref out_degree
    //! out-edges.
    //!
    //! This function returns \c true if every node in the range defined by \p
    //! first_node and \p last_node is complete, meaning that every such node is
    //! the source of an edge with every possible label.
    //!
    //! \tparam Node the type of the nodes in the word graph.
    //!
    //! \tparam Iterator1 the type of \p first_node.
    //!
    //! \tparam Iterator2 the type of \p last_node.
    //!
    //! \param wg the word graph.
    //!
    //! \param first_node iterator pointing to the first node in the range.
    //!
    //! \param last_node iterator pointing one beyond the last node in the
    //! range.
    //!
    //! \returns
    //! Whether or not the word graph is complete on the given range of nodes.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is the number of nodes in the range and \c n is
    //! out_degree().
    //!
    //! \warning No checks are performed on the arguments.
    template <typename Node, typename Iterator1, typename Iterator2>
    [[nodiscard]] bool is_complete_no_checks(WordGraph<Node> const& wg,
                                             Iterator1              first_node,
                                             Iterator2              last_node);

    //! \brief Check if every node in a range has exactly \ref out_degree
    //! out-edges.
    //!
    //! This function returns \c true if every node in the range defined by \p
    //! first_node and \p last_node is complete, meaning that
    //! every such node is the source of an edge with every possible label.
    //!
    //! \tparam Node the type of the nodes in the word graph.
    //!
    //! \tparam Iterator1 the type of \p first_node.
    //!
    //! \tparam Iterator2 the type of \p last_node.
    //!
    //! \param wg the word graph.
    //!
    //! \param first_node iterator pointing to the first node in the range.
    //!
    //! \param last_node iterator pointing one beyond the last node in the
    //! range.
    //!
    //! \returns Whether or not the word graph is complete on the given
    //! range of nodes.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is the number of nodes in the range and \c n is
    //! out_degree().
    //!
    //! \throws LibsemigroupsException if any item in the range defined by \p
    //! first_node and \p last_node is not a node of \p wg.
    template <typename Node, typename Iterator1, typename Iterator2>
    [[nodiscard]] bool is_complete(WordGraph<Node> const& wg,
                                   Iterator1              first_node,
                                   Iterator2              last_node);

    //! \brief Check if every node has exactly \ref out_degree out-edges.
    //!
    //! This function returns \c true if a WordGraph is complete, meaning that
    //! every node is the source of an edge with every possible label.
    //!
    //! \tparam Node the type of the nodes in the word graph.
    //!
    //! \param wg the word graph.
    //!
    //! \returns
    //! Whether or not the word graph is complete.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is number_of_nodes() and \c n is out_degree().
    template <typename Node>
    [[nodiscard]] bool is_complete(WordGraph<Node> const& wg) noexcept {
      return wg.number_of_edges() == wg.number_of_nodes() * wg.out_degree();
    }

    //! \brief Check if a word graph is connected.
    //!
    //! This function returns \c true if the word graph \p wg is connected and
    //! \c false if it is not. A word graph is *connected* if for every pair of
    //! nodes \c s and \c t in the graph there exists a sequence \f$u_0 = s,
    //! \ldots, u_{n}= t\f$ for some \f$n\in \mathbb{N}\f$ such that for every
    //! \f$i\f$ there exists a label \c a such that \f$(u_i, a, u_{i + 1})\f$ or
    //! \f$(u_{i + 1}, a, u_i)\f$ is an edge in the graph.
    //!
    //! \tparam Node the type of the nodes in the word graph.
    //!
    //! \param wg the word graph.
    //!
    //! \returns
    //! Whether or not the word graph is connected.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! The function does not check it arguments. In particular, if any target
    //! in \p wg is out of bounds, then bad things will happen.
    template <typename Node>
    [[nodiscard]] bool is_connected_no_checks(WordGraph<Node> const& wg);

    //! \brief Check if a word graph is connected.
    //!
    //! This function returns \c true if the word graph \p wg is connected and
    //! \c false if it is not. A word graph is *connected* if for every pair of
    //! nodes \c s and \c t in the graph there exists a sequence \f$u_0 = s,
    //! \ldots, u_{n}= t\f$ for some \f$n\in \mathbb{N}\f$ such that for every
    //! \f$i\f$ there exists a label \c a such that \f$(u_i, a, u_{i + 1})\f$ or
    //! \f$(u_{i + 1}, a, u_i)\f$ is an edge in the graph.
    //!
    //! \tparam Node the type of the nodes in the word graph.
    //!
    //! \param wg the word graph.
    //!
    //! \returns
    //! Whether or not the word graph is connected.
    //!
    //! \throws LibsemigroupsException if any target in \p wg is out of bounds,
    //! i.e. if any target \c t is not equal to \ref UNDEFINED and not in the
    //! nodes of \p wg.
    // TODO(0) in the python bindings use the no_checks  versions because we
    // can't construct word graphs there that contain bad targets.
    template <typename Node>
    [[nodiscard]] bool is_connected(WordGraph<Node> const& wg);

    //! \brief Check if there is a path from one node to another.
    //!
    //! This function returns \c true if there is a path from the node \p source
    //! to the node \p target in the word graph \p wg.
    //!
    //! \tparam Node1 the type of the nodes in the WordGraph.
    //!
    //! \tparam Node 2 the types of \p source and \p target (must
    //! satisfy `sizeof(Node2) <= sizeof(Node1)`).
    //!
    //! \param wg the WordGraph object to check.
    //! \param source the source node.
    //! \param target the source node.
    //!
    //! \returns
    //! Whether or not the node \p target is reachable from the node \p source
    //! in the word graph \p wg.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the WordGraph::out_degree.
    //!
    //! \note
    //! If \p source and \p target are equal, then, by convention, we consider
    //! \p target to be reachable from \p source, via the empty path.
    //!
    //! \warning No checks are performed on the arguments.
    //!
    //! \par Example
    //! \code
    //! WordGraph<size_t> wg;
    //! wg.add_nodes(4);
    //! wg.add_to_out_degree(1);
    //! wg.target(0, 1, 0);
    //! wg.target(1, 0, 0);
    //! wg.target(2, 3, 0);
    //! word_graph::is_reachable_no_checks(wg, 0, 1); // returns true
    //! word_graph::is_reachable_no_checks(wg, 1, 0); // returns true
    //! word_graph::is_reachable_no_checks(wg, 1, 2); // returns false
    //! word_graph::is_reachable_no_checks(wg, 2, 3); // returns true
    //! word_graph::is_reachable_no_checks(wg, 3, 2); // returns false
    //! \endcode
    template <typename Node1, typename Node2>
    [[nodiscard]] bool is_reachable_no_checks(WordGraph<Node1> const& wg,
                                              Node2                   source,
                                              Node2                   target);

    //! \brief Check if there is a path from one node to another.
    //!
    //! This function returns \c true if there is a path from the node \p source
    //! to the node \p target in the word graph \p wg.
    //!
    //! \tparam Node1 the type of the nodes in the WordGraph.
    //!
    //! \tparam Node 2 the types of \p source and \p target (must
    //! satisfy `sizeof(Node2) <= sizeof(Node1)`).
    //!
    //! \param wg the WordGraph object to check.
    //! \param source the source node.
    //! \param target the source node.
    //!
    //! \returns
    //! Whether or not the node \p target is reachable from the node \p source
    //! in the word graph \p wg.
    //!
    //! \throws LibsemigroupsException if \p source or \p target is out of
    //! bounds.
    //! \throws LibsemigroupsException if any target in \p wg is out of bounds.
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the \ref out_degree.
    //!
    //! \note
    //! If \p source and \p target are equal, then, by convention, we consider
    //! \p target to be reachable from \p source, via the empty path.
    // TODO(0) in the python bindings use the no_checks version after checking
    // that the source and target are ok, because we can't create a word graph
    // with bad targets there, so no point checking that.
    template <typename Node1, typename Node2>
    [[nodiscard]] bool is_reachable(WordGraph<Node1> const& wg,
                                    Node2                   source,
                                    Node2                   target);

    //! \brief Check if every node is reachable from some node.
    //!
    //! This function returns \c true if there exists a node in \p wg from
    //! which every other node is reachable; and \c false otherwise.
    //! A word graph is *strictly cyclic* if there exists a node \f$v\f$ from
    //! which every node is reachable (including \f$v\f$). There must be a
    //! path of length at least \f$1\f$ from the original node \f$v\f$ to
    //! itself (i.e. \f$v\f$ is not considered to be reachable from itself by
    //! default).
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the WordGraph object to check.
    //!
    //! \returns
    //! A value of type `bool`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the WordGraph::out_degree.
    //!
    //! \warning This function does not check its arguments. In particular, if
    //! \p wg has any targets that are out of bounds, then bad things will
    //! happen.
    //!
    //! \par Example
    //! \code
    //! auto wg = to_word_graph<uint8_t>(
    //!     5, {{0, 0}, {1, 1}, {2}, {3, 3}});
    //! word_graph::is_strictly_cyclic(wg);  // returns false
    //! \endcode
    // TODO(0) again use the no_checks version in the python bindings.
    template <typename Node>
    [[nodiscard]] bool is_strictly_cyclic_no_checks(WordGraph<Node> const& wg);

    //! \brief Check if every node is reachable from some node.
    //!
    //! This function returns \c true if there exists a node in \p wg from
    //! which every other node is reachable; and \c false otherwise.
    //! A word graph is *strictly cyclic* if there exists a node \f$v\f$ from
    //! which every node is reachable (including \f$v\f$). There must be a
    //! path of length at least \f$1\f$ from the original node \f$v\f$ to
    //! itself (i.e. \f$v\f$ is not considered to be reachable from itself by
    //! default).
    //!
    //! \tparam Node the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the WordGraph object to check.
    //!
    //! \returns
    //! A value of type `bool`.
    //!
    //! \throws LibsemigroupsException if any target in \p wg is out of bounds.
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the WordGraph::out_degree.
    //!
    //! \par Example
    //! \code
    //! auto wg = to_word_graph<uint8_t>(
    //!     5, {{0, 0}, {1, 1}, {2}, {3, 3}});
    //! word_graph::is_strictly_cyclic(wg);  // returns false
    //! \endcode
    // TODO(0) should return the node that everything is reachable from
    template <typename Node>
    [[nodiscard]] bool is_strictly_cyclic_no_checks(WordGraph<Node> const& wg);

    //////////////////////////////////////////////////////////////////////////
    // FRONTIER
    //////////////////////////////////////////////////////////////////////////

    //! \brief Blurb
    //!
    //! This function
    // not noexcept because it throws an exception!
    template <typename Node1, typename Node2>
    void throw_if_node_out_of_bounds(WordGraph<Node1> const& wg, Node2 v);

    //! TODO(doc)
    template <typename Node, typename Iterator>
    void throw_if_node_out_of_bounds(WordGraph<Node> const& wg,
                                     Iterator               first,
                                     Iterator               last);

    //! TODO(doc)
    template <typename Node>
    void throw_if_any_target_out_of_bounds(WordGraph<Node> const& wg);

    //! TODO(doc)
    template <typename Node, typename Iterator>
    void throw_if_any_target_out_of_bounds(WordGraph<Node> const& wg,
                                           Iterator               first,
                                           Iterator               last);

    //! TODO(doc)
    // not noexcept because it throws an exception!
    template <typename Node>
    void throw_if_label_out_of_bounds(WordGraph<Node> const&               wg,
                                      typename WordGraph<Node>::label_type lbl);

    template <typename Node>
    void throw_if_label_out_of_bounds(WordGraph<Node> const& wg,
                                      word_type const&       word);

#ifdef LIBSEMIGROUPS_EIGEN_ENABLED
    static Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
    pow(Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> const& x,
        size_t                                                       e);
#endif

    //////////////////////////////////////////////////////////////////////////
    // WordGraph - helper functions - operations
    //////////////////////////////////////////////////////////////////////////

    // TODO(doc)
    // TODO(1) tests
    // TODO(1) version where std::unordered_set is passed by reference, or make
    // this a class that stores its stack and unordered_set, not clear why we'd
    // single out the unordered_set to be passed by reference.
    // TODO(2) version which is an iterator i.e. returns an iterator or range
    // object that allows use to step through the nodes reachable from a given
    // node
    template <typename Node1, typename Node2>
    [[nodiscard]] std::unordered_set<Node1>
    nodes_reachable_from(WordGraph<Node1> const& wg, Node2 source);

    template <typename Node1, typename Node2>
    [[nodiscard]] std::unordered_set<Node1>
    nodes_reachable_from_no_checks(WordGraph<Node1> const& wg, Node2 source);

    // TODO(doc)
    // TODO(1) tests
    template <typename Node1, typename Node2>
    [[nodiscard]] size_t
    number_of_nodes_reachable_from(WordGraph<Node1> const& wg, Node2 source) {
      return nodes_reachable_from(wg, source).size();
    }

    template <typename Node1, typename Node2>
    [[nodiscard]] size_t
    number_of_nodes_reachable_from_no_checks(WordGraph<Node1> const& wg,
                                             Node2                   source) {
      return nodes_reachable_from_no_checks(wg, source).size();
    }

    //! Returns the last node on the path labelled by a word and an iterator to
    //! the position in the word reached.
    //!
    //! \tparam T the node type of the word graph
    //! \tparam S the type of the iterators into a word
    //!
    //! \param wg a word graph
    //! \param from the source node
    //! \param first iterator into a word
    //! \param last iterator into a word.
    //!
    //! \exception
    //! \noexcept
    //!
    //! \returns A pair consisting of WordGraph::node_type and \p S.
    //!
    //! \complexity
    //! At worst the distance from \p first to \p last.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed.
    template <typename Node1, typename Node2, typename Iterator>
    [[nodiscard]] std::pair<Node1, Iterator>
    last_node_on_path_no_checks(WordGraph<Node1> const& wg,
                                Node2                   from,
                                Iterator                first,
                                Iterator                last) noexcept;

    //! TODO(doc)
    template <typename Node1, typename Node2>
    std::pair<Node1, word_type::const_iterator>
    last_node_on_path(WordGraph<Node1> const& wg,
                      Node2                   from,
                      word_type const&        w);

    //! Returns the last node on the path labelled by a word and an iterator to
    //! the position in the word reached.
    //!
    //! \tparam T the node type of the word graph
    //! \tparam S the type of the iterators into a word
    //!
    //! \param wg a word graph
    //! \param from the source node
    //! \param first iterator into a word
    //! \param last iterator into a word.
    //!
    //! \throws LibsemigroupsException if any of the letters in word `[first,
    //! last)` is out of bounds.
    //!
    //! \returns A pair consisting of WordGraph::node_type and \p S.
    //!
    //! \complexity
    //! At worst the distance from \p first to \p last.
    template <typename Node1, typename Node2, typename Iterator>
    [[nodiscard]] std::pair<Node1, Iterator>
    last_node_on_path(WordGraph<Node1> const& wg,
                      Node2                   from,
                      Iterator                first,
                      Iterator                last);

    //! Returns the nodes of the word graph in topological order (see below) if
    //! possible.
    //!
    //! If it is not empty, the returned vector has the property that if an
    //! edge from a node \c n points to a node \c m, then \c m occurs before
    //! \c n in the vector.
    //!
    //! \tparam T the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the WordGraph object to check.
    //!
    //! \returns
    //! A std::vector<WordGraph<T>::node_type> that contains the nodes of
    //! \p wg in topological order (if possible) and is otherwise empty.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! WordGraph \p wg and \f$n\f$ is the number of edges. Note that for
    //! WordGraph objects the number of edges is always at most \f$mk\f$
    //! where \f$k\f$ is the WordGraph::out_degree.
    template <typename Node>
    [[nodiscard]] std::vector<Node> topological_sort(WordGraph<Node> const& wg);

    //!
    //! Returns the nodes of the word graph reachable from a given node in
    //! topological order (see below) if possible.
    //!
    //! If it is not empty, the returned vector has the property that
    //! if an edge from a node \c n points to a node \c m, then \c m occurs
    //! before \c n in the vector, and the last item in the vector is \p
    //! source.
    //!
    //! \tparam T the type used as the template parameter for the
    //! WordGraph.
    //!
    //! \param wg the WordGraph object to check.
    //! \param source the source node.
    //!
    //! \returns
    //! A std::vector<WordGraph<T>::node_type> that contains the nodes of
    //! \p wg in topological order (if possible) and is otherwise empty.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! At worst \f$O(m + n)\f$ where \f$m\f$ is the number of nodes in the
    //! subword graph of those nodes reachable from \p source
    //! and \f$n\f$ is the number of edges.
    template <typename Node1, typename Node2>
    [[nodiscard]] std::vector<Node1>
    topological_sort(WordGraph<Node1> const& wg, Node2 source);

    //////////////////////////////////////////////////////////////////////////
    // WordGraph - helper functions - properties
    //////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////
    // WordGraph - helper functions - modifiers
    //////////////////////////////////////////////////////////////////////////

    //! TODO(doc)
    // Return value indicates whether or not the graph was modified.
    // Not nodiscard because sometimes we just don't want the output
    template <typename Graph>
    bool standardize(Graph& wg, Forest& f, Order val);

    //! TODO(doc)
    // Not nodiscard because sometimes we just don't want the output
    template <typename Graph>
    std::pair<bool, Forest> standardize(Graph& wg, Order val = Order::shortlex);

    template <typename Node1, typename Node2>
    void spanning_tree_no_checks(WordGraph<Node1> const& wg,
                                 Node2                   root,
                                 Forest&                 f);

    template <typename Node1, typename Node2>
    void spanning_tree(WordGraph<Node1> const& wg, Node2 root, Forest& f);

    //! TODO(doc)
    template <typename Node1, typename Node2>
    [[nodiscard]] Forest spanning_tree(WordGraph<Node1> const& wg, Node2 root);

    template <typename Node1, typename Node2>
    [[nodiscard]] Forest spanning_tree_no_checks(WordGraph<Node1> const& wg,
                                                 Node2                   root);
  }  // namespace word_graph

  //////////////////////////////////////////////////////////////////////////
  // WordGraph - non-member functions
  //////////////////////////////////////////////////////////////////////////

  //! \ingroup word_graph_group
  //! Output the edges of a wordGraph to a stream.
  //!
  //! This function outputs the word graph \p wg to the stream \p os.
  //! The word graph is represented by the out-neighbours of each node ordered
  //! according to their labels. The symbol `-` is used to denote that an
  //! edge is not defined. For example, the word graph with 1 nodes,
  //! out-degree 2, and a single loop labelled 1 from node 0 to 0 is
  //! represented as
  //! `{{-, 0}}`.
  //!
  //! \param os the ostream
  //! \param wg the word graph
  //!
  //! \returns
  //! The first parameter \p os.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  template <typename Node>
  std::ostream& operator<<(std::ostream& os, WordGraph<Node> const& wg);

  //! \ingroup word_graph_group
  //! Constructs a word graph from a number of nodes and an \c
  //! initializer_list.
  //!
  //! This function constructs a wordGraph from its arguments whose
  //! out-degree is specified by the length of the first \c initializer_list
  //! in the 2nd parameter.
  //!
  //! \tparam Node the type of the nodes of the word graph
  //!
  //! \param num_nodes_reachable_from_root the number of nodes in the word
  //! graph. \param il the out-targets of the word graph.
  //!
  //! \returns A value of type WordGraph.
  //!
  //! \throws LibsemigroupsException
  //! if WordGraph<Node>::target throws when adding edges from \p il.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the length of \p il and \f$n\f$ is the
  //! parameter \p num_nodes_reachable_from_root.
  //!
  //! \par Example
  //! \code
  //! // Construct a word graph with 5 nodes and 10 edges (7 specified)
  //! to_word_graph<uint8_t>(5, {{0, 0}, {1, 1}, {2}, {3, 3}});
  //! \endcode
  template <typename Node>
  [[nodiscard]] WordGraph<Node>
  to_word_graph(size_t num_nodes_reachable_from_root,
                std::initializer_list<std::vector<Node>> il);

  //! \ingroup word_graph_group
  //! TODO(doc)
  template <typename Node>
  [[nodiscard]] WordGraph<Node>
  to_word_graph(size_t num_nodes_reachable_from_root,
                std::vector<std::vector<Node>> const& v);

  namespace detail {
    template <typename Subclass>
    class JoinerMeeterCommon {
     private:
      template <typename Node1, typename Node2>
      void throw_if_bad_args(WordGraph<Node1> const& x,
                             Node2                   xroot,
                             WordGraph<Node1> const& y,
                             Node2                   yroot);

     public:
      template <typename Node>
      void call_no_checks(WordGraph<Node>&       xy,
                          WordGraph<Node> const& x,
                          Node                   xroot,
                          WordGraph<Node> const& y,
                          Node                   yroot);

      template <typename Node>
      void call_no_checks(WordGraph<Node>&       xy,
                          WordGraph<Node> const& x,
                          WordGraph<Node> const& y) {
        return call_no_checks(
            xy, x, static_cast<Node>(0), y, static_cast<Node>(0));
      }

      template <typename Node, typename... Args>
      [[nodiscard]] auto call_no_checks(WordGraph<Node> const& x,
                                        Args&&... args)
          -> std::enable_if_t<sizeof...(Args) % 2 == 1, WordGraph<Node>> {
        // The versions of this function changing the 1st argument in-place
        // always have an odd number of arguments, so we check that it's even
        // here (the argument x and an odd number of further arguments).
        WordGraph<Node> xy;
        static_cast<Subclass&>(*this).call_no_checks(
            xy, x, std::forward<Args>(args)...);
        return xy;
      }

      // There's no operator() with the number of nodes reachable from the
      // roots as arguments (7 args in total) because we'd have to check that
      // they were valid, and the only way to do this is to recompute them.

      // TODO(doc)
      template <typename Node>
      void operator()(WordGraph<Node>&       xy,
                      WordGraph<Node> const& x,
                      Node                   xroot,
                      WordGraph<Node> const& y,
                      Node                   yroot) {
        throw_if_bad_args(x, xroot, y, yroot);
        call_no_checks(xy, x, xroot, y, yroot);
      }

      template <typename Node>
      void operator()(WordGraph<Node>&       xy,
                      WordGraph<Node> const& x,
                      WordGraph<Node> const& y) {
        return operator()(xy, x, static_cast<Node>(0), y, static_cast<Node>(0));
      }

      template <typename Node, typename... Args>
      [[nodiscard]] auto operator()(WordGraph<Node> const& x, Args&&... args)
          -> std::enable_if_t<sizeof...(Args) % 2 == 1, WordGraph<Node>> {
        // The versions of this function changing the 1st argument in-place
        // always have an odd number of arguments, so we check that it's even
        // here (the argument x and an odd number of further arguments).
        WordGraph<Node> xy;
        operator()(xy, x, std::forward<Args>(args)...);
        return xy;
      }

      template <typename Node1, typename Node2>
      bool is_subrelation_no_checks(WordGraph<Node1> const& x,
                                    Node2                   xroot,
                                    WordGraph<Node1> const& y,
                                    Node2                   yroot);

      template <typename Node>
      bool is_subrelation_no_checks(WordGraph<Node> const& x,
                                    WordGraph<Node> const& y) {
        return is_subrelation_no_checks(
            x, static_cast<Node>(0), y, static_cast<Node>(0));
      }

      // There's no subrelation with the number of nodes reachable from the
      // roots as arguments (6 args in total) because we'd have to check that
      // they were valid, and the only way to do this is to recompute them.

      template <typename Node1, typename Node2>
      bool is_subrelation(WordGraph<Node1> const& x,
                          Node2                   xroot,
                          WordGraph<Node1> const& y,
                          Node2                   yroot) {
        throw_if_bad_args(x, xroot, y, yroot);
        return is_subrelation_no_checks(x, xroot, y, yroot);
      }

      template <typename Node>
      bool is_subrelation(WordGraph<Node> const& x, WordGraph<Node> const& y) {
        return is_subrelation(x, static_cast<Node>(0), y, static_cast<Node>(0));
      }
    };  // JoinerMeeterCommon
  }  // namespace detail

  //! \ingroup word_graph_group
  //! TODO(doc)
  // This class is intentionally not a template so that we don't have to
  // specify the types of the nodes when constructing one of these objects.
  // Instead every member function has a template parameter Node, which is
  // deduced from the argument.
  class Joiner : public detail::JoinerMeeterCommon<Joiner> {
   private:
    detail::Duf<>                             _uf;
    std::stack<std::pair<uint64_t, uint64_t>> _stck;

    template <typename Node>
    [[nodiscard]] Node find(WordGraph<Node> const& x,
                            size_t xnum_nodes_reachable_from_root,
                            WordGraph<Node> const&               y,
                            uint64_t                             n,
                            typename WordGraph<Node>::label_type a) const;

    template <typename Node>
    void run(WordGraph<Node> const& x,
             size_t                 xnum_nodes_reachable_from_root,
             Node                   xroot,
             WordGraph<Node> const& y,
             size_t                 ynum_nodes_reachable_from_root,
             Node                   yroot);

   public:
    //! \brief Default constructor.
    //!
    //! Default constructor.
    Joiner() = default;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    Joiner(Joiner const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    Joiner(Joiner&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    Joiner& operator=(Joiner const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    Joiner& operator=(Joiner&&) = default;

    ~Joiner() = default;

    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        size_t                 xnum_nodes_reachable_from_root,
                        Node                   xroot,
                        WordGraph<Node> const& y,
                        size_t                 ynum_nodes_reachable_from_root,
                        Node                   yroot);

    // is x a subrelation of y?
    // TODO(doc)
    template <typename Node1, typename Node2>
    [[nodiscard]] bool
    is_subrelation_no_checks(WordGraph<Node1> const& x,
                             size_t xnum_nodes_reachable_from_root,
                             Node2  xroot,
                             WordGraph<Node1> const& y,
                             size_t ynum_nodes_reachable_from_root,
                             Node2  yroot);

    using detail::JoinerMeeterCommon<Joiner>::call_no_checks;
    using detail::JoinerMeeterCommon<Joiner>::operator();
    using detail::JoinerMeeterCommon<Joiner>::is_subrelation_no_checks;
    using detail::JoinerMeeterCommon<Joiner>::is_subrelation;
  };  // Joiner

  //! \ingroup word_graph_group
  // Class for forming the meet of two word graphs
  class Meeter : public detail::JoinerMeeterCommon<Meeter> {
   private:
    using node_type = std::pair<uint64_t, uint64_t>;

    std::unordered_map<node_type, uint64_t, Hash<node_type>> _lookup;
    std::vector<node_type>                                   _todo;
    std::vector<node_type>                                   _todo_new;

   public:
    //! \brief Default constructor.
    //!
    //! Default constructor.
    Meeter() = default;

    //! \brief Default copy constructor.
    //!
    //! Default copy constructor.
    Meeter(Meeter const&) = default;

    //! \brief Default move constructor.
    //!
    //! Default move constructor.
    Meeter(Meeter&&) = default;

    //! \brief Default copy assignment operator.
    //!
    //! Default copy assignment operator.
    Meeter& operator=(Meeter const&) = default;

    //! \brief Default move assignment operator.
    //!
    //! Default move assignment operator.
    Meeter& operator=(Meeter&&) = default;

    ~Meeter() = default;

    // TODO(doc)
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        size_t                 xnum_nodes_reachable_from_root,
                        Node                   xroot,
                        WordGraph<Node> const& y,
                        size_t                 ynum_nodes_reachable_from_root,
                        Node                   yroot);

    // is x a subrelation of y
    template <typename Node1, typename Node2>
    [[nodiscard]] bool
    is_subrelation_no_checks(WordGraph<Node1> const& x,
                             size_t xnum_nodes_reachable_from_root,
                             Node2  xroot,
                             WordGraph<Node1> const& y,
                             size_t ynum_nodes_reachable_from_root,
                             Node2  yroot);

    using detail::JoinerMeeterCommon<Meeter>::call_no_checks;
    using detail::JoinerMeeterCommon<Meeter>::operator();
    using detail::JoinerMeeterCommon<Meeter>::is_subrelation_no_checks;
    using detail::JoinerMeeterCommon<Meeter>::is_subrelation;
  };  // class Meeter

  template <typename Node>
  [[nodiscard]] std::string to_human_readable_repr(WordGraph<Node> const& wg);

}  // namespace libsemigroups

#include "word-graph.tpp"

#endif  // LIBSEMIGROUPS_WORD_GRAPH_HPP_
