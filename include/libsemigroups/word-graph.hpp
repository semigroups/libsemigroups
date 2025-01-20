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
#include "types.hpp"      // for word_type, enable_if_is_same

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

  //! \defgroup word_graph_group Word graphs and related functionality
  //!
  //! This page collects various classes and functions related to word graphs
  //! in `libsemigroups`.
  //!
  //! \sa \ref gabow_group, \ref paths_group

  //! \ingroup word_graph_group
  //!
  //! \brief Class for representing word graphs.
  //!
  //! Defined in `word-graph.hpp`.
  //!
  //! Instances of this class represent word graphs. If the word graph has \c n
  //! nodes, they are represented by the numbers \f$\{0, ..., n - 1\}\f$, and
  //! every node has the same number \c m of out-edges (edges with source that
  //! node and target any other node). The number \c m is referred to as the
  //! *out-degree* of the word graph, or any of its nodes.
  //!
  //! \tparam Node the type of the nodes in the word graph, must be an unsigned
  //! integer type.
  template <typename Node>
  class WordGraph {
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

    //! The type of an iterator pointing to the targets of a node.
    using const_iterator_targets =
        typename detail::DynamicArray2<Node>::const_iterator;

    ////////////////////////////////////////////////////////////////////////
    // WordGraph - constructors + destructor - public
    ////////////////////////////////////////////////////////////////////////

    //! \brief Construct from number of nodes and out degree.
    //!
    //! This function constructs a word graph with \p m nodes and where the
    //! maximum out-degree of any node is \p n. There are no edges in the
    //! defined word graph.
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

    //! \brief Re-initialize the word graph to have \p m nodes and out-degree \p
    //! n.
    //!
    //! This function puts a word graph into the state that it would have been
    //! in if it had just been newly constructed with the same parameters \p m
    //! and \p n.
    //!
    //! \param m the number of nodes in the word graph.
    //! \param n the out-degree of every node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(mn)\f$ where \f$m\f$ is the number of nodes, and \f$n\f$ is the
    //! out-degree of the word graph.
    WordGraph& init(size_type m = 0, size_type n = 0);

    //! Default copy constructor
    WordGraph(WordGraph const&);

    //! \brief Construct from WordGraph with another node type.
    //!
    //! This function can be used to construct a WordGraph<Node> as a copy of a
    //! WordGraph<OtherNode> so long as `sizeof(OtherNode) <= sizeof(Node)`.
    //!
    //! \param that the word graph to copy.
    //!
    //! \note Any edge with target \ref UNDEFINED in \p that will have target
    //! `static_cast<Node>(UNDEFINED)` in the constructed word graph.
    template <typename OtherNode>
    WordGraph(WordGraph<OtherNode> const& that);

    //! \brief Re-initialize the word graph contain a copy of another.
    //!
    //!
    //! This function puts a word graph into the same state that it would have
    //! been in if it had just been newly copy constructed with the same
    //! parameter \p that.
    //!
    //! \param that the word graph to copy.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \par Complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes, and \p n is the
    //! out-degree of the word graph.
    //!
    //! \note Any edge with target \ref UNDEFINED in \p that will have target
    //! `static_cast<Node>(UNDEFINED)` in `*this`.
    template <typename OtherNode>
    WordGraph& init(WordGraph<OtherNode> const& that);

    //! Default move constructor
    WordGraph(WordGraph&&);

    //! Default copy assignment constructor
    WordGraph& operator=(WordGraph const&);

    //! Default move assignment constructor
    WordGraph& operator=(WordGraph&&);

    ~WordGraph();

    //! \brief Construct a random word graph from number of nodes and
    //! out-degree.
    //!
    //! This function constructs a random word graph with \p number_of_nodes
    //! nodes and out-degree \p out_degree.
    //!
    //! \param number_of_nodes the number of nodes.
    //! \param out_degree the out-degree of every node.
    //! \param mt a std::mt19937 used as a random source (defaults to:
    //! std::mt19937(std::random_device()())).
    //!
    //! \returns A random word graph.
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

    //! \brief Ensures that the word graph has capacity for a given number of
    //! nodes, and out-degree.
    //!
    //!  This function ensures that the word graph has capacity for \p m nodes
    //!  and \p n labels.
    //!
    //! \param m the number of nodes.
    //! \param n the out-degree.
    //!
    //! \returns A reference to `*this`.
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
    //!
    //! \note
    //! Does not modify number_of_nodes() or out_degree().
    // Not noexcept because DynamicArray2::add_cols isn't.
    WordGraph& reserve(size_type m, size_type n);

    ////////////////////////////////////////////////////////////////////////
    // WordGraph - modifiers - public
    ////////////////////////////////////////////////////////////////////////

    //! \brief Add a number of new nodes.
    //!
    //! This function modifies a word graph in-place so that it has \p nr new
    //! nodes added.
    //!
    //! \param nr the number of nodes to add.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //! \strong_guarantee
    //!
    //! \complexity
    //! Linear in `(number_of_nodes() + nr) * out_degree()`.
    //!
    //! \iterator_validity
    //! \iterator_invalid
    // Not noexcept because DynamicArray2::add_rows isn't.
    WordGraph& add_nodes(size_type nr);

    //! \brief Add to the out-degree of every node.
    //!
    //! This function modifies a word graph in-place so that the out-degree is
    //! increased by \p nr.
    //!
    //! \param nr the number of new out-edges for every node.
    //!
    //! \returns A reference to `*this`.
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

    //! \brief Add an edge from one node to another with a given label.
    //!
    //! If \p s and \p t are nodes in \c this, and \p a is in the range `[0,
    //! out_degree())`, then this function adds an edge from \p a to \p b
    //! labelled \p a.
    //!
    //! \param s the source node.
    //! \param a the label of the edge.
    //! \param t the range node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p s, \p a, or \p t is not valid.
    //! \strong_guarantee
    //!
    //! \complexity
    //! Constant.
    // Not noexcept because throw_if_node_out_of_bounds/label aren't
    // return *this by reference.
    WordGraph& target(node_type s, label_type a, node_type t);

    //! \brief Add an edge from one node to another with a given label.
    //!
    //! This function adds an edge from the node \p s to the node \p t with
    //! label \p a.
    //!
    //! \param s the source node.
    //! \param a the label of the edge from \p s to \p t.
    //! \param t the target node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    inline WordGraph& target_no_checks(node_type s, label_type a, node_type t) {
      _dynamic_array_2.set(s, a, t);
      return *this;
    }

    //! \brief Remove an edge from a node with a given label.
    //!
    //! This function removes the edge with source node \p s labelled by \p a.
    //!
    //! \param s the source node.
    //! \param a the label of the edge from \p s.
    //!
    //! \returns A reference to `*this`.
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

    //! \brief Remove an edge from a node with a given label.
    //!
    //! This function removes the edge with source node \p s labelled by \p a.
    //!
    //! \param s the source node.
    //! \param a the label of the edge from \p s.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p s or \p a is out of range.
    //!
    //! \complexity
    //! Constant.
    WordGraph& remove_target(node_type s, label_type a);

    //! \brief Removes a given label from the word graph.
    //!
    //! This function removes the label \p a from a WordGraph
    //! object in-place. This reduces the out-degree by \c 1.
    //!
    //! \param a the label to remove.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p a is out of range.
    WordGraph& remove_label(label_type a);

    //! \brief Removes a given label from the word graph.
    //!
    //! This function removes the label \p a from a WordGraph
    //! object in-place. This reduces the out-degree by \c 1.
    //!
    //! \param a the label to remove.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning No checks are performed on the argument, and, in particular,
    //! it assumed that \p a is not out of range.
    WordGraph& remove_label_no_checks(label_type a);

    //! \brief Remove all of the edges in the word graph.
    //!
    //! Set every target of every source with every possible label to \ref
    //! UNDEFINED.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \p m is the number of nodes and \p n is the
    //! out-degree.
    inline WordGraph& remove_all_targets() {
      std::fill(_dynamic_array_2.begin(), _dynamic_array_2.end(), UNDEFINED);
      return *this;
    }

    //! \brief Swaps the edge with specified label from one node with another.
    //!
    //! This function swaps the target of the edge from the node \p m labelled
    //! \p a with the target of the edge from the node \p n labelled \p a.
    //!
    //! \param m the first node.
    //! \param n the second node.
    //! \param a the label.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    // swap m - a - > m' and n - a -> n'
    WordGraph& swap_targets_no_checks(node_type m, node_type n, label_type a) {
      _dynamic_array_2.swap(m, a, n, a);
      return *this;
    }

    //! \brief Swaps the edge with specified label from one node with another.
    //!
    //! This function swaps the target of the edge from the node \p m labelled
    //! \p a with the target of the edge from the node \p n labelled \p a.
    //!
    //! \param m the first node.
    //! \param n the second node.
    //! \param a the label.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p m, \p n, or \p a are out of range.
    //!
    //! \complexity
    //! Constant
    WordGraph& swap_targets(node_type m, node_type n, label_type a);

    //! \brief Check if two word graphs are equal.
    //!
    //! This function returns \c true if `*this` and \p that are equal, and \c
    //! false if not.
    //!
    //! \param that the word graph for comparison.
    //!
    //! \returns
    //! Whether or not `*this` and \c that are equal.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    [[nodiscard]] bool operator==(WordGraph const& that) const {
      return _dynamic_array_2 == that._dynamic_array_2;
    }

    //! \brief Check if two word graphs are inequal.
    //!
    //! This function returns \c true if `*this` and \p that are not equal, and
    //! \c false if they are equal.
    //!
    //! \param that the word graph for comparison.
    //!
    //! \returns
    //! Whether or not `*this` and \c that are not equal.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    [[nodiscard]] bool operator!=(WordGraph const& that) const {
      return !operator==(that);
    }

    //! \brief Check if a word graph is less than another.
    //!
    //! This function returns \c true if `*this` is less than \p that. This
    //! operator defines a linear order on word graphs.
    //!
    //! \param that the word graph for comparison.
    //!
    //! \returns
    //! Whether or not `*this` is strictly less than \c that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    [[nodiscard]] bool operator<(WordGraph const& that) const {
      return _dynamic_array_2 < that._dynamic_array_2;
    }

    //! \brief Check if a word graph is less than or equal to another.
    //!
    //! This function returns \c true if `*this` is less or equal to \p that.
    //! This operator defines a linear order on word graphs.
    //!
    //! \param that the word graph for comparison.
    //!
    //! \returns
    //! Whether or not `*this` is less than or equal to \c that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    [[nodiscard]] bool operator<=(WordGraph const& that) const {
      return _dynamic_array_2 <= that._dynamic_array_2;
    }

    //! \brief Check if a word graph is greater than another.
    //!
    //! This function returns \c true if `*this` is greater than \p that.
    //! This operator defines a linear order on word graphs.
    //!
    //! \param that the word graph for comparison.
    //!
    //! \returns
    //! Whether or not `*this` is strictly greater than \c that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    [[nodiscard]] bool operator>(WordGraph const& that) const {
      return _dynamic_array_2 > that._dynamic_array_2;
    }

    //! \brief Check if a word graph is greater than or equal to another.
    //!
    //! This function returns \c true if `*this` is greater or equal to \p that.
    //! This operator defines a linear order on word graphs.
    //!
    //! \param that the word graph for comparison.
    //!
    //! \returns
    //! Whether or not `*this` is greater than or equal to \c that.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(nm)\f$ where \f$n\f$ is the number of nodes and \f$m\f$
    //! is the out-degree of the word graph.
    [[nodiscard]] bool operator>=(WordGraph const& that) const {
      return _dynamic_array_2 > that._dynamic_array_2;
    }

    //! \brief Returns a hash value.
    //!
    //! \returns
    //! A value of type \c size_t.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Linear in \ref number_of_nodes times \ref out_degree.
    // not noexcept because Hash<T>::operator() isn't
    [[nodiscard]] size_t hash_value() const {
      return std::hash<decltype(_dynamic_array_2)>()(_dynamic_array_2);
    }

    ////////////////////////////////////////////////////////////////////////
    // WordGraph - nodes, targets, etc - public
    ////////////////////////////////////////////////////////////////////////

    //! \brief Get the target of the edge with given source node and label.
    //!
    //! This function returns the target of the edge with source node \p source
    //! and label \p a.
    //!
    //! \param source the node.
    //! \param a the label.
    //!
    //! \returns
    //! Returns the node adjacent to \p source via the edge labelled \p a, or
    //! \ref UNDEFINED; both are values of type \ref node_type.
    //!
    //! \complexity
    //! Constant.
    //!
    //! \throws LibsemigroupsException if \p source or \p a is not
    //! valid.
    // Not noexcept because throw_if_node_out_of_bounds/label aren't
    [[nodiscard]] node_type target(node_type source, label_type a) const;

    //! \brief Get the target of the edge with given source node and label.
    //!
    //! This function returns the the target of the edge with source node \p s
    //! and label \p a.
    //!
    //! \param s the node.
    //! \param a the label.
    //!
    //! \returns
    //! Returns the node adjacent to \p s via the edge labelled \p a, or
    //! \ref UNDEFINED; both are values of type \ref node_type.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! Constant.
    //!
    //! \warning This function does not verify \p s or \p a is valid.
    // Not noexcept because DynamicArray2::get is not
    [[nodiscard]] node_type inline target_no_checks(node_type  s,
                                                    label_type a) const {
      return _dynamic_array_2.get(s, a);
    }

    //! \brief Get the next target of an edge incident to a given node that
    //! doesn't equal \ref UNDEFINED.
    //!
    //! This function returns the next target of an edge with label greater
    //! than or equal to \p a that is incident to the node \p s.
    //!
    //! \param s the source node.
    //! \param a the label.
    //!
    //! \returns
    //! Returns a std::pair \c x where:
    //! 1. \c x.first is adjacent to \p s via an edge labelled \c x.second;
    //! 2. \c x.second is the minimum value in the range \f$[a, n)\f$ such that
    //!    `target(s, x.second)` is not equal to \ref UNDEFINED (where \f$n\f$
    //!    is the return value of \ref out_degree); and
    //! If no such value exists, then `{UNDEFINED, UNDEFINED}` is returned.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ equals out_degree().
    //!
    //! \warning This function does not check its arguments, in particular it is
    //! not verified that the parameter \p s represents a node of \c this, or
    //! that \p a is a valid label.
    //!
    //! \sa next_label_and_target.
    // Not noexcept because DynamicArray2::get is not
    [[nodiscard]] std::pair<label_type, node_type>
    next_label_and_target_no_checks(node_type s, label_type a) const;

    //! \brief Get the next target of an edge incident to a given node that
    //! doesn't equal \ref UNDEFINED.
    //!
    //! This function returns the next target of an edge with label greater
    //! than or equal to \p a that is incident to the node \p s.
    //!
    //! If `target(s, b)` equals \ref UNDEFINED for every value \c b in the
    //! range \f$[a, n)\f$, where \f$n\f$ is the return value of out_degree()
    //! then \c x.first and \c x.second equal \ref UNDEFINED.
    //!
    //! \param s the node.
    //! \param a the label.
    //!
    //! \returns
    //! Returns a std::pair
    //! \c x where:
    //! 1. \c x.first is adjacent to \p s via an edge labelled
    //!    \c x.second; and
    //! 2. \c x.second is the minimum value in the range \f$[a, n)\f$ such that
    //!    `target(s, x.second)` is not equal to \ref UNDEFINED
    //!    where \f$n\f$ is the return value of out_degree();
    //! If no such value exists, then `{UNDEFINED, UNDEFINED}` is returned.
    //!
    //! \throws LibsemigroupsException if \p s does not represent a node in \c
    //! this, or \p a is not a valid edge label.
    //!
    //! \complexity
    //! At worst \f$O(n)\f$ where \f$n\f$ equals out_degree().
    //!
    //! \sa next_label_and_target_no_checks.
    // Not noexcept because next_label_and_target_no_checks is not
    [[nodiscard]] std::pair<label_type, node_type>
    next_label_and_target(node_type s, label_type a) const;

    //! \brief Returns the number of nodes.
    //!
    //! This function returns the number of nodes in the word graph.
    //!
    //! \returns
    //! The number of nodes in the word graph.
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

    //! \brief Returns the number of edges.
    //!
    //! This function returns the total number of edges (i.e. values \c s and \c
    //! a such that `target(s, a)` is not UNDEFINED) in the word graph.
    //!
    //! \returns
    //! The total number of edges, a value of type \c size_type.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is number_of_nodes() and \c n is out_degree().
    // Not noexcept because std::count isn't
    [[nodiscard]] size_type number_of_edges() const;

    //! \brief Returns the number of edges with given source node.
    //!
    //! This function returns the number of edges incident to the given source
    //! node \p s.
    //!
    //! \param s the node.
    //!
    //! \returns
    //! A value of type \c size_type.
    //!
    //! \throws LibsemigroupsException if \p s is not a node.
    //!
    //! \complexity
    //! \f$O(n)\f$ where \c n is out_degree().
    [[nodiscard]] size_type number_of_edges(node_type s) const;

    //! \brief Returns the number of edges with given source node.
    //!
    //! This function returns the number of edges incident to the given source
    //! node \p s.
    //!
    //! \param s the node.
    //!
    //! \returns
    //! A value of type \c size_type.
    //!
    //! \complexity
    //! \f$O(n)\f$ where \c n is out_degree().
    //!
    //! \warning No checks are performed that the argument \p s is actually a
    //! node in the word graph.
    [[nodiscard]] size_type number_of_edges_no_checks(node_type s) const;

    //! \brief Returns the out-degree.
    //!
    //! This function returns the number of edge labels in the word graph.
    //!
    //! \returns
    //! The number of edge labels, a value of type \c size_type.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] size_type out_degree() const noexcept {
      return _degree;
    }

    //! \brief Returns a random access iterator pointing at the first node of
    //! the word graph.
    //!
    //! This function returns a random access iterator pointing at the first
    //! node on the word graph.
    //!
    //! \returns
    //! An \ref const_iterator_nodes.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    const_iterator_nodes cbegin_nodes() const noexcept {
      return detail::IntRange<node_type>(0, number_of_nodes()).cbegin();
    }

    //! \brief Returns a random access iterator pointing one-past-the-last node
    //! of the word graph.
    //!
    //! This function returns a random access iterator pointing one beyond the
    //! last node in the word graph.
    //!
    //! \returns
    //! An \ref const_iterator_nodes.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] const_iterator_nodes cend_nodes() const noexcept {
      return detail::IntRange<node_type>(0, number_of_nodes()).cend();
    }

    //! \brief Returns a random access iterator pointing at the target of
    //! the edge with label \p 0 incident to a given node.
    //!
    //! This function returns a random access iterator pointing at the
    //! target of the edge with label \c 0 incident to the source node \p
    //! source. This target might equal \ref UNDEFINED.
    //!
    //! \param source the source node in the word graph.
    //!
    //! \returns
    //! A \ref const_iterator_targets.
    //!
    //! \throws LibsemigroupsException if \p source is out of range (i.e.
    //! greater than or equal to \ref number_of_nodes).
    //!
    //! \complexity
    //! Constant.
    // Not noexcept because throw_if_node_out_of_bounds isn't
    [[nodiscard]] const_iterator_targets cbegin_targets(node_type source) const;

    //! \brief Returns a random access iterator pointing at the target of
    //! the edge with label \p 0 incident to a given node.
    //!
    //! This function returns a random access iterator pointing at the
    //! target of the edge with label \c 0 incident to the source node \p
    //! source. This target might equal \ref UNDEFINED.
    //!
    //! \param source a node in the word graph.
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
    cbegin_targets_no_checks(node_type source) const noexcept {
      return _dynamic_array_2.cbegin_row(source);
    }

    //! \brief Returns a random access iterator pointing one beyond the target
    //! of the edge with label `out_degree() - 1` incident to a given node.
    //!
    //! This function returns a random access iterator pointing one beyond the
    //! target of the edge with label `out_degree() - 1` incident to the source
    //! node \p source. This target might equal \ref UNDEFINED.
    //!
    //! \param source a node in the word graph.
    //!
    //! \returns
    //! A \ref const_iterator_targets.
    //!
    //! \throws LibsemigroupsException if \p source is out of range (i.e.
    //! greater than or equal to \ref number_of_nodes).
    //!
    //! \complexity
    //! Constant.
    // Not noexcept because throw_if_node_out_of_bounds isn't
    [[nodiscard]] const_iterator_targets cend_targets(node_type source) const;

    //! \brief Returns a random access iterator pointing one beyond the target
    //! of the edge with label `out_degree() - 1` incident to a given node.
    //!
    //! This function returns a random access iterator pointing one beyond the
    //! target of the edge with label `out_degree() - 1` incident to the source
    //! node \p source. This target might equal \ref UNDEFINED.
    //!
    //! \param source a node in the word graph.
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
    cend_targets_no_checks(node_type source) const noexcept {
      return _dynamic_array_2.cbegin_row(source) + _degree;
    }

    //! \brief Returns a range object containing all nodes in a word graph.
    //!
    //! This function returns a range object containing all the nodes in a
    //! word graph.
    //!
    //! \returns A range object.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] auto nodes() const noexcept {
      return rx::seq<node_type>() | rx::take(number_of_nodes());
    }

    //! \brief Returns a range object containing all labels of edges in a word
    //! graph.
    //!
    //! This function returns a range object containing all the labels of edges
    //! in a word graph.
    //!
    //! \returns A range object.
    //!
    //! \exceptions
    //! \noexcept
    [[nodiscard]] auto labels() const noexcept {
      return rx::seq<label_type>() | rx::take(out_degree());
    }

    //! \brief Returns a range object containing all the targets of edges with
    //! a given source.
    //!
    //! This function returns a range object containing all the targets of
    //! edges with source \p source.
    //!
    //! \param source the source node.
    //!
    //! \returns A range object.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \warning
    //! This function performs no checks whatsoever and assumes that \p source
    //! is a valid node of the word graph (i.e. it is not greater than or equal
    //! to \ref number_of_nodes).
    [[nodiscard]] rx::iterator_range<const_iterator_targets>
    targets_no_checks(node_type source) const noexcept {
      return rx::iterator_range(cbegin_targets_no_checks(source),
                                cend_targets_no_checks(source));
    }

    //! \brief Returns a range object containing all the targets of edges with
    //! a given source.
    //!
    //! This function returns a range object containing all the targets of
    //! edges with source \p source.
    //!
    //! \param source the source node.
    //!
    //! \returns A range object.
    //!
    //! \throws LibsemigroupsException if \p source is out of range (i.e. it is
    //! greater than or equal to \ref number_of_nodes)
    [[nodiscard]] rx::iterator_range<const_iterator_targets>
    targets(node_type source) const;

    //! \brief Returns a range object containing pairs consisting of edge
    //! labels and target nodes.
    //!
    //! This function returns a range object containing all the edge labels and
    //! targets of edges with source \p source.
    //!
    //! \param source the source node.
    //!
    //! \returns A range object.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \warning
    //! This function performs no checks whatsoever and assumes that \p source
    //! is a valid node of the word graph (i.e. it is not greater than or equal
    //! to \ref number_of_nodes).
    [[nodiscard]] auto
    labels_and_targets_no_checks(node_type source) const noexcept {
      return rx::enumerate(targets_no_checks(source));
    }

    //! \brief Returns a range object containing pairs consisting of edge
    //! labels and target nodes.
    //!
    //! This function returns a range object containing all the edge labels and
    //! targets of edges with source \p source.
    //!
    //! \param source the source node.
    //!
    //! \returns A range object.
    //!
    //! \throws LibsemigroupsException if \p source is out of bounds.
    [[nodiscard]] auto labels_and_targets(node_type source) const;

    //! \brief Modify in-place to contain the subgraph induced by a range of
    //! nodes.
    //!
    //! This function modifies a WordGraph object in-place to contain its
    //! subgraph induced by the range of nodes \p first to \p last.
    //!
    //! \param first the first node.
    //! \param last one more than the last node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! This function performs no checks whatsoever and will result in a
    //! corrupted word graph if there are any edges from the nodes \f$0, \ldots,
    //! n - 1\f$ to nodes larger than \f$n - 1\f$.
    WordGraph& induced_subgraph_no_checks(node_type first, node_type last);

    //! \brief Modify in-place to contain the subgraph induced by a range of
    //! nodes.
    //!
    //! This function modifies a WordGraph object in-place to contain its
    //! subgraph induced by the range of nodes \p first to \p last.
    //!
    //! \param first the first node.
    //! \param last one more than the last node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p first or \p last is out of range.
    //! \throws LibsemigroupsException if any edge with source in the range \p
    //! first to \p last has target outside the range \p first to \p last.
    WordGraph& induced_subgraph(node_type first, node_type last);

    //! \brief Modify in-place to contain the subgraph induced by a range of
    //! nodes.
    //!
    //! This function modifies a WordGraph object in-place to contain its
    //! subgraph induced by the range of nodes \p first to \p last.
    //!
    //! \tparam Iterator the type of \p first and \p last (should be
    //! iterators).
    //!
    //! \param first the first node.
    //!
    //! \param last one more than the last node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning
    //! This function performs no checks whatsoever and will result in a
    //! corrupted word graph if there are any edges from the nodes \f$0, \ldots,
    //! n - 1\f$ to nodes larger than \f$n - 1\f$.
    template <typename Iterator,
              typename = std::enable_if_t<detail::IsIterator<Iterator>::value>>
    WordGraph& induced_subgraph_no_checks(Iterator first, Iterator last);

    //! \brief Modify in-place to contain the subgraph induced by a range of
    //! nodes.
    //!
    //! This function modifies a WordGraph object in-place to contain its
    //! subgraph induced by the range of nodes \p first to \p last.
    //!
    //! \tparam Iterator the type of \p first and \p last (should be
    //! iterators).
    //!
    //! \param first iterator pointing at the first node.
    //!
    //! \param last iterator pointing one beyond the last node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if any value in the range specified by \p
    //! first and \p last is not a node of the word graph.
    //!
    //! \throws LibsemigroupsException if any target of any edge with
    //! source node in the range specified by \p first and \p last does not
    //! belong to the same range.
    template <typename Iterator,
              typename = std::enable_if_t<detail::IsIterator<Iterator>::value>>
    WordGraph& induced_subgraph(Iterator first, Iterator last);

    //! \brief Unites a word graph in-place.
    //!
    //! This function changes a WordGraph object in-place to contain the
    //! disjoint union of itself and \p that. The node \c n of \p that is mapped
    //! to `this->number_of_nodes() + n`.
    //!
    //! \param that the word graph to unite.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \warning This function does not check its arguments, and it is assumed
    //! that `*this` and \p that have equal out degree.
    WordGraph& disjoint_union_inplace_no_checks(WordGraph<Node> const& that);

    //! \brief Unites a word graph in-place.
    //!
    //! This function changes a WordGraph object in-place to contain the
    //! disjoint union of itself and \p that. The node \c n of \p that is mapped
    //! to `this->number_of_nodes() + n`.
    //!
    //! \param that the word graph to unite.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if `*this` and `that` do not have the
    //! same out-degree.
    WordGraph& disjoint_union_inplace(WordGraph<Node> const& that);

#ifndef PARSED_BY_DOXYGEN
    // These are currently undocumented, because they are hard to use correctly,
    // shouldn't p and q be actual permutation objects?

    // requires access to apply_row_permutation so can't be helper
    WordGraph& permute_nodes_no_checks(std::vector<node_type> const& p,
                                       std::vector<node_type> const& q,
                                       size_t                        m);

    WordGraph& permute_nodes_no_checks(std::vector<node_type> const& p,
                                       std::vector<node_type> const& q) {
      return permute_nodes_no_checks(p, q, p.size());
    }
#endif

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
    // TODO(1) remove when WordGraphView is implemented
    size_type                           _num_active_nodes;
    mutable detail::DynamicArray2<Node> _dynamic_array_2;
  };

  //! \ingroup word_graph_group
  //!
  //! \brief Namespace containing helper functions for the \ref WordGraph
  //! class.
  //!
  //! Defined in `word-graph.hpp`.
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
    //! \tparam Node  the type of the nodes of the WordGraph.
    //!
    //! \tparam Iterator the type of an iterator pointing to nodes of a word
    //! graph.
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
    //! \tparam Node  the type of the nodes of the WordGraph.
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
    //! type of the returned matrix depends on whether or not `libsemigroups`
    //! is compiled with [eigen][] enabled. The returned matrix has the number
    //! of edges with source \c s and target \c t in the `(s, t)`-entry.
    //!
    //! \tparam Node  the type of the nodes of the WordGraph.
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
    //! \tparam Node  the type of the nodes of the WordGraph.
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
    //! \tparam Node  the type of the nodes of the WordGraph.
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
    //! \tparam Node  the type of the nodes of the WordGraph.
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
    //! \tparam Node1  the type of the nodes of the WordGraph
    //! \p wg.
    //!
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \tparam Iterator the type of \p first and \p last.
    //!
    //! \param wg a word graph.
    //! \param source the starting node.
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
                                    Node2                   source,
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
    //! \tparam Node1  the type of the nodes of the WordGraph
    //! \p wg.
    //!
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
    // TODO(2) example
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
    //! \tparam Node1  the type of the nodes of the WordGraph
    //! \p wg.
    //!
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \param wg a word graph.
    //! \param from the source node.
    //! \param first iterator into a word.
    //! \param last iterator into a word.
    //!
    //! \exceptions
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
    //! \tparam Node1  the type of the nodes of the WordGraph
    //! \p wg.
    //!
    //! \tparam Node2 the type of the node \p from (must satisfy `sizeof(Node2)
    //! <= sizeof(Node1)`).
    //!
    //! \param wg a word graph.
    //! \param from the source node.
    //! \param path the word.
    //!
    //! \exceptions
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
    //! \tparam Node  the type of the nodes of the WordGraph.
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
    //! \tparam Node1  the type of the nodes of the WordGraph
    //! \p wg.
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
    //! \tparam Node1  the type of the nodes of the WordGraph
    //! \p wg.
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
    //! \tparam Node  the type of the nodes of the WordGraph
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
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note This function ignores out of bound targets in \p wg (if any).
    //!
    //! \warning
    //! No checks on the arguments of this function are performed.
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
    //! \tparam Node  the type of the nodes of the WordGraph
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
    //! \throws LibsemigroupsException if any of the nodes in the range between
    //! \p first_node and \p last_node does not belong to \p wg (i.e. is greater
    //! than or equal to WordGraph::number_of_nodes).
    //!
    //! \throws LibsemigroupsException if any of the rules in the range between
    //! \p first_rule and \p last_rule contains an invalid label (i.e. one
    //! greater than or equal to WordGraph::out_degree).
    //!
    //! \note This function ignores out of bound targets in \p wg (if any).
    template <
        typename Node,
        typename Iterator1,
        typename Iterator2,
        typename Iterator3,
        std::enable_if_t<!std::is_same_v<std::decay_t<Iterator3>, word_type>>>
    [[nodiscard]] bool is_compatible(WordGraph<Node> const& wg,
                                     Iterator1              first_node,
                                     Iterator2              last_node,
                                     Iterator3              first_rule,
                                     Iterator3              last_rule);

    //! \brief Check if a word graph is compatible with a pair of words for a
    //! range of nodes.
    //!
    //! This function returns \c true if the word graph \p wg is compatible
    //! with the words \p lhs and \p rhs at every node in the range from \p
    //! first_node to \p last_node. This means that the paths with given sources
    //! that are labelled by \p lhs leads to the same node as the path labelled
    //! by \p rhs.
    //!
    //! \tparam Node the type of the nodes of the WordGraph \p wg.
    //!
    //! \tparam Iterator1 the type of \p first_node.
    //!
    //! \tparam Iterator2 the type of \p last_node.
    //!
    //! \param wg the word graph.
    //!
    //! \param first_node iterator pointing at the first node.
    //!
    //! \param last_node iterator pointing at one beyond the last node.
    //!
    //! \param lhs the first rule.
    //!
    //! \param rhs the second rule.
    //!
    //! \return Whether or not the word graph is compatible with the given rules
    //! at each one of the given nodes.
    //!
    //! \note This function ignores out of bound targets in \p wg (if any).
    //!
    //! \warning This function does not check that its arguments are valid.
    template <typename Node, typename Iterator1, typename Iterator2>
    bool is_compatible_no_checks(WordGraph<Node> const& wg,
                                 Iterator1              first_node,
                                 Iterator2              last_node,
                                 word_type const&       lhs,
                                 word_type const&       rhs);

    //! \brief Check if a word graph is compatible with a pair of words for a
    //! range of nodes.
    //!
    //! This function returns \c true if the word graph \p wg is compatible
    //! with the words \p lhs and \p rhs at every node in the range from \p
    //! first_node to \p last_node. This means that the paths with given sources
    //! that are labelled by \p lhs leads to the same node as the path labelled
    //! by \p rhs.
    //!
    //! \tparam Node the type of the nodes of the WordGraph \p wg.
    //!
    //! \tparam Iterator1 the type of \p first_node.
    //!
    //! \tparam Iterator2 the type of \p last_node.
    //!
    //! \param wg the word graph.
    //!
    //! \param first_node iterator pointing at the first node.
    //!
    //! \param last_node iterator pointing at one beyond the last node.
    //!
    //! \param lhs the first rule.
    //!
    //! \param rhs the second rule.
    //!
    //! \return Whether or not the word graph is compatible with the given rules
    //! at each one of the given nodes.
    //!
    //! \throws LibsemigroupsException if any of the nodes in the range between
    //! \p first_node and \p last_node does not belong to \p wg (i.e. is greater
    //! than or equal to WordGraph::number_of_nodes).
    //!
    //! \throws LibsemigroupsException if any of the rules in the range between
    //! \p first_rule and \p last_rule contains an invalid label (i.e. one
    //! greater than or equal to WordGraph::out_degree).
    template <typename Node, typename Iterator1, typename Iterator2>
    bool is_compatible(WordGraph<Node> const& wg,
                       Iterator1              first_node,
                       Iterator2              last_node,
                       word_type const&       lhs,
                       word_type const&       rhs);

    //! \brief Check if every node in a range has exactly WordGraph::out_degree
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

    //! \brief Check if every node in a range has exactly WordGraph::out_degree
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
    //! \throws LibsemigroupsException if any item in the range defined by \p
    //! first_node and \p last_node is not a node of \p wg.
    //!
    //! \complexity
    //! \f$O(mn)\f$ where \c m is the number of nodes in the range and \c n is
    //! out_degree().
    template <typename Node, typename Iterator1, typename Iterator2>
    [[nodiscard]] bool is_complete(WordGraph<Node> const& wg,
                                   Iterator1              first_node,
                                   Iterator2              last_node);

    //! \brief Check if every node has exactly WordGraph::out_degree out-edges.
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
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
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
    //! \note This function ignores out of bound targets in \p wg (if any).
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
    //! where \f$k\f$ is the WordGraph::out_degree.
    //!
    //! \note
    //! If \p source and \p target are equal, then, by convention, we consider
    //! \p target to be reachable from \p source, via the empty path.
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
    //! \tparam Node the type of the nodes of the WordGraph.
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
    //! auto wg = make<WordGraph<uint8_t>>(
    //!     5, {{0, 0}, {1, 1}, {2}, {3, 3}});
    //! word_graph::is_strictly_cyclic(wg);  // returns false
    //! \endcode
    // TODO(1) should have a version that returns the node that everything is
    // reachable from
    template <typename Node>
    [[nodiscard]] bool is_strictly_cyclic(WordGraph<Node> const& wg);

    //! \brief Returns the last node on the path labelled by a word and an
    //! iterator to the position in the word reached.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //! \tparam Iterator the type of the iterators into a word.
    //!
    //! \param wg a word graph.
    //! \param source the source node.
    //! \param first iterator into a word.
    //! \param last iterator into a word.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \returns A pair consisting of the last node reached and an iterator
    //! pointing at the last letter in the word labelling an edge.
    //!
    //! \complexity
    //! At worst the distance from \p first to \p last.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed, it is
    //! assumed that \p source is a node in the word graph \p wg; and that the
    //! letters in the word described by \p first and \p last belong to the
    //! range \c 0 to WordGraph::out_degree.
    template <typename Node1, typename Node2, typename Iterator>
    [[nodiscard]] std::pair<Node1, Iterator>
    last_node_on_path_no_checks(WordGraph<Node1> const& wg,
                                Node2                   source,
                                Iterator                first,
                                Iterator                last) noexcept;

    //! \brief Returns the last node on the path labelled by a word and an
    //! iterator to the position in the word reached.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //! \tparam Iterator the type of the iterators into a word.
    //!
    //! \param wg a word graph.
    //! \param source the source node.
    //! \param first iterator into a word.
    //! \param last iterator into a word.
    //!
    //! \returns A pair consisting of WordGraph::node_type and \p S.
    //!
    //! \throws LibsemigroupsException if \p source is out of bounds.
    //!
    //! \complexity
    //! At worst the distance from \p first to \p last.
    //!
    //! \note If any value in \p wg or in the word described by \p first and \p
    //! last is out of bounds (greater than or equal to
    //! WordGraph::number_of_nodes), the path labelled by the word exits the
    //! word graph, which is reflected in the result value of this function,
    //! but does not cause an exception to be thrown.
    template <typename Node1, typename Node2, typename Iterator>
    [[nodiscard]] std::pair<Node1, Iterator>
    last_node_on_path(WordGraph<Node1> const& wg,
                      Node2                   source,
                      Iterator                first,
                      Iterator                last);

    //! \brief Returns the last node on the path labelled by a word and an
    //! iterator to the position in the word reached.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg a word graph.
    //! \param source the source node.
    //! \param w the word.
    //!
    //! \returns A pair consisting of the last node reached and an iterator
    //! pointing at the last letter in the word labelling an edge.
    //!
    //! \complexity
    //! At worst the distance from `w.size()`.
    //!
    //! \warning
    //! No checks on the arguments of this function are performed, it is
    //! assumed that \p source is a node in the word graph \p wg; and that the
    //! letters in the word described by \p first and \p last belong to the
    //! range \c 0 to WordGraph::out_degree.
    template <typename Node1, typename Node2>
    std::pair<Node1, word_type::const_iterator>
    last_node_on_path_no_checks(WordGraph<Node1> const& wg,
                                Node2                   source,
                                word_type const&        w);

    //! \brief Returns the last node on the path labelled by a word and an
    //! iterator to the position in the word reached.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg a word graph.
    //! \param source the source node.
    //! \param w the word.
    //!
    //! \returns A pair consisting of the last node reached and an iterator
    //! pointing at the last letter in the word labelling an edge.
    //!
    //! \complexity
    //! At worst the distance from `w.size()`.
    //!
    //! \note If any value in \p wg or in the word described by \p first and \p
    //! last is out of bounds (greater than or equal to
    //! WordGraph::number_of_nodes), the path labelled by the word exits the
    //! word graph, which is reflected in the result value of this function,
    //! but does not cause an exception to be thrown.
    template <typename Node1, typename Node2>
    std::pair<Node1, word_type::const_iterator>
    last_node_on_path(WordGraph<Node1> const& wg,
                      Node2                   source,
                      word_type const&        w);

    //! \brief Returns the std::unordered_set of nodes reachable from a given
    //! node in a word graph.
    //!
    //! This function returns a std::unordered_set consisting of all the nodes
    //! in the word graph \p wg that are reachable from \p source.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg the word graph.
    //! \param source the source node.
    //!
    //! \returns A std::unordered_set consisting of all the nodes in the word
    //! graph \p wg that are reachable from \p source.
    //!
    //! \throws LibsemigroupsException if \p source is out of bounds (greater
    //! than or equal to WordGraph::number_of_nodes).
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
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

    //! \brief Returns the std::unordered_set of nodes reachable from a given
    //! node in a word graph.
    //!
    //! This function returns a std::unordered_set consisting of all the nodes
    //! in the word graph \p wg that are reachable from \p source.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg the word graph.
    //! \param source the source node.
    //!
    //! \returns A std::unordered_set consisting of all the nodes in the word
    //! graph \p wg that are reachable from \p source.
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    //!
    //! \warning The arguments are not checked, and in particular it is assumed
    //! that \p source is a node of \p wg (i.e. less than
    //! WordGraph::number_of_nodes).
    template <typename Node1, typename Node2>
    [[nodiscard]] std::unordered_set<Node1>
    nodes_reachable_from_no_checks(WordGraph<Node1> const& wg, Node2 source);

    //! \brief Returns the number of nodes reachable from a given node in a
    //! word graph.
    //!
    //! This function returns the number of nodes in the word graph \p wg that
    //! are reachable from \p source.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg the word graph.
    //! \param source the source node.
    //!
    //! \returns The number of nodes in the word graph \p wg that are reachable
    //! from \p source.
    //!
    //! \throws LibsemigroupsException if \p source is out of bounds (greater
    //! than or equal to WordGraph::number_of_nodes).
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    template <typename Node1, typename Node2>
    [[nodiscard]] size_t
    number_of_nodes_reachable_from(WordGraph<Node1> const& wg, Node2 source) {
      return nodes_reachable_from(wg, source).size();
    }

    //! \brief Returns the number of nodes reachable from a given node in a
    //! word graph.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg the word graph.
    //! \param source the source node.
    //!
    //! \returns The number of nodes in the word graph \p wg that are reachable
    //! from \p source.
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    //!
    //! \warning The arguments are not checked, and in particular it is assumed
    //! that \p source is a node of \p wg (i.e. less than
    //! WordGraph::number_of_nodes).
    template <typename Node1, typename Node2>
    [[nodiscard]] size_t
    number_of_nodes_reachable_from_no_checks(WordGraph<Node1> const& wg,
                                             Node2                   source) {
      return nodes_reachable_from_no_checks(wg, source).size();
    }

    //! \brief Construct a random connected acyclic word graph with given number
    //! of nodes, and out-degree.
    //!
    //! This function constructs a random acyclic connected word graph with \p
    //! number_of_nodes nodes, and out-degree \p out_degree. This function
    //! implements the Markov chain algorithm given in \cite Carnino2011.
    //!
    //! \param number_of_nodes the number of nodes.
    //! \param out_degree the out-degree of every node.
    //! \param mt a std::mt19937 used as a random source (defaults to:
    //! std::mt19937(std::random_device()())).
    //!
    //! \returns A random connected acyclic word graph.
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p number_of_nodes is less than \c 2
    //! * \p out_degree is less than \c 2
    //!
    //! \par Complexity
    //! The complexity of the implementation is \f$O(n^2)\f$ where \p n is the
    //! number of nodes.
    template <typename Node>
    WordGraph<Node> random_acyclic(size_t       number_of_nodes,
                                   size_t       out_degree,
                                   std::mt19937 mt
                                   = std::mt19937(std::random_device()()));

    //! \brief Replace the contents of a Forest by a spanning tree of the nodes
    //! reachable from a given node in a word graph.
    //!
    //! This function replaces the content of the Forest \p f with a spanning
    //! tree of the nodes reachable from \p root in the word graph \p wg.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p root.
    //!
    //! \param wg the word graph.
    //! \param root the source node.
    //! \param f the Forest object to hold the result.
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    //!
    //! \warning The arguments are not checked, and in particular it is assumed
    //! that \p root is a node of \p wg (i.e. less than
    //! WordGraph::number_of_nodes).
    template <typename Node1, typename Node2>
    void spanning_tree_no_checks(WordGraph<Node1> const& wg,
                                 Node2                   root,
                                 Forest&                 f);

    //! \brief Replace the contents of a Forest by a spanning tree of the nodes
    //! reachable from a given node in a word graph.
    //!
    //! This function replaces the content of the Forest \p f with a spanning
    //! tree of the nodes reachable from \p root in the word graph \p wg.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p root.
    //!
    //! \param wg the word graph.
    //! \param root the source node.
    //! \param f the Forest object to hold the result.
    //!
    //! \throws LibsemigroupsException if \p root is out of bounds, i.e.
    //! greater than or equal to WordGraph::number_of_nodes.
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    template <typename Node1, typename Node2>
    void spanning_tree(WordGraph<Node1> const& wg, Node2 root, Forest& f);

    //! \brief Returns a Forest containing a spanning tree of the nodes
    //! reachable from a given node in a word graph.
    //!
    //! This function returns a Forest containing a spanning
    //! tree of the nodes reachable from \p root in the word graph \p wg.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p root.
    //!
    //! \param wg the word graph.
    //! \param root the source node.
    //!
    //! \returns A Forest object containing a spanning tree.
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    //!
    //! \warning The arguments are not checked, and in particular it is assumed
    //! that \p root is a node of \p wg (i.e. less than
    //! WordGraph::number_of_nodes).
    template <typename Node1, typename Node2>
    [[nodiscard]] Forest spanning_tree_no_checks(WordGraph<Node1> const& wg,
                                                 Node2                   root);

    //! \brief Returns a Forest containing a spanning tree of the nodes
    //! reachable from a given node in a word graph.
    //!
    //! This function returns a Forest containing a spanning
    //! tree of the nodes reachable from \p root in the word graph \p wg.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p root.
    //!
    //! \param wg the word graph.
    //! \param root the source node.
    //!
    //! \returns A Forest object containing a spanning tree.
    //!
    //! \throws LibsemigroupsException if \p root is out of bounds, i.e.
    //! greater than or equal to WordGraph::number_of_nodes.
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    template <typename Node1, typename Node2>
    [[nodiscard]] Forest spanning_tree(WordGraph<Node1> const& wg, Node2 root);

    //! \brief Standardizes a word graph in-place.
    //!
    //! This function standardizes the word graph \p wg according to the
    //! reduction order specified by \p val, and replaces the contents of the
    //! Forest \p f with a spanning tree rooted at \c 0 for the node reachable
    //! from \c 0. The spanning tree corresponds to the order \p val.
    //!
    //! \tparam Graph the type of the word graph \p wg.
    //!
    //! \param wg the word graph.
    //! \param f the Forest object to store the spanning tree.
    //! \param val the order to use for standardization.
    //!
    //! \returns
    //! This function returns \c true if the word graph \p wg is modified by
    //! this function (i.e. it was not standardized already), and \c false
    //! otherwise.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    // Not nodiscard because sometimes we just don't want the output
    template <typename Graph>
    bool standardize(Graph& wg, Forest& f, Order val);

    //! \brief Standardizes a word graph in-place.
    //!
    //! This function standardizes the word graph \p wg according to the
    //! reduction order specified by \p val, and returns a Forest object
    //! containing a spanning tree rooted at \c 0 for the node reachable from
    //! \c 0. The spanning tree corresponds to the order \p val.
    //!
    //! \tparam Graph the type of the word graph \p wg.
    //!
    //! \param wg the word graph.
    //! \param val the order to use for standardization.
    //!
    //! \returns
    //! A std::pair the first entry of which is \c true if the word graph \p wg
    //! is modified by this function (i.e. it was not standardized already), and
    //! \c false otherwise. The second entry is a Forest object containing a
    //! spanning tree for \p wg.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //!
    //! \note If any target of any edge in the word graph \p wg that is out of
    //! bounds, then this is ignored by this function.
    // Not nodiscard because sometimes we just don't want the output
    template <typename Graph>
    std::pair<bool, Forest> standardize(Graph& wg, Order val = Order::shortlex);

    //! \brief Throws if the target of any edge is out of bounds.
    //!
    //! This function throws if any target of any edge in \p wg is out of
    //! bounds (i.e. is greater than or equal to WordGraph::number_of_nodes,
    //! and not equal to \ref UNDEFINED).
    //!
    //! \tparam Node the type of the nodes in \p wg.
    //!
    //! \param wg the word graph to check.
    //!
    //! \throws LibsemigroupsException if any target of any edge in \p wg is
    //! greater than or equal to WordGraph::number_of_nodes and not equal to
    //! \ref UNDEFINED.
    template <typename Node>
    void throw_if_any_target_out_of_bounds(WordGraph<Node> const& wg);

    //! \brief Throws if the target of any edge with source in a given range is
    //! out of bounds.
    //!
    //! This function throws if any target of any edge in \p wg whose source is
    //! in the range defined by \p first and \p last is out of bounds (i.e. is
    //! greater than or equal to WordGraph::number_of_nodes, and not equal to
    //! \ref UNDEFINED).
    //!
    //! \tparam Node the type of the nodes in \p wg.
    //! \tparam Iterator the type of the 2nd and 3rd arguments.
    //!
    //! \param wg the word graph to check.
    //! \param first iterator pointing at the first node to check.
    //! \param last iterator pointing one beyond the last node to check.
    //!
    //! \throws LibsemigroupsException if any target of any edge in \p wg with
    //! source in the range \p first to \p last is greater than or equal to
    //! WordGraph::number_of_nodes and not equal to \ref UNDEFINED.
    //!
    //! \throws LibsemigroupsException if any node in the range \p first to \p
    //! last is out of bounds (i.e. not a node of \p wg).
    template <typename Node, typename Iterator>
    void throw_if_any_target_out_of_bounds(WordGraph<Node> const& wg,
                                           Iterator               first,
                                           Iterator               last);

    //! \brief Throws if a label is out of bounds.
    //!
    //! This function throws if the label \p a is out of bounds, i.e. it is
    //! greater than or equal to `wg.out_degree()`.
    //!
    //! \tparam Node the type of the nodes in \p wg.
    //!
    //! \param wg the word graph.
    //! \param a the label to check.
    //!
    //! \throws LibsemigroupsException if the label \p a is out of bounds.
    // not noexcept because it throws an exception!
    template <typename Node>
    void throw_if_label_out_of_bounds(WordGraph<Node> const&               wg,
                                      typename WordGraph<Node>::label_type a);

    //! \brief Throws if a label is out of bounds.
    //!
    //! This function throws if any of the letters in \p word are out of bounds,
    //! i.e. if they are greater than or equal to `wg.out_degree()`.
    //!
    //! \tparam Node the type of the nodes in \p wg.
    //!
    //! \param wg the word graph.
    //! \param word the word to check.
    //!
    //! \throws LibsemigroupsException if any value in \p word is out of
    //! bounds.
    template <typename Node>
    void throw_if_label_out_of_bounds(WordGraph<Node> const& wg,
                                      word_type const&       word);

    //! \brief Throws if a label is out of bounds.
    //!
    //! This function throws if any of the letters in the word defined by \p
    //! first and \p last is out of bounds, i.e. if they are greater than or
    //! equal to `wg.out_degree()`.
    //!
    //! \tparam Node the type of the nodes in \p wg.
    //! \tparam Iterator the type of the arguments \p first and \p last.
    //!
    //! \param wg the word graph.
    //! \param first iterator pointing at the first letter to check.
    //! \param last iterator pointing one beyond the last letter to check.
    //!
    //! \throws LibsemigroupsException if any value in the word word defined by
    //! \p first and \p last is out of bounds.
    template <typename Node, typename Iterator>
    void throw_if_label_out_of_bounds(WordGraph<Node> const& wg,
                                      Iterator               first,
                                      Iterator               last);

    //! \brief Throws if a node is out of bounds.
    //!
    //! This function throws if the node \p n is out of bounds
    //! i.e. if it is greater than or equal to `wg.number_of_nodes()`.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p n.
    //!
    //! \param wg the word graph.
    //! \param n the node to check.
    //!
    //! \throws LibsemigroupsException if \p n is out of bounds.
    // not noexcept because it throws an exception!
    template <typename Node1, typename Node2>
    void throw_if_node_out_of_bounds(WordGraph<Node1> const& wg, Node2 n);

    //! \brief Throws if any node in a range is out of bounds.
    //!
    //! This function throws if any node in the range from \p first to \p last
    //! is out of bounds i.e. if they are greater than or equal to
    //! `wg.number_of_nodes()`.
    //!
    //! \tparam Node the node type of the word graph.
    //! \tparam Iterator the type of the parameters \p first and \p last.
    //!
    //! \param wg the word graph.
    //! \param first an iterator pointing at the first node to check.
    //! \param last an iterator pointing one beyond the last node to check.
    //!
    //! \throws LibsemigroupsException if any node in the range \p first to \p
    //! last is out of bounds.
    // not noexcept because it throws an exception!
    template <typename Node, typename Iterator1, typename Iterator2>
    void throw_if_node_out_of_bounds(WordGraph<Node> const& wg,
                                     Iterator1              first,
                                     Iterator2              last);

    //! \brief Returns the nodes of the word graph in topological order (see
    //! below) if possible.
    //!
    //! If it is not empty, the returned vector has the property that if an
    //! edge from a node \c n points to a node \c m, then \c m occurs before
    //! \c n in the vector.
    //!
    //! \tparam Node  the type of the nodes of the WordGraph.
    //!
    //! \param wg the word graph.
    //!
    //! \returns
    //! A std::vector of Node types that contains the nodes of
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

    //! Returns the nodes of the word graph reachable from a given node in
    //! topological order (see below) if possible.
    //!
    //! If it is not empty, the returned vector has the property that
    //! if an edge from a node \c n points to a node \c m, then \c m occurs
    //! before \c n in the vector, and the last item in the vector is \p
    //! source.
    //!
    //! \tparam Node1 the node type of the word graph.
    //! \tparam Node2 the type of the node \p source.
    //!
    //! \param wg the WordGraph object to check.
    //! \param source the source node.
    //!
    //! \returns
    //! A std::vector of Node types that contains the nodes reachable from \p
    //! source in \p wg in topological order (if possible) and is otherwise
    //! empty.
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

  }  // namespace word_graph

  namespace detail {

    template <typename T>
    struct IsWordGraphHelper : std::false_type {};

    template <typename Node>
    struct IsWordGraphHelper<WordGraph<Node>> : std::true_type {};
  }  // namespace detail

  //////////////////////////////////////////////////////////////////////////
  // WordGraph - non-member functions
  //////////////////////////////////////////////////////////////////////////

  //! \ingroup word_graph_group
  //! \brief Helper variable template.
  //!
  //! The value of this variable is \c true if the template parameter \p T is
  //! \ref WordGraph for any template parameters.
  //!
  //! \tparam T a type.
  template <typename T>
  static constexpr bool IsWordGraph = detail::IsWordGraphHelper<T>::value;

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
  //! \param os the ostream.
  //! \param wg the word graph.
  //!
  //! \returns
  //! The first parameter \p os.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  template <typename Node>
  std::ostream& operator<<(std::ostream& os, WordGraph<Node> const& wg);

  //! \relates WordGraph
  //!
  //! \brief Constructs a word graph from a number of nodes and targets.
  //!
  //! This function constructs a word graph from its arguments whose
  //! out-degree is specified by the length of the first item
  //! in the 2nd parameter.
  //!
  //! \tparam Return the return type. Must satisfy \ref IsWordGraph<Return>.
  //!
  //! \param num_nodes the number of nodes in the word graph.
  //! \param targets the targets of the word graph.
  //!
  //! \returns A value of type WordGraph.
  //!
  //! \throws LibsemigroupsException
  //! if WordGraph<Node>::target throws when adding edges from \p targets.
  //!
  //! \complexity
  //! \f$O(mn)\f$ where \f$m\f$ is the length of \p targets and \f$n\f$ is the
  //! parameter \p num_nodes.
  //!
  //! \par Example
  //! \code
  //! // Construct a word graph with 5 nodes and 10 edges (7 specified)
  //! make<WordGraph<uint8_t>>(5, {{0, 0}, {1, 1}, {2}, {3, 3}});
  //! \endcode
  // Passing the 2nd parameter "targets" by value disambiguates it from the
  // other make<WordGraph>.
  template <typename Return>
  [[nodiscard]] std::enable_if_t<IsWordGraph<Return>, Return>
  make(size_t                                                         num_nodes,
       std::initializer_list<std::vector<typename Return::node_type>> targets);

  //! \relates WordGraph
  //!
  //! \copydoc make
  template <typename Return>
  [[nodiscard]] std::enable_if_t<IsWordGraph<Return>, Return>
  make(size_t                                                      num_nodes,
       std::vector<std::vector<typename Return::node_type>> const& targets);

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
  }     // namespace detail

  //! \ingroup word_graph_group
  //! \brief Class for taking joins of word graphs.
  //!
  //! This class exists for its call operators which can be used to find the
  //! join of two word graphs with the same WordGraph::out_degree. This class
  //! implements the Hopcroft-Karp algorithm \cite Hop71 for computing a
  //! finite state automata recognising the union of the languages accepted by
  //! two given automata.
  //!
  //! The input word graphs need not be complete, and the root nodes can also
  //! be specified.
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

    //! \brief Replace the contents of a word graph with the join/meet of two
    //! given word graphs with respect to given root vertices.
    //!
    //! This function replaces the contents of the word graph \p xy with the
    //! join/meet of the word graphs \p x and \p y.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param xy the word graph to store the result.
    //!
    //! \param x the first word graph to join/meet.
    //!
    //! \param xnum_nodes_reachable_from_root the number of nodes reachable in
    //! \p x from the node \p xroot (for the circumstance where this number is
    //! known apriori, and does not have to be recomputed).
    //!
    //! \param xroot the node to use as a root in \p x.
    //!
    //! \param y the second word graph to join/meet.
    //!
    //! \param ynum_nodes_reachable_from_root the number of nodes reachable in
    //! \p y from the node \p yroot.
    //!
    //! \param yroot the node to use as a root in \p y.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        size_t                 xnum_nodes_reachable_from_root,
                        Node                   xroot,
                        WordGraph<Node> const& y,
                        size_t                 ynum_nodes_reachable_from_root,
                        Node                   yroot);

    //! \brief Check if the language accepted by one word graph is contained in
    //! that accepted by another word graph.
    //!
    //! This function returns \c true if the language accepted by \p x with
    //! initial node \p xroot and accept state every node, is a subset of the
    //! corresponding language in \p y.
    //!
    //! \tparam Node1 the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \tparam Node2 the type of the nodes to use as roots.
    //!
    //! \param x the word graph whose language we are checking might be a
    //! subset.
    //!
    //! \param xnum_nodes_reachable_from_root the number of nodes reachable in
    //! \p x from the node \p xroot (for the circumstance where this number is
    //! known apriori, and does not have to be recomputed).
    //!
    //! \param xroot the node to use as the initial state in \p x.
    //!
    //! \param y the word graph whose language we are checking might be a
    //! superset.
    //!
    //! \param ynum_nodes_reachable_from_root the number of nodes reachable in
    //! \p y from the node \p yroot.
    //!
    //! \param yroot the node to use as an initial state in \p y.
    //!
    //! \returns Whether or not \p x is a subrelation of \p y.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    // Is x a subrelation of y?
    template <typename Node1, typename Node2>
    [[nodiscard]] bool
    is_subrelation_no_checks(WordGraph<Node1> const& x,
                             size_t xnum_nodes_reachable_from_root,
                             Node2  xroot,
                             WordGraph<Node1> const& y,
                             size_t ynum_nodes_reachable_from_root,
                             Node2  yroot);
#ifdef PARSED_BY_DOXYGEN
    //! \brief Replace the contents of a word graph with the join/meet of two
    //! given word graphs with respect to given root vertices.
    //!
    //! This function replaces the contents of the word graph \p xy with the
    //! join/meet of the word graphs \p x and \p y. This function is the same as
    //! the 7-argument variant but it computes the number of nodes reachable
    //! from \p xroot and \p yroot.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param xy the word graph to store the result.
    //!
    //! \param x the first word graph to join/meet.
    //!
    //! \param xroot the node to use as a root in \p x.
    //!
    //! \param y the second word graph to join/meet.
    //!
    //! \param yroot the node to use as a root in \p y.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        Node                   xroot,
                        WordGraph<Node> const& y,
                        Node                   yroot);

    //! \brief Replace the contents of a word graph with the join/meet of two
    //! given word graphs with respect to given root vertices.
    //!
    //! This function replaces the contents of the word graph \p xy with the
    //! join/meet of the word graphs \p x and \p y. This function is the same as
    //! the 5-argument variant but it uses \c 0 as the root node in both \p x
    //! and \p y.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param xy the word graph to store the result.
    //!
    //! \param x the first word graph to join/meet.
    //!
    //! \param y the second word graph to join/meet.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        WordGraph<Node> const& y);

    //! \brief Returns a word graph containing the join/meet of two given word
    //! graphs.
    //!
    //! This function returns a word graph  containing the join/meet of the word
    //! graphs \p x and \p y. If \c n is the number of arguments, then this
    //! function constructs a word graph to contain the result, forwards this
    //! and the other arguments to the overload of `call_no_checks` with `n +
    //! 1` parameters, then returns the word graph containing the result.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \tparam Args parameter pack for the remaining arguments.
    //!
    //! \param x the first word graph to join/meet.
    //! \param args the remaining arguments.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    template <typename Node, typename... Args>
    [[nodiscard]] auto call_no_checks(WordGraph<Node> const& x, Args&&... args);

    //! \brief Replace the contents of a word graph with the join/meet of two
    //! given word graphs with respect to given root vertices.
    //!
    //! This function replaces the contents of the word graph \p xy with the
    //! join/meet of the word graphs \p x and \p y. This function is the same as
    //! the 5-argument overload of \c call_no_checks but it throws if its
    //! arguments aren't valid.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param xy the word graph to store the result.
    //!
    //! \param x the first word graph to join/meet.
    //!
    //! \param xroot the node to use as a root in \p x.
    //!
    //! \param y the second word graph to join/meet.
    //!
    //! \param yroot the node to use as a root in \p y.
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p xroot isn't a node in \p x;
    //! * \p yroot isn't a node in \p y;
    //! * `x.out_degree() != y.out_degree()`.
    template <typename Node>
    void operator()(WordGraph<Node>&       xy,
                    WordGraph<Node> const& x,
                    Node                   xroot,
                    WordGraph<Node> const& y,
                    Node                   yroot);

    //! \brief Replace the contents of a word graph with the join/meet of two
    //! given word graphs with respect to given root vertices.
    //!
    //! This function replaces the contents of the word graph \p xy with the
    //! join/meet of the word graphs \p x and \p y. This function is the same as
    //! the 3-argument overload of \c call_no_checks but it throws if its
    //! arguments aren't valid.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param xy the word graph to store the result.
    //!
    //! \param x the first word graph to join/meet.
    //!
    //! \param y the second word graph to join/meet.
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p x has no nodes;
    //! * \p y has no nodes;
    //! * `x.out_degree() != y.out_degree()`.
    template <typename Node>
    void operator()(WordGraph<Node>&       xy,
                    WordGraph<Node> const& x,
                    WordGraph<Node> const& y);

    //! \brief Returns a word graph containing the join/meet of two given word
    //! graphs.
    //!
    //! This function is the same as the overload of \ref call_no_checks with
    //! the same signature, the difference being that this function throws if
    //! the arguments are invalid.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \tparam Args parameter pack for the remaining arguments.
    //!
    //! \param x the first word graph to join/meet.
    //! \param args the remaining arguments.
    //!
    //! \throws LibsemigroupsException if the arguments aren't valid. See the
    //! relevant `operator()` for more details.
    template <typename Node, typename... Args>
    [[nodiscard]] auto operator()(WordGraph<Node> const& x, Args&&... args);

    //! \brief Check if the language accepted by one word graph is contained in
    //! that defined by another word graph.
    //!
    //! This function returns \c true if the language accepted by \p x with
    //! initial node \p xroot and accept state every node, is a subset of the
    //! corresponding language in \p y. This version of the function is similar
    //! to the 6-argument overload, except that here we must compute the number
    //! of nodes in \p x and \p y reachable from \p xroot and \p yroot,
    //! respectively.
    //!
    //! \tparam Node1 the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \tparam Node2 the type of the nodes to use as roots.
    //!
    //! \param x the word graph whose language we are checking might be a
    //! subset.
    //!
    //! \param xroot the node to use as the initial state in \p x.
    //!
    //! \param y the word graph whose language we are checking might be a
    //! superset.
    //!
    //! \param yroot the node to use as an initial state in \p y.
    //!
    //! \returns Whether or not \p x is a subrelation of \p y.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    // Is x a subrelation of y?
    template <typename Node1, typename Node2>
    bool is_subrelation_no_checks(WordGraph<Node1> const& x,
                                  Node2                   xroot,
                                  WordGraph<Node1> const& y,
                                  Node2                   yroot);

    //! \brief Check if the language accepted by one word graph is contained in
    //! that defined by another word graph.
    //!
    //! This function returns \c true if the language accepted by \p x with
    //! initial node \p xroot and accept state every node, is a subset of the
    //! corresponding language in \p y. This version of the function is similar
    //! to the 4-argument overload, except that \c 0 is used as the root node
    //! in both \p x and \p y.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param x the word graph whose language we are checking might be a
    //! subset.
    //!
    //! \param y the word graph whose language we are checking might be a
    //! superset.
    //!
    //! \returns Whether or not \p x is a subrelation of \p y.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    template <typename Node>
    bool is_subrelation_no_checks(WordGraph<Node> const& x,
                                  WordGraph<Node> const& y);

    // There's no subrelation with the number of nodes reachable from the
    // roots as arguments (6 args in total) because we'd have to check that
    // they were valid, and the only way to do this is to recompute them.

    //! \brief Check if the language accepted by one word graph is contained in
    //! that defined by another word graph.
    //!
    //! This function returns \c true if the language accepted by \p x with
    //! initial node \p xroot and accept state every node, is a subset of the
    //! corresponding language in \p y. This version of the function is the
    //! same as the 4-argument overload of \c is_subrelation_no_checks, except
    //! that this function throws if its arguments are invalid.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param x the word graph whose language we are checking might be a
    //! subset.
    //!
    //! \param xroot the node to use as the initial state in \p x.
    //!
    //! \param y the word graph whose language we are checking might be a
    //! superset.
    //!
    //! \param yroot the node to use as an initial state in \p y.
    //!
    //! \returns Whether or not \p x is a subrelation of \p y.
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p xroot isn't a node in \p x;
    //! * \p yroot isn't a node in \p y;
    //! * `x.out_degree() != y.out_degree()`.
    template <typename Node1, typename Node2>
    bool is_subrelation(WordGraph<Node1> const& x,
                        Node2                   xroot,
                        WordGraph<Node1> const& y,
                        Node2                   yroot);

    //! \brief Check if the language accepted by one word graph is contained in
    //! that defined by another word graph.
    //!
    //! This function returns \c true if the language accepted by \p x with
    //! initial node \c 0 and accept state every node, is a subset of the
    //! corresponding language in \p y. This version of the function is the
    //! same as the 2-argument overload of \c is_subrelation_no_checks, except
    //! that this function throws if its arguments are invalid.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \param x the word graph whose language we are checking might be a
    //! subset.
    //!
    //! \param y the word graph whose language we are checking might be a
    //! superset.
    //!
    //! \returns Whether or not \p x is a subrelation of \p y.
    //!
    //! \throws LibsemigroupsException if any of the following hold:
    //! * \p x has no nodes;
    //! * \p y has no nodes;
    //! * `x.out_degree() != y.out_degree()`.
    template <typename Node>
    bool is_subrelation(WordGraph<Node> const& x, WordGraph<Node> const& y);

#else
    using detail::JoinerMeeterCommon<Joiner>::call_no_checks;
    using detail::JoinerMeeterCommon<Joiner>::operator();
    using detail::JoinerMeeterCommon<Joiner>::is_subrelation_no_checks;
    using detail::JoinerMeeterCommon<Joiner>::is_subrelation;
#endif
  };  // Joiner

  //! \ingroup word_graph_group
  //!
  //! \brief Class for taking meets of word graphs.
  //!
  //! This class exists for its call operators which can be used to find the
  //! meet of two word graphs with the same WordGraph::out_degree. This class
  //! implements the same algorithm as that used for computing a
  //! finite state automata recognising the intersection of the languages
  //! accepted by two given automata.
  //!
  //! The input word graphs need not be complete, and the root nodes can also
  //! be specified.
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

    //! \copydoc Joiner::call_no_checks(WordGraph<Node>&, WordGraph<Node>
    //! const&, size_t, Node, WordGraph<Node> const&, size_t, Node)
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        size_t                 xnum_nodes_reachable_from_root,
                        Node                   xroot,
                        WordGraph<Node> const& y,
                        size_t                 ynum_nodes_reachable_from_root,
                        Node                   yroot);

    //! \copydoc Joiner::is_subrelation_no_checks(WordGraph<Node1> const&,
    //! size_t, Node2, WordGraph<Node1> const&, size_t, Node2)
    // is x a subrelation of y
    template <typename Node1, typename Node2>
    [[nodiscard]] bool
    is_subrelation_no_checks(WordGraph<Node1> const& x,
                             size_t xnum_nodes_reachable_from_root,
                             Node2  xroot,
                             WordGraph<Node1> const& y,
                             size_t ynum_nodes_reachable_from_root,
                             Node2  yroot);

#ifdef PARSED_BY_DOXYGEN
    //! \copydoc Joiner::call_no_checks(WordGraph<Node>&, WordGraph<Node>
    //! const&, Node, WordGraph<Node> const&, Node)
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        Node                   xroot,
                        WordGraph<Node> const& y,
                        Node                   yroot);

    //! \copydoc Joiner::call_no_checks(WordGraph<Node>&, WordGraph<Node>
    //! const&, WordGraph<Node> const&)
    template <typename Node>
    void call_no_checks(WordGraph<Node>&       xy,
                        WordGraph<Node> const& x,
                        WordGraph<Node> const& y);

    //! \brief Returns a word graph containing the join/meet of two given word
    //! graphs.
    //!
    //! This function returns a word graph containing the join/meet of the word
    //! graphs \p x and \p y. If \c n is the number of arguments, then this
    //! function constructs a word graph to contain the result, forwards this
    //! and the other arguments to the overload of `call_no_checks` with `n +
    //! 1` parameters, then returns the word graph containing the result.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \tparam Args parameter pack for the remaining arguments.
    //!
    //! \param x the first word graph to join/meet.
    //! \param args the remaining parameters.
    //!
    //! \warning
    //! No checks whatsoever on the validity of the arguments are performed.
    template <typename Node, typename... Args>
    [[nodiscard]] auto call_no_checks(WordGraph<Node> const& x, Args&&... args);

    //! \copydoc Joiner::operator()(WordGraph<Node>&, WordGraph<Node>
    //! const&, Node, WordGraph<Node> const&, Node)
    template <typename Node>
    void operator()(WordGraph<Node>&       xy,
                    WordGraph<Node> const& x,
                    Node                   xroot,
                    WordGraph<Node> const& y,
                    Node                   yroot);

    //! \copydoc Joiner::operator()(WordGraph<Node>&, WordGraph<Node>
    //! const&, WordGraph<Node> const&)
    template <typename Node>
    void operator()(WordGraph<Node>&       xy,
                    WordGraph<Node> const& x,
                    WordGraph<Node> const& y);

    //! \brief Returns a word graph containing the join/meet of two given word
    //! graphs.
    //!
    //! This function is the same as the overload of \ref call_no_checks with
    //! the same signature, the difference being that this function throws if
    //! the arguments are invalid.
    //!
    //! \tparam Node the type of the nodes in the word graphs which are
    //! parameters to this function.
    //!
    //! \tparam Args parameter pack for the remaining arguments.
    //!
    //! \param x the first word graph to join/meet.
    //! \param args the remaining arguments.
    //!
    //! \throws LibsemigroupsException if the arguments aren't valid. See the
    //! relevant `operator()` for more details.
    template <typename Node, typename... Args>
    [[nodiscard]] auto operator()(WordGraph<Node> const& x, Args&&... args);

    //! \copydoc Joiner::is_subrelation_no_checks(WordGraph<Node1> const&,
    //! Node2, WordGraph<Node1> const&, Node2)
    template <typename Node1, typename Node2>
    bool is_subrelation_no_checks(WordGraph<Node1> const& x,
                                  Node2                   xroot,
                                  WordGraph<Node1> const& y,
                                  Node2                   yroot);

    //! \copydoc Joiner::is_subrelation_no_checks(WordGraph<Node> const&,
    //! WordGraph<Node> const&)
    template <typename Node>
    bool is_subrelation_no_checks(WordGraph<Node> const& x,
                                  WordGraph<Node> const& y);

    //! \copydoc Joiner::is_subrelation(WordGraph<Node1> const&,
    //! Node2, WordGraph<Node1> const&, Node2)
    template <typename Node1, typename Node2>
    bool is_subrelation(WordGraph<Node1> const& x,
                        Node2                   xroot,
                        WordGraph<Node1> const& y,
                        Node2                   yroot);

    //! \copydoc Joiner::is_subrelation(WordGraph<Node> const&,
    //! WordGraph<Node> const&)
    template <typename Node>
    bool is_subrelation(WordGraph<Node> const& x, WordGraph<Node> const& y);
#else
    using detail::JoinerMeeterCommon<Meeter>::call_no_checks;
    using detail::JoinerMeeterCommon<Meeter>::operator();
    using detail::JoinerMeeterCommon<Meeter>::is_subrelation_no_checks;
    using detail::JoinerMeeterCommon<Meeter>::is_subrelation;
#endif
  };  // class Meeter

  //! \ingroup word_graph_group
  //!
  //! \brief Return a human readable representation of a WordGraph object.
  //!
  //! Return a human readable representation of a WordGraph object.
  //!
  //! \tparam Node the type of the nodes in the underlying WordGraph.
  //!
  //! \param wg the WordGraph object.
  //!
  //! \returns A string containing the representation.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  template <typename Node>
  [[nodiscard]] std::string to_human_readable_repr(WordGraph<Node> const& wg);

  //! \ingroup word_graph_group
  //!
  //! \brief Return a human readable representation of a Meeter object.
  //!
  //! Return a human readable representation of a Meeter object.
  //!
  //! \param meet the Meeter object.
  //!
  //! \returns A string containing the representation.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] static inline std::string
  to_human_readable_repr(Meeter const& meet) {
    (void) meet;
    return "<Meeter of word graphs>";
  }

  //! \ingroup word_graph_group
  //!
  //! \brief Return a human readable representation of a Joiner object.
  //!
  //! Return a human readable representation of a Joiner object.
  //!
  //! \param join the Joiner object.
  //!
  //! \returns A string containing the representation.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] static inline std::string
  to_human_readable_repr(Joiner const& join) {
    (void) join;
    return "<Joiner of word graphs>";
  }

  //! \ingroup word_graph_group
  //!
  //! \brief Return a string that can be used to recreate a word graph.
  //!
  //! This function returns a std::string containing the input required to
  //! construct a copy of the argument \p wg.
  //!
  //! \tparam Node the type of the nodes of \p wg.
  //! \param wg the word graph.
  //! \param prefix a prefix for the returned string (defaults to an empty
  //! string).
  //! \param braces the braces to use in the string (defaults to `"{}"`).
  //! \param suffix a suffix for the returned string (defaults to an empty
  //! string).
  //!
  //! \returns A string containing the input required to recreate \p wg.
  //!
  //! \throws LibsemigroupsException if the argument \p braces is not of length
  //! \c 2.
  template <typename Node>
  [[nodiscard]] std::string to_input_string(WordGraph<Node> const& wg,
                                            std::string const&     prefix = "",
                                            std::string const& braces = "{}",
                                            std::string const& suffix = "");
}  // namespace libsemigroups

#include "word-graph.tpp"

#endif  // LIBSEMIGROUPS_WORD_GRAPH_HPP_
