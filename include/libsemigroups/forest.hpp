//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2019-2024 Finn Smith + James Mitchell
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

#ifndef LIBSEMIGROUPS_FOREST_HPP_
#define LIBSEMIGROUPS_FOREST_HPP_

#include <algorithm>         // for fill
#include <cstddef>           // for size_t
#include <initializer_list>  // for initializer_list
#include <iterator>          // for begin, end
#include <vector>            // for vector, allocator, operator==

#include "constants.hpp"  // for Undefined, Max, UNDEFINED, operator!=
#include "debug.hpp"      // for LIBSEMIGROUPS_ASSERT
#include "exception.hpp"  // for LIBSEMIGROUPS_EXCEPTION
#include "types.hpp"      // for word_type, enable_if_is_same

namespace libsemigroups {
  //! \ingroup word_graph_group
  //!
  //! \brief Class representing a collection of spanning trees of a word graph.
  //!
  //! Defined in `forest.hpp`.
  //!
  //! This class represents the collection of spanning trees of the strongly
  //! connected components of a word graph.
  class Forest {
    std::vector<size_t> _edge_label;
    std::vector<size_t> _parent;

   public:
    //! Alias for the type of nodes in a forest
    using node_type = size_t;

    //! Alias for the type of edge labels in a forest.
    using label_type = size_t;

    //! \brief Constructs a forest with \p n nodes.
    //!
    //! Constructs a forest with \p n nodes, that is initialised so that the
    //! parent() and label() of every node is UNDEFINED.
    //!
    //! \param n the number of nodes, defaults to \c 0.
    explicit Forest(size_t n = 0)
        : _edge_label(n, static_cast<size_t>(UNDEFINED)),
          _parent(n, static_cast<size_t>(UNDEFINED)) {}

    //! \brief Reinitialize an existing Forest object.
    //!
    //! This function reinitializes an existing Forest object so that it is in
    //! the same state as if it had just be constructed as `Forest(n)`.
    //!
    //! \param n the number of nodes, defaults to \c 0.
    //!
    //! \returns A reference to \c this.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    Forest& init(size_t n = 0);

    //! Default copy constructor
    Forest(Forest const&) = default;

    //! Default move constructor
    Forest(Forest&&) = default;

    //! Default copy assignment constructor
    Forest& operator=(Forest const&) = default;

    //! Default move assignment constructor
    Forest& operator=(Forest&&) = default;

    ~Forest();

    //! \brief Compare Forest objects for equality.
    //!
    //! This function returns \c true if and only if the parent and label of
    //! every node in \c this coincides with the parent and label in \p that.
    //!
    //! \param that the Forest for comparison.
    //!
    //! \returns whether or not the Forest objects are equal.
    [[nodiscard]] bool operator==(Forest const& that) const {
      return _parent == that._parent && _edge_label == that._edge_label;
    }

    //! \brief Compare Forest objects for inequality.
    //!
    //! This function returns the negation of `operator==(that)`
    //!
    //! \param that the Forest for comparison.
    //!
    //! \returns whether or not the Forest objects are not equal.
    [[nodiscard]] bool operator!=(Forest const& that) const {
      return !(*this == that);
    }

    //! \brief Add nodes to the Forest.
    //!
    //! This function adds \p n nodes to the forest, but no edges.
    //!
    //! \param n the number of nodes to add.
    //!
    //! \exceptions
    //! \no_libsemigroups_except
    //! \strong_guarantee
    //!
    //! \complexity
    //! At most linear in `number_of_nodes() + n`.
    //!
    //! \iterator_validity
    //! \iterator_invalid
    Forest& add_nodes(size_t n);

    //! \brief Check if there are any nodes in the forest.
    //!
    //! This function returns \c true if there are 0 nodes in the forest, and \c
    //! false otherwise.
    //!
    //! \returns
    //! A value of type `bool`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    [[nodiscard]] bool empty() const noexcept {
      return _parent.empty();
    }

    //! \brief Set the parent and edge label for a node.
    //!
    //! This function defines the parent and the label of the edge from that
    //! parent to the node \p node to be \p parent and \p gen respectively.
    //!
    //! \param node   the node whose parent and label to set.
    //! \param parent the parent node.
    //! \param gen    the label of the edge from \p parent to \p node.
    //!
    //! \throws LibsemigroupsException if \p node or \p parent exceeds
    //! number_of_nodes().
    //!
    //! \complexity
    //! Constant
    //!
    //! \note The value of \p gen can be arbitrary, and so this argument is
    //! never checked.
    //!
    //! \warning No checks are performed on the arguments of this function. In
    //! particular, if \p node or \p parent is greater than or equal to \ref
    //! number_of_nodes, then bad things may happen.
    Forest& set_parent_and_label_no_checks(node_type  node,
                                           node_type  parent,
                                           label_type gen) {
      _parent[node]     = parent;
      _edge_label[node] = gen;
      return *this;
    }

    //! \brief Set the parent and edge label for a node.
    //!
    //! This function sets the parent of \p node to be \p parent, and the
    //! associated edge-label to be \p gen.
    //!
    //! \param node   the node whose parent and label to set.
    //! \param parent the parent node.
    //! \param gen    the label of the edge from \p parent to \p node.
    //!
    //! \returns A reference to `*this`.
    //!
    //! \throws LibsemigroupsException if \p node or \p parent exceeds
    //! number_of_nodes().
    //!
    //! \complexity
    //! Constant
    //!
    //! \note The value of \p gen can be arbitrary, and so this argument is
    //! never checked.
    // not noexcept because std::vector::operator[] isn't.
    Forest& set_parent_and_label(node_type  node,
                                 node_type  parent,
                                 label_type gen) {
      throw_if_node_out_of_bounds(node);
      throw_if_node_out_of_bounds(parent);
      return set_parent_and_label_no_checks(node, parent, gen);
    }

    //! \brief Returns the number of nodes in the forest.
    //!
    //! Returns the number of nodes in the forest.
    //!
    //! \returns
    //! A `size_t`.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant
    [[nodiscard]] size_t number_of_nodes() const noexcept {
      return _parent.size();
    }

    //! \brief Returns the parent of a node.
    //!
    //! This function returns the parent of the node \p i in the forest. The
    //! value \ref UNDEFINED is returned if \p i has no parent, i.e. when \p i
    //! is a root node.
    //!
    //! \param i the node whose parent is sought.
    //!
    //! \returns
    //! A \ref node_type.
    //!
    //! \throws LibsemigroupsException if \p i exceeds \p number_of_nodes().
    //!
    //! \complexity
    //! Constant
    // not noexcept because std::vector::operator[] isn't.
    [[nodiscard]] node_type parent(node_type i) const {
      throw_if_node_out_of_bounds(i);
      return _parent[i];
    }

    //! \brief Returns the parent of a node.
    //!
    //! This function returns the parent of the node \p i in the forest. The
    //! value \ref UNDEFINED is returned if \p i has no parent, i.e. when \p i
    //! is a root node.
    //!
    //! \param i the node whose parent is sought.
    //!
    //! \returns
    //! A \ref node_type.
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning No checks are performed on the arguments of this function.
    [[nodiscard]] node_type parent_no_checks(node_type i) const {
      return _parent[i];
    }

    //! \brief Returns the label of the edge from a node to its parent.
    //!
    //! This function returns the label of the edge from the parent of \p i to
    //! the node \p i. The value \ref UNDEFINED is returned if \p i has no
    //! parent, i.e. when \p i is a root node.
    //!
    //! \param i the node whose label is sought.
    //!
    //! \returns
    //! A \ref label_type.
    //!
    //! \throws LibsemigroupsException if \p i exceeds \p number_of_nodes().
    //!
    //! \complexity
    //! Constant
    // not noexcept because std::vector::operator[] isn't.
    [[nodiscard]] label_type label(node_type i) const {
      throw_if_node_out_of_bounds(i);
      return _edge_label[i];
    }

    //! \brief Returns the label of the edge from a node to its parent.
    //!
    //! This function returns the label of the edge from the parent of \p i to
    //! the node \p i. The value \ref UNDEFINED is returned if \p i has no
    //! parent, i.e. when \p i is a root node.
    //!
    //! \param i the node whose label is sought.
    //!
    //! \returns
    //! A \ref label_type.
    //!
    //! \complexity
    //! Constant
    //!
    //! \warning No checks are performed on the arguments of this function.
    [[nodiscard]] label_type label_no_checks(node_type i) const {
      return _edge_label[i];
    }

    //! \brief Returns a const reference to the vector of parents.
    //!
    //! Returns a const reference to the vector of parents in the Forest. The
    //! value in position \c i of this vector is the parent of node \c i. If the
    //! parent equals UNDEFINED, then node \c i is a root node.
    //!
    //! \returns
    //! A std::vector<node_type>.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] std::vector<node_type> const& parents() const noexcept {
      return _parent;
    }

    //! \brief Returns a const reference to the vector of edge labels.
    //!
    //! Returns a const reference to the vector of edge labels in the Forest.
    //! The value in position \c i of this vector is the label of the edge from
    //! the parent of node \c i to \c i. If the parent equals UNDEFINED, then
    //! node \c i is a root node.
    //!
    //! \returns
    //! A std::vector<label_type>.
    //!
    //! \exceptions
    //! \noexcept
    //!
    //! \complexity
    //! Constant.
    [[nodiscard]] std::vector<label_type> const& labels() const noexcept {
      return _edge_label;
    }

    //! \brief Modifies \p w to contain the labels of the edges on the path
    //! to a root node from \p i.
    //!
    //! This function modifies its first argument \p w in-place to contain the
    //! labels of the edges on the path to a root node from node \p i.
    //!
    //! \param w value to contain the result.
    //! \param i the node.
    //!
    //! \warning No checks are performed on the arguments of this function.
    // TODO(0) to helper
    void path_to_root_no_checks(word_type& w, node_type i) const;

    //! \brief Store the labels of the edges on the path to a root node from \p
    //! i.
    //!
    //! This function writes labels of the edges on the path to a root node from
    //! node \p i to the iterator \p d_first.
    //!
    //! \tparam Iterator The type of the parameter, and the return type.
    //!
    //! \param d_first the output iterator.
    //! \param i the node.
    //!
    //! \returns An \c Iterator pointing one beyond the last letter inserted
    //! into \p d_first.
    //!
    //! \warning No checks are performed on the arguments of this function.
    template <typename Iterator>
    Iterator path_to_root_no_checks(Iterator d_first, node_type i) const {
      auto it = d_first;
      for (; parent_no_checks(i) != UNDEFINED; ++it) {
        *it = label_no_checks(i);
        LIBSEMIGROUPS_ASSERT(i != parent_no_checks(i));
        i = parent_no_checks(i);
      }
      return it;
    }

    //! \brief Returns a word containing the labels of the edges on the path
    //! to a root node from \p i.
    //!
    //! This function returns a word containing the labels of the edges on the
    //! path to a root node from node \p i.
    //!
    //! \param i the node.
    //!
    //! \returns The word labelling the path from a root node to \p i.
    //!
    //! \warning No checks are performed on the arguments of this function.
    // TODO(0) to helper
    [[nodiscard]] word_type path_to_root_no_checks(node_type i) const {
      word_type w;
      path_to_root_no_checks(w, i);
      return w;
    }

    //! \brief Modifies \p w to contain the labels of the edges on the path
    //! to a root node from \p i.
    //!
    //! This function modifies its first argument \p w in-place to contain the
    //! labels of the edges on the path to a root node from node \p i.
    //!
    //! \param w value to contain the result.
    //! \param i the node.
    //!
    //! \throws LibsemigroupsException if \p i is greater than or equal to \ref
    //! number_of_nodes.
    void path_to_root(word_type& w, node_type i) const {
      throw_if_node_out_of_bounds(i);
      path_to_root_no_checks(w, i);
    }

    //! \brief Returns a word containing the labels of the edges on the path
    //! to a root node from \p i.
    //!
    //! This function returns a word containing the labels of the edges on the
    //! path to a root node from node \p i.
    //!
    //! \param i the node.
    //!
    //! \returns The word labelling the path from a root node to \p i.
    //!
    //! \throws LibsemigroupsException if \p i is greater than or equal to \ref
    //! number_of_nodes.
    [[nodiscard]] word_type path_to_root(node_type i) const {
      word_type w;
      path_to_root(w, i);
      return w;
    }

   private:
    void throw_if_node_out_of_bounds(node_type v) const;
  };

  //! \relates Forest
  //!
  //! \brief Construct a Forest from parents and labels.
  //!
  //! This function constructs a Forest from vector of parents and labels.
  //!
  //! \tparam Return the return type. Must satisfy
  //! `std::is_same<Return, Forest>`.
  //!
  //! \param parent the vector of parents of nodes in the Forest.
  //! \param edge_labels the vector of edge labels in the Forest.
  //!
  //! \returns A newly constructed Forest with parents \p parent and edge labels
  //! \p edge_labels.
  //!
  //! \throws LibsemigroupsException if any of the following hold:
  //! * \p parent and \p edge_labels have different sizes;
  //! * \p parent and \p edge_labels do not have the value \ref UNDEFINED in the
  //! same positions (these values indicate where the roots of the trees in the
  //! forest are located and so must coincide).
  //! * Forest::set_parent_and_label throws for `parent[i]` and `edge_labels[i]`
  //! for any value of `i`.
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, Forest>
  make(std::vector<size_t> parent, std::vector<size_t> edge_labels) {
    if (parent.size() != edge_labels.size()) {
      LIBSEMIGROUPS_EXCEPTION(
          "expected the 1st and 2nd arguments (parents and edge labels) to "
          "have equal size equal, found {} != {}",
          parent.size(),
          edge_labels.size());
    }
    size_t const num_nodes = parent.size();
    Forest       result(num_nodes);
    for (size_t i = 0; i < num_nodes; ++i) {
      auto p = *(parent.begin() + i);
      auto l = *(edge_labels.begin() + i);
      if (p != UNDEFINED && l != UNDEFINED) {
        result.set_parent_and_label(i, p, l);
      } else if (!(p == UNDEFINED && l == UNDEFINED)) {
        LIBSEMIGROUPS_EXCEPTION(
            "roots not at the same indices in the 1st and 2nd arguments "
            "(parents and edge labels), expected UNDEFINED at index {} found "
            "{} and {}",
            i,
            p,
            l);
      }
    }
    return result;
  }

  //! \relates Forest
  //!
  //! \brief Construct a Forest from parents and labels.
  //!
  //! This function constructs a Forest from initializer lists of parents and
  //! labels.
  //!
  //! \tparam Return the return type. Must satisfy
  //! `std::is_same<Return, Forest>`.
  //!
  //! \param parent the initializer list of parents of nodes in the Forest.
  //! \param edge_labels the initializer list of edge labels in the Forest.
  //!
  //! \returns A newly constructed Forest with parents \p parent and edge labels
  //! \p edge_labels.
  //!
  //! \throws LibsemigroupsException if any of the following hold:
  //! * \p parent and \p edge_labels have different sizes;
  //! * \p parent and \p edge_labels do not have the value \ref UNDEFINED in the
  //! same positions (these values indicate where the roots of the trees in the
  //! forest are located and so must coincide).
  //! * Forest::set_parent_and_label throws for `parent[i]` and `edge_labels[i]`
  //! for any value of `i`.
  template <typename Return>
  [[nodiscard]] enable_if_is_same<Return, Forest>
  make(std::initializer_list<size_t> parent,
       std::initializer_list<size_t> edge_labels) {
    return make<Forest>(std::vector<size_t>(parent), std::vector(edge_labels));
  }

  //! \relates Forest
  //!
  //! \brief Return a human readable representation of a Forest object.
  //!
  //! Return a human readable representation of a Forest object.
  //!
  //! \param f the Forest.
  //!
  //! \exceptions
  //! \no_libsemigroups_except
  [[nodiscard]] std::string to_human_readable_repr(Forest const& f);
}  // namespace libsemigroups

template <>
struct fmt::formatter<libsemigroups::Forest> : fmt::formatter<std::string> {
  template <typename FormatContext>
  auto format(libsemigroups::Forest const& f, FormatContext& ctx) const {
    return formatter<string_view>::format(
        fmt::format("{{{}, {}}}", f.parents(), f.labels()), ctx);
  }
};

#endif  // LIBSEMIGROUPS_FOREST_HPP_
