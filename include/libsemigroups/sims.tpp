//
// libsemigroups - C++ library for semigroups and monoids
// Copyright (C) 2024 Reinis Cirpons
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

namespace libsemigroups {

  ////////////////////////////////////////////////////////////////////////
  // SimsSettings - impl of template mem fns
  ////////////////////////////////////////////////////////////////////////

  template <typename Subclass>
  Subclass& SimsSettings<Subclass>::clear_pruners() {
    if (_exclude_pruner_index == UNDEFINED) {
      _pruners.clear();
    } else {
      LIBSEMIGROUPS_ASSERT(_exclude_pruner_index < _pruners.size());
      _pruners.erase(_pruners.cbegin(),
                     _pruners.cbegin() + _exclude_pruner_index);
      _exclude_pruner_index = 0;
      _pruners.erase(_pruners.cbegin() + 1, _pruners.cend());
    }
    return static_cast<Subclass&>(*this);
  }

#ifndef PARSED_BY_DOXYGEN
  template <typename Subclass>
  template <typename OtherSubclass>
  SimsSettings<Subclass>&
  SimsSettings<Subclass>::init(SimsSettings<OtherSubclass> const& that) {
    _exclude      = that.exclude();
    _include      = that.include();
    _presentation = that.presentation();

    _idle_thread_restarts = that.idle_thread_restarts();
    _longs_begin          = _presentation.rules.cbegin()
                   + std::distance(that.presentation().rules.cbegin(),
                                   that.cbegin_long_rules());
    _num_threads = that.number_of_threads();
    _stats       = that.stats();
    _pruners     = that.pruners();
    return *this;
  }
#endif /* PARSED_BY_DOXYGEN */

  template <typename Subclass>
  template <typename Word>
  Subclass& SimsSettings<Subclass>::presentation(Presentation<Word> const& p) {
    if (p.alphabet().empty()) {
      LIBSEMIGROUPS_EXCEPTION(
          "the argument (a presentation) must not have 0 generators");
    }
    // This normalises the rules in the case they are of the right type but
    // not normalised
    Presentation<word_type> p_copy;
    if constexpr (std::is_same_v<Word, word_type>) {
      p_copy = p;
      presentation::normalize_alphabet(p_copy);
    } else {
      p_copy = to_presentation<word_type>(p);
    }
    try {
      presentation::validate_rules(
          p_copy, include().cbegin(), include().cend());
      presentation::validate_rules(
          p_copy, exclude().cbegin(), exclude().cend());
    } catch (LibsemigroupsException const& e) {
      LIBSEMIGROUPS_EXCEPTION(
          "the argument (a presentation) is not compatible with include() and "
          "exclude(), the following exception was thrown:\n{}",
          e.what());
    }
    _presentation = std::move(p_copy);
    _longs_begin  = _presentation.rules.cend();
    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  template <typename Iterator>
  Subclass& SimsSettings<Subclass>::include_exclude(
      Iterator                first,
      Iterator                last,
      std::vector<word_type>& include_or_exclude) {
    if (std::distance(first, last) % 2 != 0) {
      LIBSEMIGROUPS_EXCEPTION("expected the distance between the 1st and 2nd "
                              "arguments (iterators) to be even, found {}",
                              std::distance(first, last));
    }
    for (auto it = first; it != last; ++it) {
      presentation().validate_word(it->cbegin(), it->cend());
    }
    include_or_exclude.assign(first, last);
    return static_cast<Subclass&>(*this);
  }

  template <typename Subclass>
  SimsSettings<Subclass>::~SimsSettings() = default;

  namespace sims {
    template <typename Node>
    bool is_right_congruence(Presentation<word_type> const& p,
                             WordGraph<Node> const&         wg) {
      if (p.alphabet().size() != wg.out_degree()) {
        return false;
      }
      auto const N     = wg.number_of_active_nodes();
      auto       first = wg.cbegin_nodes();
      auto       last  = wg.cbegin_nodes() + N;

      if (!word_graph::is_complete(wg, first, last)
          || !word_graph::is_compatible_no_checks(
              wg, first, last, p.rules.cbegin(), p.rules.cend())) {
        return false;
      }
      auto norf = word_graph::nodes_reachable_from(wg, 0);
      return std::all_of(
          norf.begin(), norf.end(), [&N](auto n) { return n < N; });
    }

    template <typename Node>
    bool is_two_sided_congruence_no_checks(Presentation<word_type> const& p,
                                           WordGraph<Node> const&         wg) {
      auto const N     = wg.number_of_active_nodes();
      auto       first = wg.cbegin_nodes();
      auto       last  = wg.cbegin_nodes() + N;

      for (auto const& rule : right_generating_pairs(p, wg)) {
        if (!word_graph::is_compatible_no_checks(
                wg, first, last, rule.first, rule.second)) {
          return false;
        }
      }
      return true;
    }

    template <typename Node>
    bool is_maximal_right_congruence(Presentation<word_type> const& p,
                                     WordGraph<Node> const&         wg) {
      if (!is_right_congruence(p, wg)) {
        return false;
      } else if (wg.number_of_active_nodes() == 1) {
        // Universal congruence
        return false;
      }

      ToddCoxeter tc;

      auto   tree = word_graph::spanning_tree(wg, 0);
      size_t N    = wg.number_of_active_nodes();

      for (Node x = 0; x < N - 1; ++x) {
        auto wx = tree.path_to_root(x);
        std::reverse(wx.begin(), wx.end());
        for (Node y = x + 1; y < N; ++y) {
          auto wy = tree.path_to_root(y);
          std::reverse(wy.begin(), wy.end());
          auto copy = wg;
          // TODO(2) avoid the copy here
          copy.induced_subgraph_no_checks(static_cast<Node>(0),
                                          wg.number_of_active_nodes());
          tc.init(congruence_kind::onesided, copy)
              .add_generating_pair_no_checks(
                  wx.cbegin(), wx.cend(), wy.cbegin(), wy.cend());
          LIBSEMIGROUPS_ASSERT(tc.word_graph().number_of_nodes()
                               == wg.number_of_active_nodes());
          // fmt::print("x = {}, y = {}\n", x, y);
          // fmt::print("wx = {}, wy = {}\n", wx, wy);
          // std::cout << copy << std::endl;
          // fmt::print("tc.number_of_classes() == {}\n",
          // tc.number_of_classes()); fmt::print("copy.number_of_nodes() ==
          // {}\n", copy.number_of_nodes());
          LIBSEMIGROUPS_ASSERT(tc.number_of_classes()
                               < wg.number_of_active_nodes());
          if (tc.number_of_classes() > 1) {
            return false;
          }
        }
      }
      return true;
    }

    template <typename Iterator>
    BMat<> poset(Iterator first, Iterator last) {
      using WordGraph_ = std::decay_t<decltype(*first)>;
      std::vector<WordGraph_> graphs(first, last);
      size_t const            n = graphs.size();

      Joiner hk;

      BMat<> mat1(n, n);

      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          mat1(i, j) = hk.is_subrelation_no_checks(graphs[i], graphs[j]);
        }
      }

      auto   mat2 = mat1;
      auto   mat3 = mat1;
      BMat<> zero(n, n);
      auto   acc = zero;
      while (mat2 != zero) {
        mat3.product_inplace_no_checks(mat2, mat1);
        std::swap(mat3, mat2);
        acc += mat2;
      }

      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          if (!acc(i, j) && mat1(i, j)) {
            zero(i, j) = true;
          }
        }
      }
      return zero;
    }

    template <typename Iterator>
    Dot dot_poset(Iterator first, Iterator last) {
      auto mat = poset(first, last);
      auto n   = mat.number_of_rows();

      Dot result;

      result.kind(Dot::Kind::digraph)
          .add_attr("node [shape=\"box\"]")
          .add_attr("rankdir=\"BT\"")
          .add_attr("compound=true");
      size_t index = 0;
      for (auto it = first; it != last; ++it) {
        auto copy = *it;
        copy.induced_subgraph_no_checks(0, copy.number_of_active_nodes());
        Dot dot_graph = word_graph::dot(copy);
        dot_graph.name(std::to_string(index++));
        result.add_subgraph(std::move(dot_graph));
      }
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          if (mat(i, j)) {
            result
                .add_edge(fmt::format("cluster_{}_0", i),
                          fmt::format("cluster_{}_0", j))
                .add_attr("minlen", 2.5)
                .add_attr("ltail", fmt::format("cluster_{}", i))
                .add_attr("lhead", fmt::format("cluster_{}", j));
          }
        }
      }
      result.add_attr("splines=line");
      return result;
    }
  }  // namespace sims

}  // namespace libsemigroups
