//
// libsemigroups - C++ library for semigroups and monoids
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

// This file contains a generic implementation of a class Action which
// represents the action of a semigroup on a set.

namespace libsemigroups {

  ////////////////////////////////////////////////////////////////////////
  // Action - friends
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  std::ostringstream&
  operator<<(std::ostringstream&                                      os,
             Action<Element, Point, Func, Traits, LeftOrRight> const& action) {
    os << fmt::format("<{} {} action with {} generators, {} points>",
                      action.finished() ? "complete" : "incomplete",
                      LeftOrRight,
                      action.number_of_generators(),
                      action.current_size());
    return os;
  }

  ////////////////////////////////////////////////////////////////////////
  // Action - nested classes - private
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  class Action<Element, Point, Func, Traits, LeftOrRight>::MultiplierCache {
   public:
    [[nodiscard]] element_type& operator[](index_type i) {
      return _multipliers[i].second;
    }

    [[nodiscard]] bool defined(index_type i) const {
      return (i < _multipliers.size() ? _multipliers[i].first : false);
    }

    void set_defined(index_type i) {
      LIBSEMIGROUPS_ASSERT(i < _multipliers.size());
      _multipliers[i].first = true;
    }

    void init(index_type N, element_type const& sample) {
      if (N > _multipliers.size()) {
        _multipliers.resize(N, {false, One()(sample)});
      }
    }

   private:
    std::vector<std::pair<bool, element_type>> _multipliers;
  };

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  struct Action<Element, Point, Func, Traits, LeftOrRight>::Options {
    Options() : _cache_scc_multipliers(false) {}
    Options(Options const&)            = default;
    Options(Options&&)                 = default;
    Options& operator=(Options const&) = default;
    Options& operator=(Options&&)      = default;

    Options& init() {
      _cache_scc_multipliers = false;
      return *this;
    }

    bool _cache_scc_multipliers;
  };

  ////////////////////////////////////////////////////////////////////////
  // Action - constructor + destructor - public
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  Action<Element, Point, Func, Traits, LeftOrRight>::Action()
      : _gens(),
        _graph(),
        _map(),
        _options(),
        _orb(),
        _pos(0),
        _scc(_graph),
        _tmp_point(),
        _tmp_point_init(false) {}

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  Action<Element, Point, Func, Traits, LeftOrRight>&
  Action<Element, Point, Func, Traits, LeftOrRight>::init() {
    _gens.clear();
    _graph.init();
    _map.clear();
    _options.init();
    for (auto pt : _orb) {
      this->internal_free(pt);
    }
    _orb.clear();
    _pos = 0;
    _scc.init(_graph);
    Runner::init();
    // Don't reset _tmp_point or _tmp_point_init
    return *this;
  }

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  Action<Element, Point, Func, Traits, LeftOrRight>::~Action() {
    if (_tmp_point_init) {
      this->internal_free(_tmp_point);
    }
    for (auto pt : _orb) {
      this->internal_free(pt);
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Action - modifiers - public
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  Action<Element, Point, Func, Traits, LeftOrRight>&
  Action<Element, Point, Func, Traits, LeftOrRight>::reserve(size_t val) {
    _graph.reserve(val, _gens.size());
    _map.reserve(val);
    _orb.reserve(val);
    return *this;
  }

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  Action<Element, Point, Func, Traits, LeftOrRight>&
  Action<Element, Point, Func, Traits, LeftOrRight>::add_seed(
      const_reference_point_type seed) {
    auto internal_seed = this->internal_copy(this->to_internal_const(seed));
    if (!_tmp_point_init) {
      _tmp_point_init = true;
      _tmp_point      = this->internal_copy(internal_seed);
    }
    _map.emplace(internal_seed, _orb.size());
    _orb.push_back(internal_seed);
    _graph.add_nodes(1);
    return *this;
  }

  ////////////////////////////////////////////////////////////////////////
  // Action - member functions: position, empty, size, etc - public
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  typename Action<Element, Point, Func, Traits, LeftOrRight>::index_type
  Action<Element, Point, Func, Traits, LeftOrRight>::position(
      const_reference_point_type pt) const {
    auto it = _map.find(this->to_internal_const(pt));
    if (it != _map.end()) {
      return (*it).second;
    } else {
      return UNDEFINED;
    }
  }

  ////////////////////////////////////////////////////////////////////////
  // Runner - pure virtual member functions - private
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  void Action<Element, Point, Func, Traits, LeftOrRight>::run_impl() {
    size_t old_nr_gens = _graph.out_degree();
    _graph.add_to_out_degree(_gens.size() - _graph.out_degree());
    if (started() && old_nr_gens < _gens.size()) {
      // Generators were added after the last call to run
      if (_pos > 0 && old_nr_gens < _gens.size()) {
        _scc.init(_graph);
      }
      for (size_t i = 0; i < _pos; i++) {
        for (size_t j = old_nr_gens; j < _gens.size(); ++j) {
          ActionOp()(this->to_external(_tmp_point),
                     this->to_external_const(_orb[i]),
                     _gens[j]);
          auto it = _map.find(_tmp_point);
          if (it == _map.end()) {
            _graph.add_nodes(1);
            _graph.target(i, j, _orb.size());
            _orb.push_back(this->internal_copy(_tmp_point));
            _map.emplace(_orb.back(), _orb.size() - 1);
          } else {
            _graph.target(i, j, (*it).second);
          }
        }
      }
    }
    if (_pos < _orb.size() && !_gens.empty()) {
      _scc.init(_graph);
    }

    detail::Timer t;
    for (; _pos < _orb.size() && !stopped(); ++_pos) {
      for (size_t j = 0; j < _gens.size(); ++j) {
        ActionOp()(this->to_external(_tmp_point),
                   this->to_external_const(_orb[_pos]),
                   _gens[j]);
        auto it = _map.find(_tmp_point);
        if (it == _map.end()) {
          _graph.add_nodes(1);
          _graph.target(_pos, j, _orb.size());
          _orb.push_back(this->internal_copy(_tmp_point));
          _map.emplace(_orb.back(), _orb.size() - 1);
        } else {
          _graph.target(_pos, j, (*it).second);
        }
      }
      if (report()) {
        report_default("Action: found {:>12} points so far {:>7}\n",
                       detail::group_digits(_orb.size()),
                       t);
      }
    }
    report_why_we_stopped();
  }

  ////////////////////////////////////////////////////////////////////////
  // Action - member functions - private
  ////////////////////////////////////////////////////////////////////////

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  template <bool Forward>
  typename Action<Element, Point, Func, Traits, LeftOrRight>::element_type
  Action<Element, Point, Func, Traits, LeftOrRight>::multiplier_private(
      MultiplierCache& mults,
      Forest const&    f,
      index_type       pos) {
    validate_gens();
    validate_index(pos);

    if (cache_scc_multipliers()) {
      if (mults.defined(pos)) {
        return mults[pos];
      }

      mults.init(_graph.number_of_nodes(), _gens[0]);
      index_type             i = pos;
      std::stack<index_type> visited;
      while (!mults.defined(i) && f.parent(i) != UNDEFINED) {
        visited.push(i);
        mults[i] = _gens[f.label(i)];
        i        = f.parent(i);
      }
      if (visited.empty()) {
        // if pos is the scc root, then this can happen
        mults.set_defined(pos);
      } else {
        element_type tmp = One()(_gens[0]);
        while (i != pos) {
          index_type j = visited.top();
          visited.pop();
          Swap()(tmp, mults[j]);
          if constexpr (Forward) {
            internal_product(mults[j], mults[i], tmp);
          } else {
            internal_product(mults[j], tmp, mults[i]);
          }
          mults.set_defined(j);
          i = j;
        }
      }
      return mults[pos];
    } else {
      element_type out = One()(_gens[0]);
      element_type tmp = One()(_gens[0]);
      while (f.parent(pos) != UNDEFINED) {
        Swap()(tmp, out);
        if constexpr (Forward) {
          internal_product(out, _gens[f.label(pos)], tmp);
        } else {
          internal_product(out, tmp, _gens[f.label(pos)]);
        }
        pos = f.parent(pos);
      }
      return out;
    }
  }

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  void Action<Element, Point, Func, Traits, LeftOrRight>::validate_index(
      index_type i) const {
    if (i > _orb.size()) {
      LIBSEMIGROUPS_EXCEPTION(
          "index out of range, expected value in [0, {}) but found {}",
          current_size(),
          i);
    }
  }

  template <typename Element,
            typename Point,
            typename Func,
            typename Traits,
            side LeftOrRight>
  void
  Action<Element, Point, Func, Traits, LeftOrRight>::validate_gens() const {
    if (_gens.empty()) {
      LIBSEMIGROUPS_EXCEPTION("no generators defined, this function cannot "
                              "be used until at least one generator is added")
    }
  }

}  // namespace libsemigroups
