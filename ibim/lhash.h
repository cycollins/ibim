//
//  lhash.h
//  ibim
//
//  MIT License
//
//  Copyright (c) 2019 David Collins
//
//  Permission is hereby granted, free of charge, to any person obtaining a copy
//  of this software and associated documentation files (the "Software"), to deal
//  in the Software without restriction, including without limitation the rights
//  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//  copies of the Software, and to permit persons to whom the Software is
//  furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included in all
//  copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//
//  A linear ordered hash has the characteristic that the entire keyspace maps into [0,1) using
//  a hash function that preserves h(k1) <= h(k2) => k1 <= k2. To be useful, the space of keys should
//  map as uniformly as possible into the [0,1) range, but generally, this can be imposed through
//  a creative formulation of a hash function to smooth out statistical "clumping" in the key
//  space, as long as the h(k1) <= h(k2) => k1 <= k2 condition is preserved. This can be difficult to
//  do a priori, so linear ordered hashing is not perfectly general-purpose, but when it is applicable,
//  it has the useful characteristic that range queries can be accomplished in nearly constant time. By
//  adding the interpolation-based index maintenance (IBIM), it is possible to perform multi-dimensional
//  hashing and multi-dimensional range queries.
//

#ifndef lhash_h
#define lhash_h

#include <vector>
#include <list>
#include <utility>
#include <cassert>
#include <algorithm>
#include <iterator>

namespace ibim
{
static const double not_quite_one = std::acos(-1.0) / 3.2;
template<typename T> using interval = std::tuple<const T &, const T &, bool, bool>;
template<typename Key> using hash_function = std::function<double(const Key &)>;

template<typename DatumType, typename ... Keys> class lhash
{
  friend class iterator;
  
public:
  static constexpr int kMinimumOrder = 3;
  static constexpr int kDefaultInitialOrder = 5;
  static constexpr std::size_t kDefaultBurstThreshold = 16;
  static constexpr int kDefaultMaximumOrder = 24; // seems like enough, right?
  static constexpr int kMaximumIndexSpaceOrder = 32; // this is the maximum power of two for any index space
  static constexpr std::size_t key_count = sizeof ... (Keys);
  static constexpr int kActualMaxOrder = int((double(kMaximumIndexSpaceOrder) / key_count) + 0.5);
  static constexpr int kActualMinOrder = std::min(kMinimumOrder, kActualMaxOrder);
  static constexpr bool multidimensional = key_count > 1;
  
  class iterator;

private:
  
  typedef lhash<DatumType, Keys ... > linear_hash_type;
  
  int maximum_order;
  int order;
  std::size_t burst_threshold;
  std::size_t magnitude;
  int current_shift;
  std::tuple<hash_function<Keys> ... > hash_funcs;
  
  using key_set_type = std::tuple<Keys && ... >;
  using key_set_storage_type = std::tuple<Keys ... >;

  struct datum_record
  {
    DatumType datum;
    std::size_t full_index;
    key_set_storage_type keys;
    
    // to support the move semantics version of push_back when rehashing
    datum_record& operator=(datum_record &&from)
    {
      full_index = from.full_index;
      datum = std::move(from.datum);
      keys = std::move(from.keys);
    }
    
    datum_record(datum_record &&from)
      : datum(std::move(from.datum))
      , full_index(from.full_index)
      , keys(std::move(from.keys))
    {
    }

   datum_record(DatumType &&in_datum, std::size_t in_full_index, key_set_type &&in_keys)
      : datum(std::move(in_datum))
      , full_index(in_full_index)
      , keys(std::move(in_keys))
    {
    }
 };
  
  using bucket_type = std::list<datum_record>;
  template<typename T> using indexed_table_type = std::vector<T>;
  using table_type = indexed_table_type<bucket_type>;
  table_type table;
  using interval_set_type = std::tuple<interval<Keys> ... >;
  using key_index_interval = std::pair<std::size_t, std::size_t>;

  struct index_helper
  {
    std::size_t sub_indices[key_count];
    
    index_helper()
      : sub_indices{0}
    {
    }
    
    index_helper(double magnitude, double interpolants[key_count])
    {
      if constexpr(multidimensional)
      {
        double *i_ptr = interpolants;
        std::size_t *s_ptr = sub_indices;
        std::ptrdiff_t count = key_count;

        while (--count >= 0)
        {
          *s_ptr++ = magnitude * (*i_ptr++);
        }
      }
      else
      {
        sub_indices[0] = magnitude * interpolants[0];
      }
    }
    
    //  performs the bit-reversal permutation, combinging the sub_index of the multi-dimensional hash keys into a single
    //  index into the master table by interleaving their bits together and simultaneously reversing the interleaved pattern.
    //  The reversal helps to evenly spread the indices around the table by making sure that the lowest-precision bits of
    //  of the pre-reversal value (which vary in a seemingly chaotic way with respect to the input data they are based on)
    //  have the most influence on creating the final index.

    std::size_t bit_reversal(int order)
    {
      std::size_t index = 0;
      std::size_t bit_check_mask = 0x1;

      if constexpr (multidimensional)
      {
        while (--order >= 0)
        {
          int kc = key_count;
          std::size_t *s_ptr = sub_indices;

          while (--kc >= 0)
          {
            index <<= 1;
            std::size_t pattern = *s_ptr++;
            
            if ((pattern & bit_check_mask) != 0)
            {
              index |= 1;
            }
          }

          bit_check_mask <<= 1;
        }
      }
      else
      {
        std::size_t pattern = sub_indices[0];
        
        while (--order >= 0)
        {
          index <<= 1;
          
          if ((pattern & bit_check_mask) != 0)
          {
            index |= 1;
          }
          
          bit_check_mask <<= 1;
        }
      }
      
      return index;
    }
    
    
  };
  
  template<std::size_t ... Is> std::size_t create_index(bool full, std::index_sequence<Is ... >, const key_set_type &key_set)
  {
    int effective_order = full ? maximum_order : order;
    double multiplicand = double(1 << effective_order);
    double interpolants[key_count]{ std::get<Is>(hash_funcs)(std::get<Is>(key_set)) ... };
    index_helper ih(multiplicand, interpolants);
    return ih.bit_reversal(effective_order);
  }
  
  datum_record *find_record(bucket_type &bucket, const key_set_type &key_set)
  {
    for (datum_record &d : bucket)
    {
      if (d.keys == key_set)
      {
        return &d;
      }
    }

    return nullptr;
  }
  
  // one of the cool things about IBIM is that you never have to call the original hash functions after the keys are reduced to their
  // index bit pattern. Rehashing can be accomplished by simply nasking off the bits not needed for the current order of the
  // index space. This makes bursting very efficient (and contracting, though that logic isn't in yet), regardless of how expensive
  // the hash functions may be. Typically, the hash functions are pretty straightforward, though they do often require some floating
  // point multiplies - not a big deal, but more expensive typically than some logical, bit-wise operations. The only down-size to
  // the rehash operation is that it's very non-cache-coherent - lots of random access. Maybe there are refinements that could help,
  // I think this is a problem endemic to any rehashing scheme. The original paper suggested a kind of in-place rehash, starting with
  // what would now look like a table.resize() without the move to temp. In theory, this would in general result in a lot of move
  // operations for the bucket_type objects in the table. I odn't think there are any guarantees about how chache-coherent such an
  // operation would be, so burst_table first uses move semantics to assign the table contents to a temporary before the resize. The
  // move should be O(1) because it only needs to move the backing store from the old table to the temp. Then we move the items as
  // necessary. The original idea from the paper was that for a one-dimensional IBIM table, only half the items would need to be
  // moved, beause if moving the mask revealed a 0 bit, the item would be in the same location in the resized table as it would after
  // the resize. However, this would only be a big win in case where a resize on the heap did not require a block move. All that
  // was also written in the 80's, before any consideration of an array with elements that were C++ objects. In general, a std::vector
  // or similar container resize, is going to involve a lot of moves, so the implementatiopn below simply controls the moves and
  // attempts to do them only once. I believe this is preferable for a modern implementation, but experimentation will tell. Another
  // comsideration arguing against the approach from the paper is that for IBIM tables with dimension N, N > 1, the advantages drop
  // off quickly, such that only about 1/(2^N) items would be expected to remain in place, and IBIM shows most of its advantages for
  // such tables.
  
  void rehash(table_type& old_table, int new_shift)
  {
    for (bucket_type &b : old_table)
    {
      for (datum_record &d : b)
      {
        std::size_t new_index = d.full_index >> new_shift;
        table[new_index].push_back(std::move(d));
      }
    }
  }
  
  void burst_table()
  {
    table_type temp = std::move(table); // should be fast - hopefully just moves the underlying data to temp
    int new_order = order + 1;
    std::size_t new_magnitude = magnitude << 1;
    int new_shift = current_shift - key_count;

    std::size_t new_table_size = temp.size() << key_count;
    table.resize(new_table_size);
    
    rehash(temp, new_shift);
    
    order = new_order;
    magnitude = new_magnitude;
    current_shift = new_shift;
  }
  
  bool size_maintenance(const bucket_type &bucket)
  {
    if ((bucket.size() == burst_threshold) && (order < maximum_order))
    {
      burst_table();
      return true;
    }
    
    return false;
  }
  
  typename bucket_type::iterator find_common(const Keys & ... keys)
  {
  }
  
  void insert_common(DatumType &datum, key_set_type &key_set)
  {
    // the "full_index" here means the final index this set of keys will generate when the table is burst to its
    // highest level. In theory we could just calculate the index for the current burst level, but the one is
    // easily calculated by shifting the full index appropriately to the right. If we insert a new datum_record
    // into a bucket, we store with it the full index, so that it can be repositioned in the table when a burst
    // takes place without calling the hash functions again. This lessens the performance impact of the hash
    // functions in the eent that they are expensive. They generally are not, though they virtually always involve
    // a floating point multiply and divide, which is still more expensive than the shift, not to mention iterating
    // through the tuple of functions and keys and doing the interleave and reversal. The "full" parameter on
    // create_index will be false when we calculate a test index when retrieving data to save a little time, but here
    // in the insert routine, we want both the test index and the full index in case we create a new data record.
    
    std::size_t full_index = create_index(true, std::index_sequence_for<Keys ... >{}, key_set);
    std::size_t index = full_index >> current_shift;
    
    assert(index <= table.size());
    
    bucket_type &bucket = table[index];
    datum_record *found_record = nullptr;
    
    if (!bucket.empty())
    {
      found_record = find_record(bucket, key_set);
    }
    
    if (found_record)
    {
      found_record->datum = std::move(datum);
    }
    else
    {
      if (size_maintenance(bucket))
      {
        insert_common(datum, key_set);
        return;
      }
      
      bucket.emplace_back(std::move(datum), full_index, std::move(key_set));
    }
  }
  
  template<std::size_t ... Is> iterator generate_iterator(std::index_sequence<Is ... >
                                                          , const interval_set_type &interval_set
                                                          , const bool all_in[key_count]
                                                          ) const
  {
    double multiplicand = double(1 << order);
    double lower_interpolants[key_count]{ all_in[Is] ? 0.0 : std::get<Is>(hash_funcs)(std::get<0>(std::get<Is>(interval_set))) ... };
    index_helper l_ih(multiplicand, lower_interpolants);
    double upper_interpolants[key_count]{ all_in[Is] ? 1.0 : std::get<Is>(hash_funcs)(std::get<1>(std::get<Is>(interval_set))) ... };
    index_helper u_ih(multiplicand, upper_interpolants);
    
    std::array<key_index_interval, key_count> orthogonal_pri{ key_index_interval(l_ih.sub_indices[Is], u_ih.sub_indices[Is]) ... };
    
    iterator range_iterator(this, interval_set, orthogonal_pri, l_ih);
    return range_iterator;
  }
  
  const iterator end_iter();

public:
  
  lhash(int in_maximum_order, int initial_order, std::size_t in_burst_threshold, hash_function<Keys> ... in_hash_funcs)
    : maximum_order(std::min(in_maximum_order, kActualMaxOrder))
    , order(std::min<int>(initial_order, maximum_order))
    , burst_threshold(in_burst_threshold)
    , magnitude(1 << order)
    , table(1 << (order * key_count))
    , current_shift((maximum_order - 1) * key_count)
    , hash_funcs(in_hash_funcs ... )
  {
  }
  
  lhash(std::function<double(const Keys &)> ... in_hash_funcs)
    : lhash(kDefaultMaximumOrder
      , kDefaultInitialOrder
      , kDefaultBurstThreshold
      , in_hash_funcs ... )
  {
  }
  
  void insert(DatumType datum, Keys ... in_keys)
  {
    key_set_type key_set(static_cast<Keys &&>(in_keys) ... );
    insert_common(datum, key_set);
  }

  // observation: all the really involved data structures hide their iterator implementations entirely. std::vector and std::array
  // do not because for efficiency, they are aliases of simple pointers. Given how difficult it is to get templated code to work
  // in anything other than a header, I wonder how that works in a case like this, where the parameterized elements of ibim::lhash
  // figure so prominanently in the implementation of its iterator.
  using const_iterator = const iterator;
  class iterator : public std::iterator<std::forward_iterator_tag, DatumType>
  {
  private:
    const linear_hash_type &table;
    interval_set_type filter_intervals;
    std::array<key_index_interval, key_count> pre_reverse_index_intervals;
    index_helper reverse_merge_helper;
    std::size_t current_bucket;
    typename bucket_type::iterator current_record;
    bool ranged;
    bool current_bucket_needs_filtering;
    bool at_beginning;
    bool at_end;
    
    template<typename T> static bool passes(const T &key_val, const interval<T> &check_interval)
    {
      const T &lower_bound(std::get<0>(check_interval));
      const T &upper_bound(std::get<1>(check_interval));
      
      if ((key_val < lower_bound) || (key_val > upper_bound))
      {
        return false;
      }
      
      bool closed_left = std::get<2>(check_interval);
      bool closed_right = std::get<3>(check_interval);
      
      if ((!closed_left && (key_val == lower_bound)) ||
          (!closed_right && (key_val == upper_bound)))
      {
        return false;
      }
      
      return true;
    }
    
    template<std::size_t ... Is> bool current_record_passes_filter(std::index_sequence<Is ... > ) const
    {
      return (passes(std::get<Is>(current_record.keys), std::get<Is>(filter_intervals)) && ... );
    }
    
    void update_needs_filtering(const std::size_t *sub_index_ptr, const key_index_interval *intervals) const
    {
      if (!ranged)
      {
        return;
      }
      
      int i = key_count;
      while (--i >= 0)
      {
        std::size_t dimension_index = *sub_index_ptr++;
        if ((dimension_index == std::get<0>(*intervals)) ||
            (dimension_index == std::get<1>(*intervals)))
        {
          current_bucket_needs_filtering = true;
          return;
        }
        
        ++intervals;
      }
      
      current_bucket_needs_filtering = false;
    }
    
    bool current_record_record_has_valid_data()
    {
      return (!current_bucket_needs_filtering || current_record_passes_filter(std::index_sequence_for<Keys ... >{}));
    }
    
    void advance_one_record()
    {
      ++current_record;
      while (current_record == table.table[current_bucket].end())
      {
        bool still_iterating = false;
        
        for (int i = 0; i < key_count; ++i)
        {
          std::size_t &dimension_sub_index(reverse_merge_helper.sub_indices[i]);
          if (dimension_sub_index == std::get<1>(pre_reverse_index_intervals[i]))
          {
            dimension_sub_index = std::get<0>(pre_reverse_index_intervals[i]);
          }
          else
          {
            still_iterating = true;
            ++dimension_sub_index;
            break;
          }
        }
        
        if (!still_iterating)
        {
          at_end = true;
          return;
        }
        
        update_needs_filtering(reverse_merge_helper.sub_indices, pre_reverse_index_intervals.const_data());
        current_bucket = reverse_merge_helper.bit_reversal(table.order);
        current_record = table.table[current_bucket].begin();
      }
    }

    void ensure_valid_starting_record()
    {
      if ((current_record != table.table[current_bucket].end()) && !current_record_record_has_valid_data())
      {
        ++(*this);
      }
    }
    
  public:
    iterator(const linear_hash_type &in_table
      , const interval_set_type &interval_set
      , const std::array<key_index_interval, key_count> &in_pre_reverse_index_intervals
      , const index_helper &in_reverse_merge_helper
    )
      : table(in_table)
      , filter_intervals(std::move(interval_set))
      , pre_reverse_index_intervals(in_pre_reverse_index_intervals)
      , reverse_merge_helper(in_reverse_merge_helper)
      , current_bucket(reverse_merge_helper.bit_reversal(table.order))
      , current_record(table.table[current_bucket].begin())
      , current_bucket_needs_filtering(true)
      , ranged(true)
      , at_beginning(true)
      , at_end(false)
    {
      ensure_valid_starting_record();
    }
    
    // an iterator used for begin()
    iterator(const linear_hash_type &in_table)
      : table(in_table)
      , pre_reverse_index_intervals{ key_index_interval(0, table.magnitude - 1) } // all buckets in all indices
      , current_bucket(0)
      , current_record(table.table[current_bucket].begin())
      , current_bucket_needs_filtering(false)
      , ranged(false)
      , at_beginning(true)
      , at_end(false)
    {
      ensure_valid_starting_record();
    }
    
    // an iterator constructor for end()
    iterator()
      : current_bucket_needs_filtering(false)
      , ranged(false)
      , at_beginning(false)
      , at_end(true)
    {
    }
    
    iterator(const iterator &) = default;
    iterator(iterator &&) = default;

    iterator &operator++()
    {
      while (!at_end)
      {
        advance_one_record();
        if (current_record_record_has_valid_data())
        {
          break;
        }
      }
      return *this;
    }
    
    inline const_iterator &operator++() const
    {
      return (*this)++;
    }

    iterator operator++(int)
    {
      iterator retval = *this;
      ++(*this);
      return retval;
    }
    
    const_iterator operator++(int) const
    {
      const_iterator retval = *this;
      ++(*this);
      return retval;
    }
    
    bool operator==(const iterator &other) const
    {
      if (&other == &end_iter)
      {
        return at_end;
      }
      
      return (table == other.table) && (**this == *other);
    }
    
    bool operator!=(iterator other) const
    {
      return !(*this == other);
    }
    
    typename iterator::reference operator*()
    {
      return current_record.datum;
    }
    
    inline typename iterator::const_reference operator*() const
    {
      return *(*this);
    }
  };
  
  iterator &end() const
  {
    return end_iter;
  }
  
  const_iterator &cend() const
  {
    return end();
  }
  
  iterator begin()
  {
    return iterator(*this);
  }
  
  const_iterator cbegin() const
  {
    return begin();
  }

  iterator range_query(const interval<Keys> & ... intervals)
  {
    // create a tuple from the intervals, to avoid marshalling the parameters to helper functions all over again
    interval_set_type interval_set(intervals ... );
    bool all_in[key_count]{ is_void(intervals) ... };
    iterator range_iterator(std::move(generate_iterator(std::index_sequence_for<Keys ... >{}, interval_set, all_in)));
    return range_iterator;
  }
  
//  DatumType &operator[](const lookup_type &lookup)
//  {
//    return DatumType();
//  }
};
  
}

#endif /* lhash_h */
