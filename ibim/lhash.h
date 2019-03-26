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
  
template<typename DatumType, typename ... Keys> class lhash
{
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

private:
  
  typedef lhash<DatumType, Keys ... > linear_hash_type;
  
  int maximum_order;
  int order;
  std::size_t burst_threshold;
  std::size_t magnitude;
  int current_shift;
  std::tuple<std::function<double(Keys)> ... > hash_funcs;
  
  struct datum_record
  {
    DatumType datum;
    std::size_t full_index;
    std::tuple<Keys ... > keys;
    
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

    datum_record(const DatumType &in_datum, std::size_t in_full_index, std::tuple<Keys ... > &&in_keys)
      : datum(in_datum)
      , full_index(in_full_index)
      , keys(std::move(in_keys))
    {
    }
    
    datum_record(DatumType &&in_datum, std::size_t in_full_index, std::tuple<Keys ... > &&in_keys)
      : datum(std::move(in_datum))
      , full_index(in_full_index)
      , keys(std::move(in_keys))
    {
    }
 };
  
  typedef std::list<datum_record> bucket_type;
  template<typename T> using indexed_table_type = std::vector<T>;
  typedef indexed_table_type<bucket_type> table_type;
  table_type table;
  
  struct index_helper
  {
    std::size_t sub_indices[key_count];
    
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
  
  template<std::size_t ... Is> std::size_t create_index(bool full, std::index_sequence<Is ... >, const std::tuple<Keys ...>& key_set)
  {
    int effective_order = full ? maximum_order : order;
    double multiplicand = double(1 << effective_order);
    double interpolants[key_count]{ std::get<Is>(hash_funcs)(std::get<Is>(key_set)) ... };
    index_helper ih(multiplicand, interpolants);
    return ih.bit_reversal(effective_order);
  }
  
  datum_record *find_record(bucket_type &bucket, const std::tuple<Keys ... > &key_set)
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
  
  void insert_common(DatumType &&datum, std::tuple<Keys ...> &&key_set)
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
        insert_common(std::move(datum), std::move(key_set));
        return;
      }
      
      datum_record new_record(std::move(datum), full_index, std::move(key_set));
      bucket.push_back(std::move(new_record));
    }
  }
  
public:
  
  lhash(int in_maximum_order, int initial_order, std::size_t in_burst_threshold, std::function<double(const Keys &)> ... in_hash_funcs)
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
  
  //  so the following is a bit of an experiment here's the idea: the original insert looked something like:
  //
  //  void insert(const DatumType &datum, const Keys & ... in_keys)
  //  {
  //    std::tuple<Keys ... > key_set(in_keys ... );
  //    insert_common(datum, std::move(key_set));
  //  }
  //
  //  in case we want to use move semantics for the datum...
  //
  //  void insert(DatumType &&datum, const Keys & ... in_keys)
  //  {
  //    std::tuple<Keys ... > key_set(in_keys ... );
  //    insert_common(std::move(datum), std::move(key_set));
  //  }
  //
  //  So that's all well and good, but the above signature ensures that the tuple that holds key set will
  //  be initialized using copy semantics. What if a copy operation for a key object is expensive and a
  //  move is less so Not likely because keys are typically pretty primitive types, but as an exercise, I
  //  wanted to know if I could make one or more of the keys use a std::move to specify that the key data
  //  should be stored in the hash entry using move semantics. So the rewrite below uses an arbitrary set
  //  of types for input, U and T... The compiler will deduce the types for each parameter, that way each
  //  one can be individually exposed to std::forward in case it was specified with std::move. Types are
  //  still locked down, because the template instantiation logic will still require a match in count
  //  for each of the input types, and only the specified types in the tuple can be constructed (move or
  //  copy). Another upside to this signature is that each key will be individually constructed with whatever
  //  inputs are supported by their respective constructors. So even if I didn't care about move vs copy
  //  under user control, this is by far a more flexible calling convention. Follow-up: The logic works,
  //  but one of the problems I note (and it seems endemic to modern C++ coding practices) is an enormous
  //  number of what seem like unnecessary moves/copies. There's generally a big difference between the
  //  two, but it's difficult to know which is preferable and for which template types. Copy semantics
  //  for std::string is optimized with reference-counted separate backing stores, which incidentally
  //  means that it's move operation is usually pretty much the same as its copy. On the other hand, I
  //  initially thought that std::tuple might be the same, with separate backing store for the aggregate
  //  object storage (making a move really quick). It isn't. A move involves a sub-move for each constituent
  //  and a copy involves a sub-copy for each. A std::array move is quick (an exchange of backing store)
  //  and a copy is potentially very slow, with a copy for each object in the array. These things are
  //  discoverable, but there is no contract for any of this, which means if they are used in an
  //  implementation, the behavior could change substantially. One of the great beauties of C and older
  //  C++ was that the relationship between source code and object code was very easy to predict. Starting
  //  with C++11 and the normalization of the use of the old STL and meta-programming, we are programming
  //  for the compiler, not the target hardware. A vast portion of the implementation is hidden, even
  //  for operations that normally would be easily inspectable. I'm not a fan. On the other hand, the
  //  power of these later enhancements is undeniable, and this implementation of IBIM is an example.
  //  In terms of functionality per line of code, the later innovaitons in the language are very
  //  impressive. But extensive profiling is requied to understand the implications of any given
  //  implementation choice, and it is necessary after every build because changes in implementation
  //  could vary whildly between run-time versions. I remain ambivalent.
  
  template<typename U, typename ... T> void insert(U datum, T ... in_keys)
  {
    std::tuple<Keys ... > key_set(std::forward<Keys>(in_keys) ... );
    insert_common(std::forward<DatumType>(datum), std::move(key_set));
  }
  
  // the interval type here is meant as an input to a multi-dimensional range query (one of the main reasons for this data structure). It
  // looks like there might be a std library interval object coming in some version of C++ post-20, but it's not here now. There is an
  // interval in each dimension, and each can be open or closed on either end... at least that's the plan right now. The claim of IBIM is
  // that it can resolve the set of buckets containing the records requested in close to O(1) time, or O(n) in the number of dimensions of
  // the key space, which is fixed for any given instance of an IBIM hash table. I haven't decided how to package the results and iterating
  // over the results will obviously be O(n) in the size of the record set returned, but the identification of the buckets containing
  // the data is unrelated to the number of records in the table.
  //
  // To express an interval, we'll use 4-tuple with an lower and upper bound followed by two bools indicating whether either end of the
  // the interval is closed.
  
  template<typename T> using interval = std::tuple<const T &, const T &, bool, bool>;
  using interval_set_type = std::tuple<interval<Keys> ... >;
  using bucket_index_collection_type = std::vector<std::size_t>;
  using key_index_interval = std::pair<std::size_t, std::size_t>;
  
  class iterator;

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
  
  class iterator : public std::iterator<std::forward_iterator_tag, DatumType>
  {
  friend linear_hash_type;
    
  private:
    const linear_hash_type &table;
    interval_set_type filter_intervals;
    std::array<key_index_interval, key_count> pre_reverse_index_intervals;
    index_helper reverse_merge_helper;
    
    std::size_t current_bucket;
    bool current_bucket_needs_filtering;
    typename bucket_type::iterator current_record;
    bool at_beginning;
    bool at_end;
    
    // I obviously could make this into one compound if-statement, but during debugging, it's best to have things proken out
    // It's kind of a strict requirement that all keys should be both less-than comparable and equals-comperable. Maybe my
    // my initial idea about an "interval" having modifyable closed_left and closed_right is too simplistic. There should
    // be template-classes for each possible combination (closed or open on either side and closed or open on both sides).
    // That way, a user can implement only the comparison operators required.
    
    template<typename T> bool passes(const T &key_val, const interval<T> &check_interval)
    {
      const T &lower_bound(std::get<0>(check_interval));
      const T &upper_bound(std::get<1>(check_interval));
      bool closed_left = std::get<2>(check_interval);
      bool closed_right = std::get<3>(check_interval);
      
      if (key_val < lower_bound)
      {
        return false;
      }
      
      if (!closed_left && (key_val == lower_bound))
      {
        return false;
      }
      
      if (key_val > upper_bound)
      {
        return false;
      }
      
      if (!closed_right && (key_val == lower_bound))
      {
        return false;
      }
      
      return true;
    }
    
    template<std::size_t ... Is> bool current_record_passes_filter(std::index_sequence<Is ... >)
    {
      return (passes(std::get<Is>(current_record.keys), std::get<Is>(filter_intervals)) && ... );
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
      , current_bucket_needs_filtering(true)
      , current_record(table.table[current_bucket].begin())
      , at_beginning(true)
      , at_end(false)
    {
    }
    
    iterator(const iterator &) = default;
    iterator(iterator &&) = default;
    
    iterator &operator++()
    {
      if (!at_end)
      {
        bucket_type &bucket(table.table[current_bucket]);
        
        ++current_record;
        if (current_bucket_needs_filtering)
        {
          for( ;current_record != bucket.end(); ++current_record)
          {
            if (current_record_passes_filter())
            {
              break;
            }
          }
        }
        
        if (current_record == bucket.end())
        {
          bool still_interating = false;
          
          for (int i = 0; i < key_count; ++i)
          {
            std::size_t &dimension_sub_index(reverse_merge_helper.sub_indices[i]);
            if (dimension_sub_index == std::get<1>(pre_reverse_index_intervals[i]))
            {
              dimension_sub_index = std::get<0>(pre_reverse_index_intervals[i]);
            }
            else
            {
              still_interating = true;
              ++dimension_sub_index;
              break;
            }
          }
          
          if (still_interating)
          {
            current_bucket = reverse_merge_helper.bit_reversal(table.order);
            current_record = table.table[current_bucket].begin();
          }
          else
          {
            at_end = true;
          }
        }
      }
      
      return *this;
    }
    
    iterator operator++(int)
    {
      iterator retval = *this;
      ++(*this);
      return retval;
    }
    
    bool operator==(iterator other) const
    {
    }
    
    bool operator!=(iterator other) const
    {
      return !(*this == other);
    }
    
    typename iterator::reference operator*() const
    {
      return DatumType();
    }
  };
  
  iterator range_query(const interval<Keys> & ... intervals)
  {
    // create a tuple from the intervals, to avoid marshalling the parameters to helper functions all over again
    interval_set_type interval_set(intervals ... );
    bool all_in[key_count]{ is_void(intervals) ... };
    iterator range_iterator(std::move(generate_iterator(std::index_sequence_for<Keys ... >{}, interval_set, all_in)));
    return range_iterator;
  }
};
  
}

#endif /* lhash_h */
