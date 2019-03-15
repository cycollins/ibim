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

namespace ibim
{
template<typename DatumType, typename ... Keys> class lhash
{
private:
  
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
    std::tuple<Keys ...> keys;
    
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

    datum_record(const DatumType &in_datum, std::size_t in_full_index, std::tuple<Keys ...> &&in_keys)
      : datum(in_datum)
      , full_index(in_full_index)
      , keys(std::move(in_keys))
    {
    }
    
    datum_record(DatumType &&in_datum, std::size_t in_full_index, std::tuple<Keys ...> &&in_keys)
      : datum(std::move(in_datum))
      , full_index(in_full_index)
      , keys(std::move(in_keys))
    {
    }
 };
  
  typedef std::list<datum_record> bucket_type;
  typedef std::vector<bucket_type> table_type;
  table_type table;
  
//  template<std::size_t ... Is> void calc_interpolants(double interpolants[sizeof ... (Keys)], std::index_sequence<Is ... >, const std::tuple<Keys ...>& key_set)
//  {
//    interpolants = (double [sizeof ... (Keys)]){ std::get<Is>(hash_funcs)(std::get<Is>(key_set)) ... }; // surprised (but pleased) this works.
//  }
  
  struct index_helper
  {
    std::size_t sub_indices[sizeof ... (Keys)];
    
    index_helper(double magnitude, double interpolants[sizeof ... (Keys)])
    {
      if constexpr(sizeof ... (Keys) > 1)
      {
        double *i_ptr = interpolants;
        std::size_t *s_ptr = sub_indices;
        std::ptrdiff_t count = sizeof ... (Keys);

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

      if constexpr (sizeof ... (Keys) > 1)
      {
        while (--order >= 0)
        {
          int key_count = sizeof ... (Keys);
          std::size_t *s_ptr = sub_indices;

          while (--key_count >= 0)
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
//    double interpolants[sizeof ... (Keys)];
//    calc_interpolants(interpolants, std::index_sequence_for<Keys ... >{}, key_set);
//    std::array<double, sizeof ... (Keys)> interpolants = calc_interpolants2(std::index_sequence_for<Keys ... >{}, key_set);

    double interpolants[sizeof ... (Keys)]{ std::get<Is>(hash_funcs)(std::get<Is>(key_set)) ... };
    index_helper ih(multiplicand, interpolants);
    return ih.bit_reversal(effective_order);
  }
  
  datum_record *find_record(bucket_type &bucket, const std::tuple<Keys ...> &key_set)
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
    int new_shift = current_shift - sizeof ... (Keys);

    std::size_t new_table_size = temp.size() << sizeof ... (Keys);
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
  
  // question: what if copy semantics is slow for the keys? The keys get copied into the tupple and then copied into the datum_record.
  // I'm sure there's a way to accomplish move semantics for the keys if called for. The goal would be to use std::forward so that any
  // individual key could be move or copy. For now, it's OK that that the tuple stores the acutal storage type (kind of unsafe to hash
  // references or pointers anyway.
  
public:
  
//  static const int kMinimumSpread = 3;
  static constexpr int kMinimumOrder = 3;
  static constexpr int kDefaultInitialOrder = 5;
  static constexpr std::size_t kDefaultBurstThreshold = 16;
  static constexpr int kDefaultMaximumOrder = 24; // seems like enough, right?
  static constexpr int kMaximumIndexSpaceOrder = 63; // this is the maximum power of two for any index space (don't worry - it will only be virtual when it gets this large)
  static constexpr int kActualMaxOrder = static_cast<int>(kMaximumIndexSpaceOrder / sizeof ... (Keys));
  static constexpr int kActualMinOrder = std::min(kMinimumOrder, kActualMaxOrder);
  
  lhash(int in_maximum_order, int initial_order, std::size_t in_burst_threshold, std::function<double(const Keys &)> ... in_hash_funcs)
    : maximum_order(std::min(in_maximum_order, kActualMaxOrder))
    , order(std::min<int>(initial_order, maximum_order))
    , burst_threshold(in_burst_threshold)
    , magnitude(1 << order)
    , table(1 << (order * sizeof ... (Keys)))
    , current_shift((maximum_order - 1) * sizeof ... (Keys))
    , hash_funcs(in_hash_funcs ... )
  {
  }
  
  lhash(std::function<double(const Keys &)> ... in_hash_funcs)
    : lhash(kDefaultMaximumOrder
      , kDefaultInitialOrder
      , kDefaultBurstThreshold
      , in_hash_funcs ...)
  {
  }
  
  void insert(const DatumType &datum, const Keys & ... in_keys)
  {
//    static_assert(std::is_integral<T>::value || std::is_same<const DatumType &, T>::value || std::is_same<DatumType &&, T>::value,
//                  "Argument 'datum' can only be integral or of type 'const DatumType &' or 'DatumType &&'");
    
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
    
    std::tuple<Keys ...> key_set(in_keys ...);
    
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
      found_record->datum = datum;
    }
    else
    {
      if (size_maintenance(bucket))
      {
        insert(datum, in_keys ...);
        return;
      }
      
      datum_record new_record(datum, full_index, std::move(key_set));
      bucket.push_back(std::move(new_record));
    }
  }
};
  
}

#endif /* lhash_h */
