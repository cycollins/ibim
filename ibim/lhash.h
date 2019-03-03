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
//  Created by David Collins on 2/22/19.
//  Copyright © 2019 David Collins. All rights reserved.
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
#include <cassert>

namespace ibim
{
static const int kDefaultInitialOrder = 3;
  
template<typename DatumType, typename ... Keys> class lhash
{
private:
  
  int order;
  std::size_t magnitude;
  std::tuple<std::function<double(Keys)> ... > hash_funcs;
  struct datum_record
  {
    DatumType datum;
    std::tuple<Keys ...> keys;
    
    datum_record(DatumType in_datum, const Keys & ... in_keys) : datum(in_datum), keys(in_keys ...) {}
    datum_record(DatumType in_datum, const std::tuple<Keys ...>& in_keys) : datum(in_datum), keys(in_keys) {}
  };
  typedef std::list<datum_record> bucket_type;
  typedef std::vector<bucket_type> table_type;
  table_type table;
  
  template<std::size_t ... Is> void calc_interpolants(double interpolants[sizeof ... (Keys)], std::index_sequence<Is ... >, const Keys & ... in_keys)
  {
    double interp[sizeof ... (Keys)] = { std::get<Is>(hash_funcs)(in_keys) ... };
    std::memcpy(interpolants, interp, sizeof(interp)); // how to avoid this copy?
  }
  
  struct index_helper
  {
    std::size_t sub_indices[sizeof ... (Keys)];
    
    index_helper(std::size_t magnitude, double interpolants[sizeof ... (Keys)])
    {
      double *i_ptr = interpolants;
      std::size_t *s_ptr = sub_indices;
      std::ptrdiff_t count = sizeof ... (Keys);

      while (--count >= 0)
      {
        *s_ptr++ = magnitude * (*i_ptr++);
      }
    }
    
    //  performs the bit-reversal permutation, combinging the sub_index of the multi-dimensional hash into a single index into the master table. This will
    //  later allow for efficient, multidimensional range queries.

    std::size_t interleave(int order)
    {
      std::size_t index = 0;
      std::size_t bit_check_mask = 0x1;

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
      
      return index;
    }
  };
  
  std::size_t find_index(const Keys & ... in_keys)
  {
    double interpolants[sizeof ... (Keys)];
    calc_interpolants(interpolants, std::index_sequence_for<Keys ... >{}, in_keys ... );
    index_helper ih(magnitude, interpolants);
    return (sizeof ... (Keys) > 1) ? ih.interleave(order) : ih.sub_indices[0];
  }
  
  datum_record* find_record(std::list<datum_record>& bucket, const std::tuple<Keys ...> &key_set)
  {
    for (auto &d : bucket)
    {
      if (d.keys == key_set)
      {
        return &d;
      }
    }

    return nullptr;
  }
  
public:
  
  lhash(int initial_order, std::function<double(const Keys &)> ... in_hash_funcs)
    : order(initial_order)
    , magnitude(1 << order)
    , table(1 << (order * sizeof ... (Keys)))
    , hash_funcs(in_hash_funcs ... )
  {
  }
  
  lhash(std::function<double(const Keys &)> ... in_hash_funcs) : lhash(kDefaultInitialOrder, in_hash_funcs ...) {}
  
  void insert(DatumType datum, const Keys & ... in_keys)
  {
    std::size_t index = find_index(in_keys ... );
    
    assert(index <= table.size());
    
    bucket_type& bucket = table[index];
    datum_record *found_record = nullptr;
    
    if (!bucket.empty())
    {
      auto key_set = std::make_tuple(in_keys ... );
      found_record = find_record(bucket, key_set);
    }
    
    if (found_record)
    {
      found_record->datum = datum;  // for this simplified example, we'll just overwrite as necessary. We could except or return an error or ignore or ...
    }
    else
    {
      found_record = new datum_record(datum, in_keys ... );
      bucket.push_back(*found_record);
    }
    
    // no burst logic yet, and certainly no contraction logic
  }
};
  
}

#endif /* lhash_h */
