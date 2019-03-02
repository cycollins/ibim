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

namespace ibim
{
static const int kDefaultInitialOrder = 3;
  
template<typename DatumType, typename ... Keys> class lhash
{
private:
  
  int order;
  size_t magnitude;
  size_t key_set_size;
  std::tuple<std::function<double(Keys)> ... > hash_funcs;
  
  struct datum_record
  {
    datum_record(DatumType in_datum, const Keys & ... in_keys) : datum(in_datum), keys(in_keys ...) {}
    DatumType datum;
    
    std::tuple<Keys ...> keys;
  };
  
  template<std::size_t ... Is> void calc_interpolants(double interpolants[], std::index_sequence<Is ... >, const Keys & ... in_keys)
  {
    double interp[sizeof ... (Keys)] = { std::get<Is>(hash_funcs)(in_keys) ... };
    std::memcpy(interpolants, interp, sizeof(interp));
  }
  
  size_t find_index(const Keys & ... in_keys)
  {
    double interpolants[sizeof ... (Keys)];
    calc_interpolants(interpolants, std::index_sequence_for<Keys ... >{}, in_keys ... );
    
    double *i_ptr = interpolants;
    double *end_ptr = interpolants + sizeof ... (Keys);
    
    while (i_ptr != end_ptr)
    {
      
    }
  }
  
  std::vector<std::list<datum_record>> table;
  
public:
  
  lhash(int initial_order, std::function<double(const Keys &)> ... in_hash_funcs)
    : order(initial_order)
    , magnitude(1 << order)
    , key_set_size(sizeof...(Keys))
    , table(magnitude * key_set_size)
    , hash_funcs(in_hash_funcs ... )
  {
  }
  
  lhash(std::function<double(const Keys &)> ... in_hash_funcs) : lhash(kDefaultInitialOrder, in_hash_funcs ...) {}
  
  void insert(DatumType datum, const Keys & ... in_keys)
  {
    std::size_t index = find_index(in_keys ... );
  }
  
};
  
}

#endif /* lhash_h */
