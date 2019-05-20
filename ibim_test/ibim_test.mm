//
//  ibim_test.m
//  ibim_test
//
//  Created by David Collins on 3/29/19.
//  Copyright © 2019 David Collins. All rights reserved.
//

#import <XCTest/XCTest.h>
#include <string>
#include <cmath>
#include <iostream>
#include <unordered_set>
#include "lhash.h"

@interface ibim_test : XCTestCase

@end

@implementation ibim_test

double case_insensitive_ascii_alphabetic_string(const std::string &key)
{
  double accumulator = 0.0;
  auto c_ptr_end = key.rend();
  for (auto c_ptr = key.rbegin(); c_ptr != c_ptr_end; c_ptr++)
  {
    int c = *c_ptr;
    if (std::isalpha(c))
    {
      int lc = std::toupper(c);
      
      if ((lc >= 'A') && (lc <= 'Z'))
      {
        double addend = double(lc - 'A');
        accumulator += addend;
        accumulator /= 26.0;
      }
    }
  }
  
  return accumulator;
}

// Integer keys are tricky in IBIM. They have the pathalogical feature that they run the risk of generating nice rational values
// in the required [0 1) real interval. In IBIM, that can make for lots of zero bits at the end of the index when the burst
// table order gets high, and those zero bits will stay as the order gets higher. The bit-reversal trick that we use to make
// the more-chaotic-looking lower-order bits dominate the calculation of final index now betrays us because a rational hash
// result may result in stable lower-order bits for a possibly-large subset of the key-space. This means that the higher
// the table order gets, the more clumping of records into small set of indices (typically, the table gets more smooth and
// homogenous as the order gets higher). So one trick is to take an irational number that is close to 1.0 and multiply all
// results by it so that our results remain in a slightly compressed sub interval of [0 1), which still meets the basic
// requirements for an IBIM H function, but it ensures that every result will be irational and therefor have more varied
// bit-pattern in its index calculation. I apply the same logic to the floating point keys, though it seems less likely
// that floating point data would result in lots of rational values in the interval. On the other hand, I can imagine
// data records that might use double data to represent say GPA's. That would involve lots of integer values and when
// divided by the range [0 4] or [0 5], there could be lots of rational results that might cause less-than-ideal
// behavior.

double int_with_min_max(int value, int min, int max)
{
  if (value < min)
    return 0.0;
  
  if (value >= max)
    return ibim::not_quite_one;
  
  return ibim::not_quite_one * double(value - min) / (max - min);
}

struct insertion_check_info
{
  mutable int count;
  std::string id;
  
  insertion_check_info(const char *in_id)
    : count(1)
    , id(in_id)
  {
  }
  
  bool operator==(const insertion_check_info& other) const
  {
    return id.compare(other.id) == 0;
  }
};

struct insertion_check_info_hash
{
  std::size_t operator()(const insertion_check_info &value) const
  {
    return std::hash<std::string>()(value.id);
  }
};

using check_set = std::unordered_set<insertion_check_info, insertion_check_info_hash>;

void register_entry(check_set &set, const char *entry)
{
  auto found = set.find(entry);
  
  if (found == set.end())
  {
    set.emplace(entry);
  }
  else
  {
    ++(found->count);
  }
}

struct test_datum_class
{
  static check_set move_set;
  static check_set copy_set;
  static check_set move_assign_set;
  static check_set copy_assign_set;

  std::string id;
  
  test_datum_class(const char *in_id)
    : id(in_id)
  {
  }
  
  test_datum_class(test_datum_class &&other)
    : id(std::move(other.id))
  {
    register_entry(move_set, id.c_str());
  }
  
  test_datum_class(const test_datum_class &other)
    : id(other.id)
  {
    register_entry(copy_set, id.c_str());
  }
  
  test_datum_class &operator=(test_datum_class &&other)
  {
    id = std::move(other.id);
    register_entry(move_assign_set, id.c_str());
    return *this;
  }
  
  test_datum_class &operator=(const test_datum_class &other)
  {
    id = other.id;
    register_entry(copy_assign_set, id.c_str());
    return *this;
  }
};

struct test_key_class
{
  static check_set move_set;
  static check_set copy_set;
  
  std::string id;
  
  test_key_class(const char* in_id)
    : id(in_id)
  {
  }
  
  test_key_class(test_key_class &&other)
    : id(std::move(other.id))
  {
    register_entry(move_set, id.c_str());
  }
  
  test_key_class(const test_key_class &other)
    : id(other.id)
  {
    register_entry(copy_set, id.c_str());
  }
  
  bool operator==(const test_key_class &other) const
  {
    return (id.compare(other.id) == 0);
  }
  
  bool operator<(const test_key_class &other) const
  {
    return (id.compare(other.id) < 0);
  }
};

struct double_with_min_max
{
  const double min;
  const double range;
  
  double operator ()(double value) const
  {
    double bias = value - min;
    
    if (value < 0.0)
      return 0.0;
    
    if (bias >= range)
      return ibim::not_quite_one;
    
    return (ibim::not_quite_one * bias) / range;
  }
  
  double_with_min_max(double min, double max)
  : min(min)
  , range(max - min)
  {
  }
};

struct test_key_class_hash
{
  double operator()(const test_key_class &key) const
  {
    return case_insensitive_ascii_alphabetic_string(key.id);
  }
} kch;

check_set test_datum_class::move_set;
check_set test_datum_class::copy_set;
check_set test_datum_class::move_assign_set;
check_set test_datum_class::copy_assign_set;
check_set test_key_class::move_set;
check_set test_key_class::copy_set;

void reset_move_copy_data()
{
  test_datum_class::move_set.clear();
  test_datum_class::copy_set.clear();
  test_datum_class::move_assign_set.clear();
  test_datum_class::copy_assign_set.clear();
  test_key_class::move_set.clear();
  test_key_class::copy_set.clear();
}
- (void)setUp
{
  reset_move_copy_data();
}

- (void)tearDown
{
  // Put teardown code here. This method is called after the invocation of each test method in the class.
}

// maybe too many things to test here, might need breaking down

bool did_datum_move_construct(const char* name, int count)
{
  auto found = test_datum_class::move_set.find(name);
  return (found != test_datum_class::move_set.end()) && (found->count == count);
}

bool did_datum_copy_construct(const char* name, int count)
{
  auto found = test_datum_class::copy_set.find(name);
  return (found != test_datum_class::copy_set.end()) && (found->count == count);
}

bool did_datum_move_assign(const char* name, int count)
{
  auto found = test_datum_class::move_assign_set.find(name);
  return (found != test_datum_class::move_assign_set.end()) && (found->count == count);
}

bool no_datum_move_assign(const char* name)
{
  auto found = test_datum_class::move_assign_set.find(name);
  return (found == test_datum_class::move_assign_set.end());
}

bool no_datum_copy_assign(const char* name)
{
  auto found = test_datum_class::copy_assign_set.find(name);
  return (found == test_datum_class::copy_assign_set.end());
}

bool no_datum_move_construct(const char* name)
{
  auto found = test_datum_class::move_set.find(name);
  return (found == test_datum_class::move_set.end());
}

bool no_datum_copy_construct(const char* name)
{
  auto found = test_datum_class::copy_set.find(name);
  return (found == test_datum_class::copy_set.end());
}

bool did_key_move_construct(const char* name, int count)
{
  auto found = test_key_class::move_set.find(name);
  return (found != test_key_class::move_set.end()) && (found->count == count);
}

bool did_key_copy_construct(const char* name, int count)
{
  auto found = test_key_class::copy_set.find(name);
  return (found != test_key_class::copy_set.end()) && (found->count == count);
}

bool no_key_move_construct(const char* name)
{
  auto found = test_key_class::move_set.find(name);
  return (found == test_key_class::move_set.end());
}

bool no_key_copy_construct(const char* name)
{
  auto found = test_key_class::copy_set.find(name);
  return (found == test_key_class::copy_set.end());
}


- (void)testInsert
{
  ibim::hash_function<std::string> ciaas(case_insensitive_ascii_alphabetic_string);
  ibim::hash_function<int> iwmm(std::bind(int_with_min_max, std::placeholders::_1, 0, 1000));
  
  double_with_min_max dwmm(30.0, 50.0);
  
  {
    reset_move_copy_data();
    ibim::lhash<test_datum_class, test_key_class, std::string, int, double> test_hash(kch, ciaas, iwmm, dwmm);
    
    test_key_class mary_key("mary");
    test_key_class edward_key("edward");
    test_datum_class angus_datum("angus");
    test_datum_class bridget_datum("bridget");
    
    test_hash.insert(angus_datum, std::move(mary_key), "fred", 700, 40.0);
    test_hash.insert(std::move(bridget_datum), edward_key, "nancy", 500, 31.0);
    
    if (!did_key_move_construct("mary", 2)
        || !no_key_copy_construct("mary"))
    {
      XCTAssert(FALSE, @"wrong number of moves for \"mary\"");
    }

    if (!did_key_copy_construct("edward", 1)
        || !did_key_move_construct("edward", 1))
    {
      XCTAssert(FALSE, @"wrong number of moves and copies for \"edward\"");
    }
 
    if (!did_datum_copy_construct("angus", 1)
        || !did_datum_move_construct("angus", 1)
        || !no_datum_move_assign("angus")
        || !no_datum_copy_assign("angus"))
    {
      XCTAssert(FALSE, @"wrong number of moves and copies for \"angus\"");
    }
    
    if (!did_datum_move_construct("bridget", 2)
        || !no_datum_copy_construct("bridget")
        || !no_datum_move_assign("bridget")
        || !no_datum_copy_assign("bridget"))
    {
      XCTAssert(FALSE, @"wrong number of moves for \"bridget\"");
    }

    reset_move_copy_data();
    
    test_datum_class dana_datum("dana");
    test_datum_class fred_datum("fred");
    test_key_class mary_key2("mary");  // old mary_key was moved, and may no longer be valid
    
    test_hash.insert(dana_datum, std::move(mary_key2), "fred", 700, 40.0);
    test_hash.insert(std::move(fred_datum), edward_key, "nancy", 500, 31.0);
    
    if (!did_key_move_construct("mary", 1)
        || !no_key_copy_construct("mary"))
    {
      XCTAssert(FALSE, @"wrong number of moves for overwrite test of \"tina\"");
    }
    
    if (!did_key_copy_construct("edward", 1)
        || !no_key_move_construct("edward"))
    {
      XCTAssert(FALSE, @"wrong number copies for overwrite test of \"edward\"");
    }
    
    if (!did_datum_copy_construct("dana", 1)
        || !did_datum_move_assign("dana", 1)
        || !no_datum_move_construct("dana")
        || !no_datum_copy_assign("dana"))
    {
      XCTAssert(FALSE, @"wrong number of copies and move assigns for overwrite test of \"dana\"");
    }

    if (!did_datum_move_construct("fred", 1)
        || !did_datum_move_assign("fred", 1)
        || !no_datum_copy_construct("fred")
        || !no_datum_copy_assign("fred"))
    {
      XCTAssert(FALSE, @"wrong number of moves and move assigns for overwrite test of \"fred\"");
    }
  }

  {
    reset_move_copy_data();
    
    test_key_class mary_key("mary");
    test_key_class edward_key("edward");
    
    ibim::lhash<int, test_key_class> test_hash(kch);
    
    test_hash.insert(4, std::move(mary_key));
    test_hash.insert(6, edward_key);
    
    if (!did_key_move_construct("mary", 2)
        || !no_key_copy_construct("mary"))
    {
      XCTAssert(FALSE, @"wrong number of moves for 1-dimensional intrinsic data test of \"mary\"");
    }
    
    if (!did_key_copy_construct("edward", 1)
        || !did_key_move_construct("edward", 1))
    {
      XCTAssert(FALSE, @"wrong number of moves and copies for 1-dimensional intrinsic data test of \"edward\"");
    }
    
    reset_move_copy_data();
    
    test_key_class mary_key2("mary");
    
    test_hash.insert(9, std::move(mary_key2));
    test_hash.insert(19, edward_key);
    
    if (!did_key_move_construct("mary", 1)
        || !no_key_copy_construct("mary"))
    {
      XCTAssert(FALSE, @"wrong number of moves for 1-dimensional intrinsic overwrite data test of \"mary\"");
    }
    
    if (!did_key_copy_construct("edward", 1)
        || !no_key_move_construct("edward"))
    {
      XCTAssert(FALSE, @"wrong number of moves and copies for 1-dimensional intrinsic overwrite data test of \"edward\"");
    }
  }
  
  {
    reset_move_copy_data();
    
    test_datum_class angus_datum("angus");
    test_datum_class bridget_datum("bridget");
    
    ibim::lhash<test_datum_class, int> test_hash(iwmm);
    
    test_hash.insert(angus_datum, 700);
    test_hash.insert(std::move(bridget_datum), 500);
    
    if (!did_datum_copy_construct("angus", 1)
        || !did_datum_move_construct("angus", 1)
        || !no_datum_move_assign("angus")
        || !no_datum_copy_assign("angus"))
    {
      XCTAssert(FALSE, @"wrong number of moves and copies for 1-dimensional user-defined datum of \"angus\"");
    }
    
    if (!did_datum_move_construct("bridget", 2)
        || !no_datum_copy_construct("bridget")
        || !no_datum_move_assign("bridget")
        || !no_datum_copy_assign("bridget"))
    {
      XCTAssert(FALSE, @"wrong number of moves for  1-dimensional user-defined datum of \"bridget\"");
    }
    
    reset_move_copy_data();
    
    test_datum_class bridget_datum2("bridget");
    
    test_hash.insert(angus_datum, 700);
    test_hash.insert(std::move(bridget_datum2), 500);
    
    if (!did_datum_copy_construct("angus", 1)
        || !no_datum_copy_assign("angus")
        || !did_datum_move_assign("angus", 1)
        || !no_datum_copy_assign("angus"))
    {
      XCTAssert(FALSE, @"wrong number of moves and copies for 1-dimensional user-defined overwrite datum of \"angus\"");
    }
    
    if (!did_datum_move_construct("bridget", 1)
        || !no_datum_copy_construct("bridget")
        || !did_datum_move_assign("bridget", 1)
        || !no_datum_copy_assign("bridget"))
    {
      XCTAssert(FALSE, @"wrong number of moves for  1-dimensional user-defined overwrite datum of \"bridget\"");
    }
  }
}

//- (void)testPerformanceExample {
//    // This is an example of a performance test case.
//    [self measureBlock:^{
//        // Put the code you want to measure the time of here.
//    }];
//}

@end
