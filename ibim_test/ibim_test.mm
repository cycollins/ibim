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
#include "lhash.h"

@interface ibim_test : XCTestCase

@end

@implementation ibim_test

- (void)setUp
{
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown
{
    // Put teardown code here. This method is called after the invocation of each test method in the class.
}

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

- (void)testExample
{
  struct double_with_min_max
  {
    double min;
    double range;
    
    double operator ()(double value)
    {
      double bias = value - min;
      
      if (value < 0.0)
        return 0.0;
      
      if (bias >= range)
        return ibim::not_quite_one;
      
      return (ibim::not_quite_one * bias) / range;
    }
    
    double_with_min_max(double min, double max) : min(min), range(max - min) {}
  };
  
  struct test_datum
  {
    std::string content;
    
    test_datum(const std::string &content)
    : content(content)
    {
    }
    
    test_datum(test_datum &&other)
    : content(std::move(other.content))
    {
      std::cout << "Move constructor of " << content << "_datum." << std::endl;
    }
    
    test_datum(const test_datum &other)
    : content(other.content)
    {
      std::cout << "Copy constructor of " << content << "_datum." << std::endl;
    }
    
    test_datum &operator=(test_datum &&) = default;
    test_datum &operator=(test_datum &) = default;
  };
  
  struct key_class
  {
    std::string id;
    
    key_class(const std::string &id)
    : id(id)
    {
    }
    
    key_class(key_class &&other)
    : id(std::move(other.id))
    {
      std::cout << "Move constructor of " << id << "_key." << std::endl;
    }
    
    key_class(const key_class &other)
    : id(other.id)
    {
      std::cout << "Copy constructor of " << id << "_key." << std::endl;
    }
    
    key_class &operator=(key_class &&) = default;
    key_class &operator=(key_class &) = default;
    
    bool operator==(const key_class& other) const
    {
      return (id == other.id);
    }
  };
  
  struct key_class_hash
  {
    double operator()(const key_class &key)
    {
      return case_insensitive_ascii_alphabetic_string(key.id);
    }
  } kk;
  
  std::function<double(const std::string&)> ciaas(case_insensitive_ascii_alphabetic_string);
  std::function<double(int)> iwmm(std::bind(int_with_min_max, std::placeholders::_1, 0, 1000));
  double_with_min_max dwmm(30.0, 50.0);
  
  ibim::lhash<test_datum, key_class, std::string, int, double> test_hash(kk, ciaas, iwmm, dwmm);
  
  key_class mary_key("mary");
  key_class edward_key("edward");
  
  test_datum angus_datum("angus");
  test_datum bridget_datum("bridget");
  
  test_hash.insert(angus_datum, std::move(mary_key), "fred", 700, 40.0);
  test_hash.insert(std::move(bridget_datum), edward_key, "nancy", 500, 31.0);
}

//- (void)testPerformanceExample {
//    // This is an example of a performance test case.
//    [self measureBlock:^{
//        // Put the code you want to measure the time of here.
//    }];
//}

@end
