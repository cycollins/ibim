//
//  main.cpp
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

#include <iostream>
#include <string>
#include "lhash.h"

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

double int_with_min_max(int value, int min, int max)
{
  if (value < min)
    return 0.0;
  
  if (value >= max)
    return 1.0;
  
  return double(value - min) / (max - min);
}

int main(int argc, const char * argv[]) {
  std::function<double(const std::string&)> ciaas(case_insensitive_ascii_alphabetic_string);
  std::function<double(int)> fwmm(std::bind(int_with_min_max, std::placeholders::_1, 0, 1000));
  
  ibim::lhash<int, std::string, std::string, int> test_hash(ciaas, ciaas, fwmm);
  
  test_hash.insert(1, "fred", "mary", 700);
  
  // insert code here...
  std::cout << "Hello, World!\n";
  return 0;
}
