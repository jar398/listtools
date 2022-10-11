#!/usr/bin/env python3

# Disjoint union, with injections and discrimination

# out takes a value in the coproduct and turns it into a pair (tag,
# val) where tag tells you which side it came from

class Coproduct:
  def __init__(self, 
               in_left_fun, in_right_fun, case_fun,
               swapped=None):
    self.in_left_fun = in_left_fun
    self.in_right_fun = in_right_fun
    self.case_fun = case_fun
    if not swapped:
      def swapped_case_fun(z, when_left, when_right):
        return case_fun(z, when_right, when_left)
      swapped = Coproduct(in_right_fun, in_left_fun,
                          swapped_case_fun,
                          swapped=self)
    self.swapped = swapped
  def in_left(self, x): return (self.in_left_fun)(x)
  def in_right(self, y): return (self.in_right_fun)(y)
  def case(self, z, when_left, when_right):
    return (self.case_fun)(z, when_left, when_right)
  def swap(self): return self.swapped

def test():
  def in_left_fun(x): return ('A', x)
  def in_right_fun(y): return ('B', y)
  def case_fun(z, when_left, when_right):
    (tag, val) = z
    if tag == 'A': return val
    if tag == 'B': return val
    assert False
  AB = Coproduct(in_left_fun, in_right_fun, case_fun)
  print (AB.in_right(198))
  print (AB.case(AB.in_left(7), lambda x:("win", x), lambda x:"lose"))
  print (AB.case(AB.in_right(199), lambda x:"loser", lambda x:("winner", x)))
  BA = AB.swap()
  print (BA.in_right(198))
  print (BA.case(BA.in_left(7), lambda x:("win", x), lambda x:"lose"))
  print (BA.case(BA.in_right(199), lambda x:"loser", lambda x:("winner", x)))

  def in_left_fun(x): return 2*x+1
  def in_right_fun(y): return 2*y
  def case_fun(z, when_left, when_right):
    if z & 1 != 0: return when_left(z//2)
    else: return when_right(z//2)
  odd2 = Coproduct(in_left_fun, in_right_fun, case_fun)
  print (odd2.case(odd2.in_left(7), lambda x:("win", x), lambda x:"lose"))
  print (odd2.case(odd2.in_right(199), lambda x:"loser", lambda x:("winner", x)))
  print (odd2.in_right(198))

if __name__ == '__main__':
  test()
