
# out takes a value in the coproduct and turns it into a pair (tag,
# val) where tag tells you which side it came from

class Coproduct:
  def __init__(self, join_fun, split_fun):
    self.AB = Side(self, join_fun, split_fun)
    def split_BA(z):
      (x, y) = split_fun(z)
      return (y, x)
    self.BA = Side(self,
                   lambda y, x: join_fun(x, y),
                   split_BA)
    self.AB.other_side = self.BA
    self.BA.other_side = self.AB

class Side:
  def __init__(self, coproduct, join_fun, split_fun):
    self.coproduct = coproduct
    self.join_fun = join_fun
    self.split_fun = split_fun
  def flip(self): return self.other_side

  def join(self, x, y): return (self.join_fun)(x, y)
  def split(self, z): return (self.split_fun)(z)
  def in_left(self, x): return self.join(x, None)
  def in_right(self, y): return self.join(None, y)
  def out_left(self, z): return self.split(z)[0]
  def out_right(self, z): return self.split(z)[1]

  def flip(self): return self.other_side

if __name__ == '__main__':
  def join1(x, y):
    if x and y: assert x == y
    if x: return ('A', x)
    else: return ('B', y)
  def split1(z):
    (which, val) = z
    if which == 'A': return (val, None)
    else: return (None, val)
  AB = Coproduct(join1, split1).AB
  print (AB.in_right(198))
  print (AB.out_left(AB.in_left(7)))
  print (AB.out_right(AB.in_right(199)))
  BA = AB.flip()
  print (AB.in_right(198))
  print (BA.out_left(BA.in_left(7)))
  print (BA.out_right(BA.in_right(199)))

  def join2(x, y):
    if x and y: assert x == y
    if x: return 2*x+1
    else: return 2*y
  def split2(z):
    if z & 1 != 0: return (z//2, None)
    else: return (None, z//2)
  odd2 = Coproduct(join2, split2).AB
  print( odd2.out_left(odd2.in_left(33)) )
  print( odd2.out_right(odd2.in_right(44)) )
