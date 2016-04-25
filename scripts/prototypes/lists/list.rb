# Author: ralf.mueller

class MyVector < Array
  def add(value)
    self << value unless self.include?(value)
    self
  end
end
class MyHash < Hash
  def add(value)
    self[value] = {} unless self.keys.include?(value)
    self
  end
end
