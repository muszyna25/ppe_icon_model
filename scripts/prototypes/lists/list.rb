# Author: ralf.mueller

class MyList < Array
  def add(value)
    self << value unless self.include?(value)
  end
end
