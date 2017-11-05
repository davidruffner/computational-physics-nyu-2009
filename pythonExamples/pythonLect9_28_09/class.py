class MyClass:

    i = 3

    def __init__(self):
        self.i = 0

    def double(self):
        self.i = self.i*2

Thingy = MyClass()
Thingy.double()
