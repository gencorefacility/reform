import unittest

class TestOrderExample(unittest.TestCase):
    def test_z_third(self):
        print("This runs third despite being defined first")
        
    def test_a_first(self):
        print("This runs first despite being defined second")
        
    def test_b_second(self):
        print("This runs second despite being defined third")
        
if __name__ == '__main__':
    unittest.main()
