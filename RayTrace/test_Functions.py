import unittest
import Functions

class TestFunctions(unittest.TestCase):

    def test_absorption(self):
        # Test that absorption does not allow zero values for ps, temp
        with self.assertRaises(ValueError):
            Functions.absorption(1,1,1,0)
        with self.assertRaises(ValueError):
            Functions.absorption(0,1,1,1)
        with self.assertRaises(ValueError):
            Functions.absorption(0,1,1,0)
        # Test that absorption does not allow negative temperature
        with self.assertRaises(ValueError):
            Functions.absorption(1,1,1,-1)
        
        #Assert Value of 0 
        self.assertEqual(Functions.absorption(1,0,1,1), 0)
        
        #I need a verifiable value for absorption to assert that the calculation is performed correctly

    def test_cross(self):
        with self.assertRaises(ValueError):
            Functions.cross([1],[1])
        with self.assertRaises(ValueError):
            Functions.cross([1,1],[1,1])

if __name__=='__main__':
    unittest.main()