import unittest


class test_build(unittest.TestCase):
    def no_op(self): # This is just to have a passing test to validate the CI pipeline
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
