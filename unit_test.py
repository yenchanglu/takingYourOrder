import unittest
import _FaceRecog

class test_eigenface(unittest.TestCase):
	def test_result(self):
            accuracy = _FaceRecog.validation()
            print('Accuracy: ' + str(accuracy * 100) + '%')
            self.assertGreaterEqual(accuracy, 0.8)
if __name__ == "__main__":
	unittest.main()
