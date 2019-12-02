import unittest
import needleman as code
class Test(unittest.TestCase):
  def test_S(self):
    #Given
      #input
    row = 1
    col = 1
    match = 5
    mismatch = -5
    seq1 = 'GAAC'
    seq2 = 'CAAGAC'
      #expected result
    expected_result = mismatch
    #When
    result = code.S(row, col, match, mismatch, seq1, seq2)
    #Then
    self.assertEqual(expected_result, result)

  def test_max_n(self):
    #Given
      #input
    diag = -13
    left = -12
    top = -3
      #expected result
    expected_result = (top, ['t'])
    #When
    result = code.max_n( diag, left, top )
    #Then
    self.assertEqual(expected_result, result)

  def test_hij(self):
    #Given
      #input
    left_val = {'back_trace': ['t'], 'value': -6}
    top_val = {'back_trace': ['t', 'l'], 'value': -6}
    diag_val = {'back_trace': ['t'], 'value': -4}
    row = 3
    col = 1
    gap = -2
    match = 5
    mismatch = -5
    seq1 = 'GAAC'
    seq2 = 'CAAAC'
      #expected result
    expected_result = { 'value': -8, 'back_trace': ['t', 'l'] }
    #When
    result = code.hij( left_val, top_val, diag_val, row, col, gap, match, mismatch, seq1, seq2 )
    #Then
    self.assertEqual(expected_result, result)

if __name__ == '__main__':
  unittest.main()

