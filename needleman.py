import argparse
import numpy as np
import copy
from pprint import pprint


def S(row, col, match, mismatch, seq1, seq2):
  if seq1[col-1] == seq2[row-1]:
    return match
  else: 
    return mismatch
 
def max_n( diag, left, top ):
  if diag > left and diag > top:
    return (diag, ['d'])
  elif left > diag and left > top:
    return (left, ['l'])
  elif top > diag and top > left:
    return (top, ['t'])
  elif diag == left == top:
    return (diag, ['d', 'l', 't'])
  elif diag == left:
    return (diag, ['d', 'l'])
  elif diag == top:
    return (diag, ['d', 't'])
  elif top == left:
    return (top, ['t', 'l'])
  else:
    print("error")
    return None
   
def hij( left_val, top_val, diag_val, row, col, gap, match, mismatch, seq1, seq2 ):
  diag = diag_val['value'] + S(row, col, match, mismatch, seq1, seq2)
  left = left_val['value'] + gap
  top = top_val['value'] + gap
  ret = max_n( diag, left, top )
  return {
    'value': ret[0],
    'back_trace': ret[1]
  }
 
def move( am, alignments, path, row, col, gap ):
  node = am[row][col]
  if node['back_trace'] is None:
    alignments.append( path )
    return None
 
  assert isinstance( node['back_trace'], list )
  assert len( node['back_trace'] ) != 0
 
  for dir in node['back_trace']:
    next_path = copy.deepcopy( path )
    next_path.append( dir )
    if dir == 't':
      next_row = row - 1
      next_col = col
    elif dir == 'l':
      next_row = row
      next_col = col - 1
    elif dir == 'd':
      next_row = row - 1
      next_col = col -1
    move( am, alignments, next_path, next_row, next_col, gap )
  return alignments

def calc_matrix_fun(seq1, seq2, gap, match, mismatch):
  calc_matrix = [ [None for j in range(len(seq1)+1)] for i in range(len(seq2)+1) ]
  calc_matrix[0][0] = {
    'value': 0,
    'back_trace': None
  }
  for i in range( 1, len(calc_matrix[0]) ):
    calc_matrix[0][i] = {
      'value': calc_matrix[0][i-1]['value'] + gap,
      'back_trace': ['l']
    }
  for i in range( 1, len(calc_matrix) ):
    calc_matrix[i][0] = {
      'value': calc_matrix[i-1][0]['value'] + gap,
      'back_trace': ['t']
    }
  
  
  for col in range( 1, len(calc_matrix[0]) ): # iterating on colums
    for row in range( 1, len(calc_matrix) ): # on rows
      left_val = calc_matrix[row][col-1]
      top_val = calc_matrix[row-1][col]
      diag_val = calc_matrix[row-1][col-1]
      calc_matrix[row][col] = hij(left_val, top_val, diag_val, row, col, gap, match, mismatch, seq1, seq2)

  return calc_matrix

def trace_back( alignment_matrix, gap ):
  alignments = []
  move( alignment_matrix, alignments, [], len(alignment_matrix)-1, len(alignment_matrix[0])-1, gap )
  return alignments
 
def parse():
  # Parse arguments
  parser = argparse.ArgumentParser(description='Calculate alignment score and generate sequence alignment')
  parser.add_argument('-a', '--seq1', required=True, help='Provide filename  of a file containing first sequence')
  parser.add_argument('-b', '--seq2', required=True, help='Provide filename  of a file containing second sequence')
  parser.add_argument('-c', '--config', required=True, help='Provide filename  of a file containing configuration')
  
  args = parser.parse_args()
  seq1_file = args.seq1
  seq2_file = args.seq2
  config_file = args.config
  
  # seq1_file = 'seq1.txt'
  # seq2_file = 'seq2.txt'
  # config_file = 'config.txt'
  
  with open( seq1_file, 'r+' ) as f:
    for idx, line in enumerate( f ):
      if line.strip()[0] == '>':
        continue  # ignore comments in FASTA
      seq1 = str(line.strip()) 
    
  with open( seq2_file, 'r+' ) as f:
    for idx, line in enumerate( f ):
      if line.strip()[0] == '>':
        continue  # ignore comments in FASTA
      seq2 =  str(line.strip()) 
    
  print( "seq1: " + seq1 )  
  print( "seq2: " + seq2 )
  
  config = {}
  with open( config_file, 'r+' ) as f:
    for idx, line in enumerate( f ):
      if line.strip()[0] == '#':
        continue  # ignore comments
      sliced = line.split( '=' )
      if len( sliced ) != 2:
        print( "WARNING: Unparsable line %d." % (idx+1) )
      # remove whitespaces and newlines
      key = sliced[0].strip()
      value = sliced[1].strip()
      if len( value ) != 0:
        try:
          value = int( value )
        except ValueError as exc:
          print( "WARNING: Unable to parse line %d (key: %s) as integer; leaving as string." % (idx+1, key) )
      else:
        print( "WARNING: Empty value on line %d for key %s." % (idx+1, key) )
      config[key] = value

  # Parameters of Needleman-Wunsch algorithm
  try:
    gap = config['GP']
    match = config['SAME']
    mismatch = config['DIFF']
  except KeyError as exc:
    print( "Your configuration file is incomplete. Exception text: " + str(exc) )
    exit()
  
  print( "gap: " + str(gap) )
  print( "match: " + str(match) )
  print( "mismatch: " + str(mismatch) )

  return [seq1, seq2, gap, match, mismatch]

def print_console_output( score, alignments ):
  print( "Score: " + str(score) )
  print( "Number of alignments found: " + str(len(alignments)) )
  print( "Alignments:" )
  for alignment in alignments:
    print( alignment )
  return None

def create_results( seq1, seq2, alignments ):
  results = []
  
  for idx, alignment in enumerate( alignments ):
    print( "%d)" % idx )
    s1 = ""
    s2 = ""
    s1_idx = 0
    s2_idx = 0
    for dir in reversed( alignment ):
      if dir == 'l':
        s1 = s1 + seq1[s1_idx]
        s2 = s2 + '_'
        s1_idx = s1_idx + 1
      elif dir == 't':
        s1 = s1 + '_'
        s2 = s2 + seq2[s2_idx]
        s2_idx = s2_idx + 1
      elif dir == 'd':
        s1 = s1 + seq1[s1_idx]
        s2 = s2 + seq2[s2_idx]
        s1_idx = s1_idx + 1
        s2_idx = s2_idx + 1
      else:
        print( 'what?' )
        exit()
  
    print( s1 )
    print( s2 )
    results.append( (s1, s2) )
  return results

def print_to_file(results, score):
  with open( 'results.txt', 'w+' ) as f:
    f.write( "Score = %d\n\n" % score )
    for result in results:
      f.write( result[0] + '\n' )
      f.write( result[1] + '\n\n' )
  return 0

def main():
  input_param = parse()

  seq1 = input_param[0]
  seq2 = input_param[1]
  gap = input_param[2]
  match = input_param[3]
  mismatch = input_param[4]

  calc_matrix = calc_matrix_fun( seq1, seq2, gap, match, mismatch )
  alignments = trace_back( calc_matrix, gap )
  score = calc_matrix[len(calc_matrix)-1][len(calc_matrix[0])-1]['value']
  
  print_console_output(score, alignments)
  
  results = create_results( seq1, seq2, alignments )
  
  print_to_file(results, score)

  return 0
  
if __name__ == '__main__':
    main()

# USAGE:
# python needleman.py -a seq1.txt -b seq2.txt -c config.txt

# #Config.txt
# GP   = -2
# SAME =  5
# DIFF = -5
# MAX_NUMBER_PATHS = 100
# MAX_SEQ_LENGTH = 100
