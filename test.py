import sys
import Levenshtein

def complement(base):
  if base == 'A':
    return 'T'
  elif base == 'C':
    return 'G'
  elif base == 'G':
    return 'C'
  elif base == 'T':
    return 'A'
  else:
    return 'N'

def reverse_complement(string):
  return ''.join(map(complement, string[::-1]))

def process_sam(filename):
  f = open(filename)
  reads = []
  for line in f:
    if line[0] == '@':
      continue
    fields = line.split()
    if 16 & int(fields[1]) > 0:
      reads.append(reverse_complement(fields[9].strip()))
    else:
      reads.append(fields[9].strip())
  f.close()
  return reads

def process_uncompressed(filename):
  f = open(filename)
  reads = []
  for line in f:
    reads.append(line.strip())
  return reads

def compare(sam, uncompressed):
  sam_reads = process_sam(sam)
  reads = process_uncompressed(uncompressed)
  assert(len(sam_reads) == len(reads))
  for i, (x, y) in enumerate(zip(sam_reads, reads)):
    assert(len(x) == len(y))
    value = Levenshtein.editops(x, y)
    dist = Levenshtein.distance(x, y)
    if dist != 0:
      print "Example: " + str(i)
      print dist
      print "Ref: " + x
      print "Att: " + y

compare(sys.argv[1], sys.argv[2])
