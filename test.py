import sys
def dist(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def process_sam(filename):
  f = open(filename)
  reads = []
  for line in f:
    if line[0] == '@':
      continue
    reads.append(line.split()[9].strip())
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
  for x, y in zip(sam_reads, reads):
    value = dist(x, y)
    if value != 0:
      print value
      print x
      print y

compare(sys.argv[1], sys.argv[2])
