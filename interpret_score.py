import sys

def parse_shifts(input):
    return [int(s) for s in input.split(":")]

bit_shifts = parse_shifts(sys.argv[1])
x = long(sys.argv[2])

results = []
while len(bit_shifts):
    results.append(x % (1 << bit_shifts[-1]))
    x >>= bit_shifts[-1]
    del bit_shifts[-1]

results.append(x)

print ", ".join(str(x) for x in reversed(results))
